"""
Complete, allele-specific scan of an analysis' gVCF against a local ClinVar release,
and reconciliation of the SQLite catalog with the result.

Why this exists
---------------
The catalog used to be seeded from a filtered MyVariant.info sample and joined to the
gVCF by *position only*. That (a) silently omitted rare pathogenic variants absent from
the seed sample, and (b) could mark a variant "carried" when the sample had a different
allele at the same position. This module replaces that with:

  1. one streaming pass over the whole gVCF; every called ALT allele is matched to
     ClinVar by normalised (chr, pos, deleted, inserted), plus a reference-free haplotype
     equivalence check for indels written differently (genomescope.alleles);
  2. a coverage tally for Pathogenic / Likely pathogenic loci, so "not found" can be told
     apart from "not callable";
  3. reconciliation of the catalog: carried findings are added, and EVERY existing
     snp-indel row is re-adjudicated allele-specifically (false carriers are corrected,
     nothing is deleted); the scan is cross-checked against that independent matcher;
  4. post-condition checks and a JSON / TSV report.

Usage (from backend/; the web pipeline runs the same steps for every new analysis)
-------------------------------------------------------------------------------------
  python -m genomescope.clinvar_scan --workspace JOB_DIR --update
        download + verify + index the newest ClinVar release, scan, update the catalog
  python -m genomescope.clinvar_scan --workspace JOB_DIR
        scan + update using the existing index
  python -m genomescope.clinvar_scan --workspace JOB_DIR --dry-run
        scan + report only; no catalog writes
"""

from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import sys
import time
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from genomescope import clinvar_release as cr
from genomescope import indexer, workspace
from genomescope.alleles import (
    DIFFERENT_ALT_STATUS,
    UNVERIFIED_ALLELE_STATUS,
    is_plain_allele,
    normalize_vcf_allele,
    spans_touch,
    vcf_alleles_equivalent,
)
from genomescope.indexer import get_db_connection, map_review_status_to_evidence, upsert_records
from genomescope.vcf_engine import (
    SNP_VCF_NAME,
    determine_variant_type,
    format_vcf_genotype_and_status,
    normalize_chrom,
    parse_info_field,
    parse_sample_field,
)

KEEP_BACKUPS = 3


def default_report_path() -> Path:
    return workspace.current().scan_report_path


def default_carried_tsv_path() -> Path:
    return workspace.current().scan_carried_path


def backup_dir() -> Path:
    return workspace.current().backup_dir

# A reference block only counts as "confirmed reference" at this depth (gVCF blocks in
# this pipeline report GQ 0, so depth is the usable signal).
MIN_CONFIRM_DP = 10
# Indel representation search window around a called indel (bp), and the longest indel
# (REF or ALT, bp) compared by haplotype equivalence. Longer ClinVar deletions / insertions
# are structural-variant territory (Manta / Canvas files), not small-variant representation
# differences, and a small call that merely lies inside one is a different event.
INDEL_WINDOW = 60
MAX_EQUIV_INDEL = 50

LOW_QUALITY_STATUS = "Low Quality Call (Unconfirmed)"
NO_RECORD_STATUS = "No Call (Low Coverage)"
SCAN_ROW_PREFIX = "CLINVAR:"
STALE_MARKER = "no germline classification for this single variant"
CARRIED_MATCHES = ("exact", "equivalent")

# coverage state per high-impact locus
_COV_NONE, _COV_REF_OK, _COV_REF_LOW, _COV_CALLED = 0, 1, 2, 3


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

@dataclass
class Finding:
    variation_id: int
    chr: str
    pos: int                 # gVCF POS of the row that carries the allele
    ref: str
    alt: str                 # the specific ALT allele matched (not the whole ALT list)
    alt_index: int
    copies: int
    ploidy: int
    gt: str
    filter: str
    match: str               # "exact" | "equivalent"
    tier: int
    stars: int
    cv_pos: int
    cv_ref: str
    cv_alt: str
    line: str                # raw gVCF line (no newline)

    @property
    def passes(self) -> bool:
        return self.filter == "PASS"


@dataclass
class ScanResult:
    findings: List[Finding] = field(default_factory=list)
    review: List[Dict[str, Any]] = field(default_factory=list)   # unresolved high-impact indel candidates
    stats: Dict[str, Any] = field(default_factory=dict)
    coverage: Dict[str, Any] = field(default_factory=dict)


def _parse_gt(gt: str) -> List[Optional[int]]:
    return [int(x) if x.isdigit() else None for x in gt.replace("|", "/").split("/")]


def _format_value(fmt: List[bytes], sample: List[bytes], key: bytes) -> Optional[bytes]:
    try:
        i = fmt.index(key)
    except ValueError:
        return None
    return sample[i] if i < len(sample) else None


def scan_gvcf(
    vcf_path: Path,
    index: cr.ClinVarIndex,
    progress: Optional[Callable[[str], None]] = None,
) -> ScanResult:
    """Stream the whole gVCF once and match every called ALT allele against ClinVar."""
    say = progress or (lambda _m: None)
    t0 = time.time()
    loci = index.high_impact_loci()
    cov: Dict[str, bytearray] = {c: bytearray(len(p)) for c, p in loci.items()}
    chrom_cache: Dict[bytes, str] = {}

    result = ScanResult()
    seen: Dict[Tuple[int, int, str, int], Finding] = {}
    lines = ref_blocks = variant_rows = alleles_checked = nonplain = 0
    benign_pass = benign_filtered = 0
    unresolved_seen: set = set()

    with open(vcf_path, "rb") as f:
        for raw in f:
            if raw[:1] == b"#":
                continue
            cols = raw.rstrip(b"\r\n").split(b"\t")
            if len(cols) < 10:
                continue
            lines += 1
            if lines % 2_000_000 == 0:
                say(f"  scanned {lines:,} gVCF lines ({len(seen):,} ClinVar matches so far) ...")
            cname = cols[0]
            chrom = chrom_cache.get(cname)
            if chrom is None:
                chrom = chrom_cache[cname] = normalize_chrom(cname.decode("ascii", "replace"))
            pos = int(cols[1])
            alt_b = cols[4]
            filt = cols[6].decode("ascii", "replace")
            arr = loci.get(chrom)

            if alt_b == b".":
                ref_blocks += 1
                if arr is None:
                    continue
                info = cols[7]
                end = int(info[4:].split(b";", 1)[0]) if info.startswith(b"END=") else pos
                lo, hi = bisect_left(arr, pos), bisect_right(arr, end)
                if lo < hi:
                    fmt, sample = cols[8].split(b":"), cols[9].split(b":")
                    dp_b = _format_value(fmt, sample, b"MIN_DP") or _format_value(fmt, sample, b"DP")
                    try:
                        dp = int(dp_b) if dp_b is not None else 0
                    except ValueError:
                        dp = 0
                    code = _COV_REF_OK if (filt == "PASS" and dp >= MIN_CONFIRM_DP) else _COV_REF_LOW
                    st = cov[chrom]
                    for i in range(lo, hi):
                        if st[i] < code:
                            st[i] = code
                continue

            variant_rows += 1
            ref = cols[3].decode("ascii", "replace").upper()
            alts = alt_b.decode("ascii", "replace").upper().split(",")
            gt = cols[9].split(b":", 1)[0].decode("ascii", "replace")
            called = _parse_gt(gt)
            ploidy = len(called)
            if arr is not None:
                lo, hi = bisect_left(arr, pos), bisect_right(arr, pos + len(ref) - 1)
                st = cov[chrom]
                for i in range(lo, hi):
                    st[i] = _COV_CALLED

            counts: Dict[int, int] = {}
            for a in called:
                if a:
                    counts[a] = counts.get(a, 0) + 1
            for k, n_copies in counts.items():
                if k > len(alts):
                    continue
                allele = alts[k - 1]
                if not (is_plain_allele(allele) and is_plain_allele(ref)):
                    nonplain += 1
                    continue
                alleles_checked += 1
                npos, nref, nalt = normalize_vcf_allele(pos, ref, allele)
                hits: List[Tuple[int, int, str, str, int, int, str]] = [
                    (v, p, r, a, t, s, "exact") for (v, p, r, a, t, s) in index.exact(chrom, npos, nref, nalt)
                ]
                if len(nref) != len(nalt) and max(len(ref), len(allele)) <= MAX_EQUIV_INDEL:
                    exact_ids = {h[0] for h in hits}
                    for (v, p, r, a, t, s) in index.indels_between(chrom, pos - INDEL_WINDOW, pos + len(ref) + INDEL_WINDOW):
                        if v in exact_ids or max(len(r), len(a)) > MAX_EQUIV_INDEL or not spans_touch(pos, ref, p, r):
                            continue
                        eq = vcf_alleles_equivalent(p, r, a, pos, ref, allele)
                        if eq:
                            hits.append((v, p, r, a, t, s, "equivalent"))
                        elif eq is None and t >= cr.TIER_LIKELY_PATH and (v, pos) not in unresolved_seen:
                            unresolved_seen.add((v, pos))
                            result.review.append({
                                "variation_id": v, "chr": chrom, "clinvar": f"{p} {r}>{a}", "gvcf": f"{pos} {ref}>{allele}",
                                "genotype": gt, "filter": filt,
                                "reason": "overlapping indel representations that cannot be compared without the reference sequence",
                            })
                for (v, p, r, a, t, s, how) in hits:
                    if t == cr.TIER_NONE:
                        continue
                    if t <= cr.TIER_BENIGN:
                        if filt == "PASS":
                            benign_pass += 1
                        else:
                            benign_filtered += 1
                        continue
                    key = (v, pos, ref, k)
                    if key in seen:
                        continue
                    seen[key] = Finding(
                        variation_id=v, chr=chrom, pos=pos, ref=ref, alt=allele, alt_index=k, copies=n_copies,
                        ploidy=ploidy, gt=gt, filter=filt, match=how, tier=t, stars=s,
                        cv_pos=p, cv_ref=r, cv_alt=a, line=raw.rstrip(b"\r\n").decode("utf-8", "replace"),
                    )

    result.findings = sorted(seen.values(), key=lambda x: (-x.tier, -x.stars, x.chr, x.pos))
    by_tier: Dict[str, Dict[str, int]] = {}
    for fd in result.findings:
        slot = by_tier.setdefault(cr.TIER_NAMES[fd.tier], {"PASS": 0, "not_PASS": 0})
        slot["PASS" if fd.passes else "not_PASS"] += 1
    by_tier[cr.TIER_NAMES[cr.TIER_BENIGN]] = {"PASS": benign_pass, "not_PASS": benign_filtered, "cataloged": 0}
    result.stats = {
        "gvcf_lines": lines,
        "reference_blocks": ref_blocks,
        "variant_rows": variant_rows,
        "alt_alleles_evaluated": alleles_checked,
        "non_standard_alts_skipped": nonplain,
        "clinvar_matches_by_tier": by_tier,
        "equivalent_representation_matches": sum(1 for x in result.findings if x.match == "equivalent"),
        "unresolved_high_impact_indels": len(result.review),
        "scan_seconds": round(time.time() - t0, 1),
    }
    tot = sum(len(a) for a in loci.values())
    counts_cov = {"confirmed_reference": 0, "reference_low_confidence": 0, "variant_called": 0, "no_gvcf_record": 0}
    names = {_COV_NONE: "no_gvcf_record", _COV_REF_OK: "confirmed_reference", _COV_REF_LOW: "reference_low_confidence", _COV_CALLED: "variant_called"}
    for c, st in cov.items():
        for code in st:
            counts_cov[names[code]] += 1
    result.coverage = {
        "high_impact_loci": tot,
        "definition": f"distinct GRCh38 positions of ClinVar Pathogenic / Likely pathogenic SNV+indel alleles; reference counted as confirmed only for PASS blocks with depth >= {MIN_CONFIRM_DP}",
        **counts_cov,
    }
    return result


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------

def _genotype_columns(finding: Finding) -> Dict[str, Any]:
    """Catalog columns that describe the sample's call for one matched ClinVar allele."""
    parts = finding.line.split("\t")
    info_map = parse_info_field(parts[7])
    sample_map = parse_sample_field(parts[8], parts[9])
    fmt = format_vcf_genotype_and_status(
        parts[3], parts[4], parts[5], parts[6], info_map, sample_map, "snp-indel", expected_risk_allele=finding.cv_alt
    )
    carried = finding.passes and all(a is not None for a in _parse_gt(finding.gt))
    if not carried:
        status = LOW_QUALITY_STATUS
    elif finding.copies >= finding.ploidy:
        status = "Detected (Homozygous)"
    else:
        status = "Carrier (Heterozygous)"
    return {
        "ref_allele": finding.ref,
        "alt_allele": finding.alt,
        "risk_version": finding.cv_alt,
        "end_position": finding.pos + len(finding.ref) - 1,
        "your_data": fmt["your_data"],
        "your_status": status,
        "genotype": fmt["gt"],
        "alleles_called": fmt["alleles_called"],
        "dp": fmt["dp"],
        "gq": fmt["gq"],
        "ad": fmt["ad"],
        "qual": fmt["qual"],
        "filter_status": finding.filter,
        "is_carried": 1 if carried else 0,
        "raw_vcf_line": finding.line,
        "clinvar_match": finding.match if carried else "low_quality",
    }


def _annotation_columns(ann: Dict[str, Any], finding: Finding, release: str) -> Dict[str, Any]:
    gene = ((ann.get("gene_symbol") or "").split(";")[0].strip()) or "Intergenic"
    rcvs = [r for r in (ann.get("rcvs") or "").split(", ") if r]
    gt_cols = format_vcf_genotype_and_status(
        finding.ref, finding.alt, "0", finding.filter, {}, parse_sample_field(finding.line.split("\t")[8], finding.line.split("\t")[9]), "snp-indel"
    )
    evidence, stars = map_review_status_to_evidence(ann.get("review_status") or "", finding.filter, gt_cols["dp"], gt_cols["gq"])
    vcv = f"VCV{finding.variation_id:09d}"
    return {
        "rcv_id": rcvs[0] if rcvs else vcv,
        "all_rcv_ids": ", ".join(rcvs),
        "gene": gene,
        "gene_name": "",
        "classification": ann.get("classification") or "Uncertain Significance",
        "review_status": ann.get("review_status") or "no classification provided",
        "evidence": f"{evidence} | ClinVar {vcv} ({release} release)",
        "evidence_stars": stars,
        "condition": ann.get("conditions") or f"{gene}-related genetic phenotype",
        "hgvs": ann.get("name") or "",
        "molecular_consequence": ann.get("type") or "",
        "variant_type": determine_variant_type(finding.ref, finding.alt, {}, "snp-indel"),
        "risk_frequency": None,
        "risk_frequency_text": "Not annotated (ClinVar scan)",
    }


def _new_row(finding: Finding, ann: Dict[str, Any], release: str) -> Dict[str, Any]:
    vcv = f"VCV{finding.variation_id:09d}"
    parts = finding.line.split("\t")
    rsid = ann.get("rsid") or (parts[2] if parts[2].startswith("rs") else None) or f"chr{finding.chr}:{finding.pos}"
    row: Dict[str, Any] = {
        "unique_key": f"{SCAN_ROW_PREFIX}{finding.chr}:{finding.pos}:{vcv}",
        "chr": finding.chr,
        "position": finding.pos,
        "variant_id": rsid,
        "source_file": SNP_VCF_NAME,
        "clinvar_variation_id": str(finding.variation_id),
        "clinvar_release": release,
    }
    row.update(_annotation_columns(ann, finding, release))
    row.update(_genotype_columns(finding))
    return row


# ---------------------------------------------------------------------------
# Catalog reconciliation
# ---------------------------------------------------------------------------

@dataclass
class CatalogChanges:
    inserted: int = 0
    linked_updated: int = 0
    adjudicated: int = 0
    unverified: int = 0
    false_carriers_corrected: List[Dict[str, Any]] = field(default_factory=list)
    carriers_added_to_existing: List[Dict[str, Any]] = field(default_factory=list)
    status_changes: int = 0
    scan_vs_adjudication_disagreements: List[Dict[str, Any]] = field(default_factory=list)
    unverified_rows: List[Dict[str, Any]] = field(default_factory=list)
    classification_differs_from_release: List[Dict[str, Any]] = field(default_factory=list)
    skipped_no_annotation: List[int] = field(default_factory=list)
    false_negatives_corrected: List[Dict[str, Any]] = field(default_factory=list)
    stale_classification_rows: List[Dict[str, Any]] = field(default_factory=list)
    linked_by_rsid: int = 0
    legacy_allele_verified: int = 0


def _row_vids(row: Dict[str, Any], rcv_map: Dict[str, int], alleles: Dict[int, Dict[str, Any]]) -> List[int]:
    """ClinVar variation id(s) a catalog row stands for: explicit link first, then via RCV."""
    out: List[int] = []
    for col in ("clinvar_variation_id", "ncbi_variation_id"):
        v = str(row.get(col) or "").strip()
        if v.isdigit() and int(v) not in out:
            out.append(int(v))
    rcvs = [str(row.get("rcv_id") or "")] + [x.strip() for x in str(row.get("all_rcv_ids") or "").split(",")]
    for r in rcvs:
        v = rcv_map.get(r.strip())
        if v and v not in out:
            out.append(v)
    # keep only variations located at this row's locus (a stray RCV must not hijack a row)
    near = [v for v in out if v in alleles and alleles[v]["chr"] == row["chr"] and abs(alleles[v]["pos"] - int(row["position"])) <= 10]
    return near


def _best(findings: Sequence[Finding]) -> Finding:
    return sorted(findings, key=lambda x: (not x.passes, x.match != "exact", -x.copies))[0]


def _adjudicate(row: Dict[str, Any], cv: Dict[str, Any]) -> Dict[str, Any]:
    """Independent, allele-specific read of the gVCF at the row's ClinVar allele
    (genomescope.ncbi_service.build_your_genotype) -> catalog genotype columns."""
    from genomescope import ncbi_service as ns

    gt = ns.build_your_genotype(cv["chr"], cv["pos"], target_vcf=(cv["pos"], cv["ref"], cv["alt"]))
    carried = bool(gt.get("is_carried")) and gt.get("carries_ncbi_alt") is True
    cols: Dict[str, Any] = {
        "ref_allele": cv["ref"],
        "alt_allele": cv["alt"],
        "risk_version": cv["alt"],
        "is_carried": 1 if carried else 0,
        "clinvar_match": "adjudicated_carried" if carried else "adjudicated_not_carried",
    }
    if not gt.get("found"):
        cols.update({
            "your_status": NO_RECORD_STATUS, "your_data": "No gVCF record at this locus", "genotype": "./.",
            "alleles_called": None, "dp": None, "gq": None, "ad": None, "qual": None, "filter_status": None,
            "raw_vcf_line": None, "clinvar_match": "adjudicated_no_record",
        })
        return cols
    parts = gt["raw_vcf_line"].split("\t")
    info_map = parse_info_field(parts[7])
    sample_map = parse_sample_field(parts[8], parts[9])
    shown_ref = ns.locus_ref_for_row(parts, None, (cv["pos"], cv["ref"], cv["alt"]))
    fmt = format_vcf_genotype_and_status(shown_ref, parts[4], parts[5], parts[6], info_map, sample_map, "snp-indel", expected_risk_allele=cv["alt"])
    status = gt["your_status"]
    if carried and gt.get("zygosity_for_ncbi_alt") == "homozygous":
        status = "Detected (Homozygous)"
    elif carried:
        status = "Carrier (Heterozygous)" if gt.get("zygosity_for_ncbi_alt") == "heterozygous" else status
    if not carried and status in ("Carrier (Heterozygous)", "Detected (Homozygous)"):
        status = UNVERIFIED_ALLELE_STATUS if gt.get("carries_ncbi_alt") is None else DIFFERENT_ALT_STATUS
    cols.update({
        "your_status": status, "your_data": fmt["your_data"], "genotype": gt.get("genotype"),
        "alleles_called": gt.get("alleles_called"), "dp": gt.get("dp"), "gq": gt.get("gq"), "ad": gt.get("ad"),
        "qual": fmt["qual"], "filter_status": gt.get("filter"), "raw_vcf_line": gt.get("raw_vcf_line"),
    })
    if gt.get("carries_ncbi_alt") is None and not carried:
        cols["clinvar_match"] = "adjudicated_undecided"
    return cols


def _own_allele(row: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """The substitution/indel a legacy row names (REF and risk allele), when usable."""
    ref, alt = str(row.get("ref_allele") or "").upper(), str(row.get("risk_version") or "").upper()
    if is_plain_allele(ref) and is_plain_allele(alt) and ref != alt:
        return {"chr": row["chr"], "pos": int(row["position"]), "ref": ref, "alt": alt, "tier": -1}
    return None


def _snapshot(row: Dict[str, Any]) -> Dict[str, Any]:
    return {k: row.get(k) for k in ("id", "variant_id", "chr", "position", "gene", "classification", "your_status", "is_carried", "ref_allele", "alt_allele", "risk_version")}


def _update_row(conn: sqlite3.Connection, row_id: int, cols: Dict[str, Any]) -> None:
    keys = list(cols)
    conn.execute(f"UPDATE variants SET {', '.join(k + ' = ?' for k in keys)} WHERE id = ?", [cols[k] for k in keys] + [row_id])


def reconcile_catalog(
    result: ScanResult,
    index: cr.ClinVarIndex,
    release: str,
    min_tier: int = cr.TIER_CONFLICTING,
    summary_gz: Optional[Path] = None,
    dry_run: bool = False,
    progress: Optional[Callable[[str], None]] = None,
) -> CatalogChanges:
    say = progress or (lambda _m: None)
    indexer.init_db()
    ch = CatalogChanges()
    conn = get_db_connection()
    try:
        # Direct VCF reads ("RAW:" rows) make no ClinVar claim and are left alone
        rows = [
            dict(r) for r in conn.execute(
                "SELECT * FROM variants WHERE source_file = ? AND unique_key NOT LIKE 'RAW:%'", (SNP_VCF_NAME,)
            )
        ]
        by_id = {r["id"]: r for r in rows}
        rcv_map = index.variations_for_rcvs(
            [str(r.get("rcv_id") or "") for r in rows] + [x for r in rows for x in str(r.get("all_rcv_ids") or "").split(",")]
        )
        linked_vids: Dict[int, List[int]] = {}      # row id -> candidate vids
        all_vids = set(rcv_map.values())
        for r in rows:
            for col in ("clinvar_variation_id", "ncbi_variation_id"):
                v = str(r.get(col) or "").strip()
                if v.isdigit():
                    all_vids.add(int(v))
        all_vids |= {f.variation_id for f in result.findings}
        alleles = index.alleles_for(all_vids)
        for r in rows:
            linked_vids[r["id"]] = _row_vids(r, rcv_map, alleles)
        # Rows citing RCV accessions that no longer exist (merged / withdrawn since the MyVariant
        # snapshot) are re-linked by dbSNP rsID + locus, choosing the allele the row names.
        unlinked = [r for r in rows if not linked_vids[r["id"]] and not str(r.get("ncbi_variation_id") or "").strip()]
        by_rs = index.alleles_by_rsid(str(r.get("variant_id") or "") for r in unlinked)
        for r in unlinked:
            cands = [c for c in by_rs.get(str(r.get("variant_id") or ""), [])
                     if c["chr"] == r["chr"] and abs(c["pos"] - int(r["position"])) <= 10]
            want = str(r.get("risk_version") or "").upper()
            pick = [c for c in cands if c["alt"] == want] or (cands if len(cands) == 1 else [])
            if len(pick) == 1 or (pick and len({(c["pos"], c["ref"], c["alt"]) for c in pick}) == 1):
                linked_vids[r["id"]] = [pick[0]["variation_id"]]
                alleles[pick[0]["variation_id"]] = {**pick[0], "variation_id": pick[0]["variation_id"]}
                ch.linked_by_rsid += 1

        rows_by_vid: Dict[int, List[int]] = {}
        for rid, vids in linked_vids.items():
            for v in vids:
                rows_by_vid.setdefault(v, []).append(rid)

        groups: Dict[int, List[Finding]] = {}
        for fd in result.findings:
            if fd.tier >= min_tier:
                groups.setdefault(fd.variation_id, []).append(fd)
        anns = index.annotations_for(groups.keys())
        missing = [v for v in groups if v not in anns]
        if missing and summary_gz is not None and Path(summary_gz).exists():
            say(f"  reading annotation for {len(missing):,} lower-tier findings from {Path(summary_gz).name} ...")
            anns.update(cr.annotations_from_summary(Path(summary_gz), missing))
        ch.skipped_no_annotation = [v for v in groups if v not in anns]

        handled_rows: set = set()
        new_rows: List[Dict[str, Any]] = []
        updates: List[Tuple[int, Dict[str, Any]]] = []

        # 1. variations the sample carries (or has a filtered call for)
        for vid, fds in groups.items():
            fd = _best(fds)
            ann = anns.get(vid)
            if ann is None:
                continue
            gcols = _genotype_columns(fd)
            existing = rows_by_vid.get(vid, [])
            if not existing:
                new_rows.append(_new_row(fd, ann, release))
                ch.inserted += 1
                continue
            for rid in existing:
                row = by_id[rid]
                cols = dict(gcols)
                cols.pop("end_position", None)
                cols["clinvar_variation_id"] = str(vid)
                cols["clinvar_release"] = release
                if str(row.get("unique_key", "")).startswith(SCAN_ROW_PREFIX) and not row.get("ncbi_synced_at"):
                    cols.update(_annotation_columns(ann, fd, release))
                if int(row.get("is_carried") or 0) == 0 and cols["is_carried"] == 1:
                    ch.carriers_added_to_existing.append(_snapshot(row) | {"new_status": cols["your_status"]})
                if row.get("your_status") != cols["your_status"]:
                    ch.status_changes += 1
                if row.get("classification") and ann.get("classification") and row["classification"] != ann["classification"] and not str(row.get("unique_key", "")).startswith(SCAN_ROW_PREFIX):
                    ch.classification_differs_from_release.append(_snapshot(row) | {"release_classification": ann["classification"]})
                updates.append((rid, cols))
                handled_rows.add(rid)
                ch.linked_updated += 1

        # 2. every other snp-indel row: independent allele-specific adjudication
        todo = [r for r in rows if r["id"] not in handled_rows]
        for n, row in enumerate(todo, 1):
            if n % 500 == 0:
                say(f"  adjudicated {n:,}/{len(todo):,} existing catalog rows ...")
            vids = linked_vids.get(row["id"], [])
            if not vids and str(row.get("ncbi_variation_id") or "").strip():
                # imported from live NCBI with its own allele-specific check, but not an SNV/indel
                # in the local index (e.g. no classification, complex record): leave as verified
                updates.append((row["id"], {"clinvar_release": release, "clinvar_match": "ncbi_verified"}))
                continue
            if not vids and int(row.get("is_carried") or 0) == 1:
                own = _own_allele(row)
                if own:
                    cols = _adjudicate(row, own)
                    if cols["is_carried"] == 1:
                        cols.update({"clinvar_release": release, "clinvar_match": "legacy_allele_verified"})
                        ch.legacy_allele_verified += 1
                        ch.unverified_rows.append(_snapshot(row) | {"note": "genotype verified against the row's own allele; no matching ClinVar record in this release"})
                        updates.append((row["id"], cols))
                        continue
            if not vids:
                cols = {"clinvar_release": release, "clinvar_match": "unverified"}
                if int(row.get("is_carried") or 0) == 1:
                    cols.update({"is_carried": 0, "your_status": UNVERIFIED_ALLELE_STATUS})
                    ch.false_carriers_corrected.append(_snapshot(row) | {"reason": "no ClinVar allele could be tied to this row in the current release; carried status cannot be verified"})
                ch.unverified += 1
                ch.unverified_rows.append(_snapshot(row))
                updates.append((row["id"], cols))
                continue
            cv = alleles[vids[0]]
            cols = _adjudicate(row, cv)
            cols["clinvar_variation_id"] = str(vids[0])
            cols["clinvar_release"] = release
            ch.adjudicated += 1
            if int(row.get("is_carried") or 0) == 1 and cols["is_carried"] == 0:
                ch.false_carriers_corrected.append(_snapshot(row) | {
                    "reason": "sample does not carry the ClinVar allele",
                    "clinvar_allele": f"{cv['chr']}:{cv['pos']} {cv['ref']}>{cv['alt']}",
                    "gvcf_line": row.get("raw_vcf_line"),
                    "new_status": cols["your_status"],
                })
            if cols["is_carried"] == 1 and cv.get("tier", 0) >= min_tier:
                # both matchers must agree for every allele the scan is responsible for
                ch.scan_vs_adjudication_disagreements.append(_snapshot(row) | {
                    "reason": "adjudication finds the allele carried but the scan produced no finding (investigate)",
                    "clinvar_allele": f"{cv['chr']}:{cv['pos']} {cv['ref']}>{cv['alt']}",
                })
            elif int(row.get("is_carried") or 0) == 0 and cols["is_carried"] == 1:
                ch.false_negatives_corrected.append(_snapshot(row) | {
                    "reason": "sample carries this allele (below the cataloged tier or not found by the legacy join)",
                    "clinvar_allele": f"{cv['chr']}:{cv['pos']} {cv['ref']}>{cv['alt']}", "clinvar_tier": cr.TIER_NAMES.get(cv.get("tier", 0)),
                })
            if cv.get("tier", 0) == cr.TIER_NONE and STALE_MARKER not in str(row.get("evidence") or ""):
                cols["evidence"] = f"{row.get('evidence') or ''} | ClinVar {release}: {STALE_MARKER} (legacy classification not confirmed)".strip(" |")
                ch.stale_classification_rows.append(_snapshot(row) | {"note": "ClinVar now holds no germline classification for this single variant"})
            if row.get("your_status") != cols.get("your_status"):
                ch.status_changes += 1
            updates.append((row["id"], cols))

        if not dry_run:
            for rid, cols in updates:
                _update_row(conn, rid, cols)
            conn.commit()
            if new_rows:
                upsert_records(new_rows)
        say(f"  catalog: {ch.inserted} added, {ch.linked_updated} updated from scan, {ch.adjudicated} re-adjudicated, {ch.unverified} unverified")
        return ch
    finally:
        conn.close()


def verify_catalog(result: ScanResult, min_tier: int = cr.TIER_CONFLICTING, skip_vids: Sequence[int] = ()) -> List[str]:
    """Post-conditions after reconciliation. Returns violations (empty = all hold)."""
    problems: List[str] = []
    conn = get_db_connection()
    try:
        by_vid: Dict[str, List[sqlite3.Row]] = {}
        for r in conn.execute("SELECT id, is_carried, clinvar_variation_id, ncbi_variation_id, clinvar_match FROM variants WHERE source_file = ?", (SNP_VCF_NAME,)):
            for col in ("clinvar_variation_id", "ncbi_variation_id"):
                if r[col]:
                    by_vid.setdefault(str(r[col]), []).append(r)
        for fd in result.findings:
            if fd.tier < min_tier or not fd.passes or fd.variation_id in skip_vids:
                continue
            rows = by_vid.get(str(fd.variation_id), [])
            if not any(r["is_carried"] == 1 for r in rows):
                problems.append(f"VCV{fd.variation_id:09d} chr{fd.chr}:{fd.pos} {fd.ref}>{fd.alt} ({cr.TIER_NAMES[fd.tier]}) is carried per the scan but has no carried catalog row")
        bad = conn.execute(
            "SELECT id, variant_id, chr, position, clinvar_match FROM variants WHERE source_file = ? AND is_carried = 1 "
            "AND unique_key NOT LIKE 'RAW:%' "
            "AND (clinvar_match IS NULL OR clinvar_match NOT IN ('exact','equivalent','adjudicated_carried','ncbi_verified','legacy_allele_verified'))", (SNP_VCF_NAME,)
        ).fetchall()
        for r in bad:
            problems.append(f"catalog row {r['id']} ({r['variant_id']} chr{r['chr']}:{r['position']}) is marked carried without a verified ClinVar allele match ({r['clinvar_match']})")
        dup = conn.execute(
            "SELECT clinvar_variation_id, COUNT(*) c FROM variants WHERE unique_key LIKE ? GROUP BY clinvar_variation_id, chr, position HAVING c > 1", (SCAN_ROW_PREFIX + "%",)
        ).fetchall()
        for r in dup:
            problems.append(f"duplicate scan rows for VCV{r[0]}")
    finally:
        conn.close()
    return problems


# ---------------------------------------------------------------------------
# Reporting / bookkeeping
# ---------------------------------------------------------------------------

def _ensure_runs_table(conn: sqlite3.Connection) -> None:
    conn.execute("CREATE TABLE IF NOT EXISTS clinvar_scan_runs (id INTEGER PRIMARY KEY AUTOINCREMENT, finished_at TEXT, summary TEXT)")


def record_run(summary: Dict[str, Any]) -> None:
    conn = get_db_connection()
    try:
        _ensure_runs_table(conn)
        conn.execute("INSERT INTO clinvar_scan_runs (finished_at, summary) VALUES (?, ?)", (summary["finished_at"], json.dumps(summary)))
        conn.commit()
    finally:
        conn.close()


def last_run() -> Optional[Dict[str, Any]]:
    try:
        conn = get_db_connection()
        try:
            _ensure_runs_table(conn)
            row = conn.execute("SELECT summary FROM clinvar_scan_runs ORDER BY id DESC LIMIT 1").fetchone()
            return json.loads(row[0]) if row else None
        finally:
            conn.close()
    except Exception:
        return None


def status(index_path: Optional[Path] = None) -> Dict[str, Any]:
    """Read-only state of the local ClinVar verification, for GET /api/clinvar/status."""
    index_path = Path(index_path) if index_path else cr.INDEX_PATH
    out: Dict[str, Any] = {"index": None, "last_scan": None, "catalog": None, "state": "index_missing", "message": ""}
    meta: Optional[Dict[str, Any]] = None
    if cr.ClinVarIndex.available(index_path):
        try:
            ix = cr.ClinVarIndex(index_path)
            try:
                meta = ix.meta()
            finally:
                ix.close()
            out["index"] = {k: meta.get(k) for k in (
                "release_date", "release_date_source", "source_md5", "built_at", "rows_indexed",
                "rows_not_vcf_representable", "indexed_by_tier",
            )}
        except cr.ClinVarError as e:
            out["message"] = str(e)
    ws = workspace.get()
    run = last_run() if ws is not None else None
    gvcf = ws.snp_vcf_path if ws is not None else None
    if run:
        out["last_scan"] = {
            "finished_at": run.get("finished_at"), "clinvar_release": run.get("clinvar_release"), "dry_run": run.get("dry_run"),
            "alt_alleles_evaluated": run.get("scan", {}).get("alt_alleles_evaluated"),
            "clinvar_matches_by_tier": run.get("scan", {}).get("clinvar_matches_by_tier"),
            "coverage": run.get("coverage_of_pathogenic_and_likely_pathogenic_loci"),
            "catalog": run.get("catalog"), "postcondition_violations": run.get("postcondition_violations"),
            "unresolved_high_impact_indels": run.get("scan", {}).get("unresolved_high_impact_indels"),
        }
    try:
        if ws is None:
            raise sqlite3.Error("no active analysis")
        conn = get_db_connection()
        try:
            conn.execute("SELECT clinvar_match FROM variants LIMIT 1")
            by_match = {(r[0] or "none"): r[1] for r in conn.execute(
                "SELECT clinvar_match, COUNT(*) FROM variants WHERE source_file = ? AND unique_key NOT LIKE 'RAW:%' GROUP BY clinvar_match", (SNP_VCF_NAME,))}
            unverified_carried = conn.execute(
                "SELECT COUNT(*) FROM variants WHERE source_file = ? AND unique_key NOT LIKE 'RAW:%' AND is_carried = 1 "
                "AND (clinvar_match IS NULL OR clinvar_match NOT IN ('exact','equivalent','adjudicated_carried','ncbi_verified','legacy_allele_verified'))", (SNP_VCF_NAME,)).fetchone()[0]
            out["catalog"] = {"rows_by_clinvar_match": by_match, "carried_rows_without_verified_allele": unverified_carried}
        finally:
            conn.close()
    except sqlite3.Error:
        pass

    if meta is None:
        out["state"] = "index_missing"
        out["message"] = out["message"] or "No local ClinVar release yet. Start a new analysis with 'Download / refresh ClinVar release' selected."
    elif ws is None:
        out["state"] = "no_analysis"
        out["message"] = f"ClinVar release {meta.get('release_date')} is available. Select an analysis to see its verification state."
    elif not run or run.get("dry_run"):
        out["state"] = "never_run"
        out["message"] = "ClinVar index present but this analysis has not been verified. Run a new analysis with ClinVar verification enabled."
    else:
        same_release = run.get("clinvar_release") == meta.get("release_date")
        same_gvcf = gvcf is not None and gvcf.exists() and run.get("gvcf", {}).get("size") == gvcf.stat().st_size and run.get("gvcf", {}).get("mtime") == int(gvcf.stat().st_mtime)
        clean = not run.get("postcondition_violations")
        if same_release and same_gvcf and clean:
            out["state"] = "verified"
            out["message"] = f"Catalog verified against ClinVar release {run.get('clinvar_release')}."
        else:
            out["state"] = "stale"
            out["message"] = "The last scan does not match the current ClinVar index / gVCF, or reported violations. Run a new analysis with ClinVar verification enabled."
    out["checked_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return out


def write_reports(summary: Dict[str, Any], result: ScanResult, ch: Optional[CatalogChanges], index: cr.ClinVarIndex, report_path: Path, tsv_path: Path) -> None:
    report = dict(summary)
    report["review_needed_unresolved_indels"] = result.review
    if ch is not None:
        report["catalog_changes_detail"] = {
            "false_carriers_corrected": ch.false_carriers_corrected,
            "false_negatives_corrected": ch.false_negatives_corrected,
            "stale_classification_rows": ch.stale_classification_rows,
            "carriers_added_to_existing_rows": ch.carriers_added_to_existing,
            "scan_vs_adjudication_disagreements": ch.scan_vs_adjudication_disagreements,
            "unverified_rows": ch.unverified_rows,
            "classification_differs_from_release": ch.classification_differs_from_release,
        }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    anns = index.annotations_for({f.variation_id for f in result.findings})
    with open(tsv_path, "w", encoding="utf-8", newline="") as out:
        out.write("tier\tclassification\tstars\tgene\tvariation_id\trsid\tchr\tpos\tref\talt\tgenotype\tcopies\tfilter\tmatch\tcondition\tname\n")
        for f in result.findings:
            a = anns.get(f.variation_id, {})
            out.write("\t".join(str(x) for x in (
                cr.TIER_NAMES[f.tier], a.get("classification", ""), f.stars, (a.get("gene_symbol") or ""), f.variation_id, a.get("rsid") or "",
                f.chr, f.pos, f.ref, f.alt, f.gt, f.copies, f.filter, f.match, (a.get("conditions") or "").replace("\t", " "), (a.get("name") or "").replace("\t", " "),
            )) + "\n")


def backup_catalog() -> Optional[Path]:
    src = indexer.db_path()
    if not Path(src).exists():
        return None
    bdir = backup_dir()
    bdir.mkdir(parents=True, exist_ok=True)
    dest = bdir / f"genetic_variants.db.pre-clinvar-scan-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    shutil.copy2(src, dest)
    old = sorted(bdir.glob("genetic_variants.db.pre-clinvar-scan-*"))
    for p in old[:-KEEP_BACKUPS]:
        p.unlink(missing_ok=True)
    return dest


def run_scan(
    gvcf: Path,
    index: cr.ClinVarIndex,
    dry_run: bool = False,
    min_tier: int = cr.TIER_CONFLICTING,
    summary_gz: Optional[Path] = None,
    backup: bool = True,
    progress: Optional[Callable[[str], None]] = print,
    report_path: Optional[Path] = None,
    tsv_path: Optional[Path] = None,
) -> Dict[str, Any]:
    say = progress or (lambda _m: None)
    report_path = Path(report_path) if report_path else default_report_path()
    tsv_path = Path(tsv_path) if tsv_path else default_carried_tsv_path()
    meta = index.meta()
    release = meta.get("release_date", "unknown")
    say(f"Scanning {gvcf.name} against ClinVar release {release} ({meta.get('rows_indexed', 0):,} GRCh38 SNV/indel alleles) ...")
    result = scan_gvcf(gvcf, index, progress=say)
    ch: Optional[CatalogChanges] = None
    backup_path: Optional[Path] = None
    problems: List[str] = []
    if not dry_run:
        if backup:
            backup_path = backup_catalog()
            if backup_path:
                say(f"Catalog backup: {backup_path}")
        ch = reconcile_catalog(result, index, release, min_tier=min_tier, summary_gz=summary_gz, progress=say)
        problems = verify_catalog(result, min_tier=min_tier, skip_vids=ch.skipped_no_annotation)
    cataloged = sum(1 for f in result.findings if f.tier >= min_tier)
    summary = {
        "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dry_run": dry_run,
        "clinvar_release": release,
        "clinvar_source_md5": meta.get("source_md5"),
        "clinvar_rows_not_vcf_representable": meta.get("rows_not_vcf_representable"),
        "clinvar_not_vcf_representable_by_type": meta.get("not_vcf_representable_by_type"),
        "gvcf": {"name": gvcf.name, "size": gvcf.stat().st_size, "mtime": int(gvcf.stat().st_mtime)},
        "scan": result.stats,
        "coverage_of_pathogenic_and_likely_pathogenic_loci": result.coverage,
        "findings_considered_for_catalog": cataloged,
        "min_tier_cataloged": cr.TIER_NAMES[min_tier],
        "catalog": None if ch is None else {
            "rows_added": ch.inserted,
            "rows_updated_from_scan": ch.linked_updated,
            "legacy_rows_readjudicated": ch.adjudicated,
            "legacy_rows_unverified": ch.unverified,
            "false_carriers_corrected": len(ch.false_carriers_corrected),
            "false_negatives_corrected": len(ch.false_negatives_corrected),
            "legacy_rows_relinked_by_rsid": ch.linked_by_rsid,
            "legacy_carried_rows_verified_by_own_allele": ch.legacy_allele_verified,
            "stale_classification_rows": len(ch.stale_classification_rows),
            "carriers_added_to_existing_rows": len(ch.carriers_added_to_existing),
            "scan_vs_adjudication_disagreements": len(ch.scan_vs_adjudication_disagreements),
            "status_changes": ch.status_changes,
            "findings_skipped_no_annotation": len(ch.skipped_no_annotation),
            "classification_differs_from_release": len(ch.classification_differs_from_release),
        },
        "postcondition_violations": problems,
        "backup": backup_path.name if backup_path else None,
        "limitations": [
            "Only SNV / indel alleles with a precise GRCh38 VCF representation in ClinVar are matched here; large CNVs, repeat expansions and translocations are not (see clinvar_not_vcf_representable_by_type). Structural variants are covered by the Manta / Canvas files, not by this scan.",
            "Multi-nucleotide variants recorded in ClinVar as one allele are matched only when the gVCF represents them identically; a pair of adjacent SNV calls is not assumed to be the same phased allele.",
            f"Indels written differently from ClinVar are matched by haplotype equivalence only for indels up to {MAX_EQUIV_INDEL} bp whose reference spans overlap; unresolved short cases are listed under review_needed_unresolved_indels. Longer ClinVar deletions / insertions are not compared with small-variant calls.",
            "Benign / likely benign matches are counted but not cataloged. ClinVar classifications are as of the release date and can change.",
            "Research-grade WGS calls are not a clinical diagnosis; confirm clinically relevant findings with an accredited laboratory.",
        ],
    }
    write_reports(summary, result, ch, index, report_path, tsv_path)
    if not dry_run:
        record_run(summary)
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_summary(s: Dict[str, Any]) -> None:
    sc = s["scan"]
    cov = s["coverage_of_pathogenic_and_likely_pathogenic_loci"]
    print("\n=== ClinVar scan summary ===")
    print(f"ClinVar release        : {s['clinvar_release']}  (md5 {s.get('clinvar_source_md5')})")
    print(f"gVCF                   : {s['gvcf']['name']}  ({sc['gvcf_lines']:,} lines, {sc['variant_rows']:,} variant rows)")
    print(f"ALT alleles evaluated  : {sc['alt_alleles_evaluated']:,}   (scan {sc['scan_seconds']} s)")
    print("ClinVar matches (PASS / filtered):")
    for tier, c in sc["clinvar_matches_by_tier"].items():
        print(f"  {tier:<58s} {c['PASS']:>6,} / {c['not_PASS']:>6,}")
    print(f"Equivalent-representation indel matches: {sc['equivalent_representation_matches']}")
    print(f"Unresolved high-impact indels (review): {sc['unresolved_high_impact_indels']}")
    print(f"Pathogenic/LP loci     : {cov['high_impact_loci']:,}  confirmed-ref {cov['confirmed_reference']:,} | low-confidence-ref {cov['reference_low_confidence']:,} | variant-called {cov['variant_called']:,} | no gVCF record {cov['no_gvcf_record']:,}")
    if s["catalog"]:
        print("Catalog:", json.dumps(s["catalog"]))
    if s["postcondition_violations"]:
        print(f"!! {len(s['postcondition_violations'])} POST-CONDITION VIOLATION(S):")
        for p in s["postcondition_violations"][:20]:
            print("   -", p)
    print(f"Report: {default_report_path()}\nCarried list: {default_carried_tsv_path()}")


def prepare_index(
    update: bool = False,
    rebuild: bool = False,
    summary_file: Optional[Path] = None,
    progress: Optional[Callable[[str], None]] = print,
) -> Path:
    """
    Make sure a usable local ClinVar index exists (downloading / rebuilding as requested).
    Returns the variant_summary.txt.gz path used (for VUS annotation look-ups).
    Raises cr.ClinVarError when no index can be provided.
    """
    say = progress or (lambda _m: None)
    summary_path = Path(summary_file) if summary_file else (cr.CLINVAR_DIR / cr.SUMMARY_NAME)
    info = cr.read_release_info(summary_path.parent)
    if update and not summary_file:
        info = cr.download_release(progress=say)
    need_index = rebuild or bool(summary_file) or not cr.ClinVarIndex.available()
    if update and not need_index:
        try:
            ix = cr.ClinVarIndex()
            try:
                need_index = ix.meta().get("source_md5") != info.get("md5")
            finally:
                ix.close()
        except cr.ClinVarError:
            need_index = True  # unreadable / old-schema index
    if need_index:
        if not summary_path.exists():
            raise cr.ClinVarError(f"{summary_path.name} has not been downloaded yet; enable the ClinVar download to fetch it.")
        say("Building the local ClinVar index (several minutes) ...")
        stats = cr.build_index(summary_path, cr.INDEX_PATH, info=info, progress=say)
        say(f"Index built: {stats['rows_indexed']:,} alleles, release {stats['release_date']}")
    return summary_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Scan an analysis' gVCF against a local ClinVar release and reconcile its catalog.")
    ap.add_argument("--workspace", type=Path, required=True, help="analysis job folder (contains inputs/ and work/)")
    ap.add_argument("--update", action="store_true", help="download/verify the newest ClinVar release and rebuild the index first")
    ap.add_argument("--rebuild-index", action="store_true", help="rebuild the index from the already-downloaded release file")
    ap.add_argument("--summary-file", type=Path, help="use this variant_summary.txt.gz instead of downloading (offline)")
    ap.add_argument("--gvcf", type=Path, default=None, help="snp-indel gVCF to scan (default: the workspace's)")
    ap.add_argument("--dry-run", action="store_true", help="scan and write the report only; do not modify the catalog")
    ap.add_argument("--include-vus", action="store_true", help="also add carried Uncertain-significance findings to the catalog (annotation is read from the downloaded release file)")
    ap.add_argument("--no-backup", action="store_true", help="skip the catalog backup")
    args = ap.parse_args(argv)

    ws = workspace.Workspace(args.workspace)
    with workspace.activate(ws):
        gvcf = args.gvcf or ws.snp_vcf_path
        if not gvcf.exists():
            print(f"gVCF not found: {gvcf}", file=sys.stderr)
            return 1
        try:
            summary_path = prepare_index(args.update, args.rebuild_index, args.summary_file)
            index = cr.ClinVarIndex()
            try:
                summary = run_scan(
                    gvcf, index, dry_run=args.dry_run,
                    min_tier=cr.TIER_VUS if args.include_vus else cr.TIER_CONFLICTING,
                    summary_gz=summary_path,
                    backup=not args.no_backup,
                )
            finally:
                index.close()
        except cr.ClinVarError as e:
            print(f"ClinVar error: {e}", file=sys.stderr)
            return 1
        _print_summary(summary)
    return 2 if summary["postcondition_violations"] else 0


if __name__ == "__main__":
    sys.exit(main())
