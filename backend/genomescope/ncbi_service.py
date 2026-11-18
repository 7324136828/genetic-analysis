"""
NCBI service layer: combines live NCBI data (ClinVar, dbSNP, Gene, PubMed) with the
user's own WGS gVCF reads and builds the objects of the /api/ncbi/* API contract.

All NCBI failures are non-fatal: they are collected into the response's "errors" list.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from genomescope import indexer
from genomescope.alleles import DIFFERENT_ALT_STATUS  # noqa: F401  (also re-exported as ncbi_service.DIFFERENT_ALT_STATUS)
from genomescope.indexer import format_frequency, get_db_connection, map_review_status_to_evidence, upsert_records
from genomescope.ncbi_client import (
    EUTILS_BASE,
    NCBI_WEB,
    SIGNIFICANCE_FILTERS,
    VARIATION_BASE,
    NcbiError,
    alleles_equivalent,
    clinvar_sort_key,
    dbsnp_spdis,
    extract_clinvar_spdi,
    gene_symbol_clause,
    get_client,
    parse_clinvar_summary,
    parse_dbsnp_summary,
    parse_gene_summary,
    parse_pubmed_summary,
    parse_refsnp_frequencies,
    record_rsids,
    risk_allele_frequency_from_dbsnp,
    risk_allele_frequency_from_refsnp,
    spdi_to_vcf,
    translate_search_term,
    utc_now_iso,
)
from genomescope.vcf_engine import (
    SNP_VCF_NAME,
    determine_variant_type,
    format_vcf_genotype_and_status,
    normalize_chrom,
    parse_info_field,
    parse_sample_field,
    vcf_engine,
)

Spdi = Tuple[Optional[str], int, str, str]
ClinVarItem = Tuple[Dict[str, Any], Dict[str, Any]]  # (raw esummary, parsed NcbiClinVarRecord)

PLACEHOLDER_CONDITIONS = {"not provided", "not specified", "none provided", "see cases", "genetic condition"}
ANCHOR_FETCH_BUDGET = 6
ANCHOR_CLUSTER_SPAN = 100_000
# ClinVar ids fetched for an rsID: esearch returns the newest first, and the simple
# (often oldest) record must not be cut off by diplotype / haplotype entries.
RSID_RETMAX = 100
RSID_RETMAX_HARD = 500
# Wall-clock budgets (seconds) for endpoint bodies when NCBI is slow / unreachable
LOOKUP_BUDGET = 45.0
SCAN_BUDGET = 60.0
IMPORT_BUDGET = 90.0
FALLBACK_BUDGET = 20.0


def _add_error(errors: List[str], msg: str) -> None:
    msg = str(msg)[:300]
    if msg not in errors:
        errors.append(msg)


# ---------------------------------------------------------------------------
# ClinVar
# ---------------------------------------------------------------------------

def fetch_clinvar(ids: Sequence[str], errors: List[str]) -> List[ClinVarItem]:
    ids = [str(i) for i in ids if str(i).strip()]
    if not ids:
        return []
    try:
        raw_map = get_client().esummary("clinvar", ids)
    except NcbiError as exc:
        _add_error(errors, f"ClinVar esummary: {exc}")
        return []
    out: List[ClinVarItem] = []
    seen = set()
    for i in ids:
        raw = raw_map.get(i)
        if raw is None or i in seen:
            continue
        seen.add(i)
        try:
            out.append((raw, parse_clinvar_summary(raw)))
        except Exception as exc:  # malformed record: skip but report
            _add_error(errors, f"ClinVar record {i} could not be parsed: {type(exc).__name__}")
    return out


def search_clinvar_ids(term: str, retmax: int, errors: List[str]) -> Tuple[int, List[str], str]:
    try:
        res = get_client().esearch("clinvar", term, retmax=retmax)
        return (res["count"], res["ids"], res["translated"] or term)
    except NcbiError as exc:
        _add_error(errors, f"ClinVar esearch: {exc}")
        return (0, [], term)


HIGH_REVIEW_FILTER = '("practice guideline"[RVST] OR "reviewed by expert panel"[RVST] OR "criteria provided, multiple submitters, no conflicts"[RVST])'


def search_clinvar_pool(term: str, pool: int, errors: List[str]) -> Tuple[int, List[str], str]:
    """
    esearch returns newest records first. When a query has more hits than the ranking pool,
    also pull records with 2-4 star review status so well-established variants
    (e.g. HFE C282Y) are not crowded out by recent single-submitter entries.
    """
    total, ids, translated = search_clinvar_ids(term, pool, errors)
    if total > len(ids) and ids:
        prioritized: List[str] = []
        extra_terms = [f"({term}) AND {HIGH_REVIEW_FILTER}"]
        if "clinsig" not in term:
            extra_terms.insert(0, f"({term}) AND {HIGH_REVIEW_FILTER} AND {SIGNIFICANCE_FILTERS['pathogenic']}")
        for extra in extra_terms:
            _t, hi_ids, _q = search_clinvar_ids(extra, max(1, pool // 2), errors)
            prioritized += [i for i in hi_ids if i not in prioritized]
        if prioritized:
            ids = (prioritized + [i for i in ids if i not in prioritized])[:pool]
    return (total, ids, translated)


def fill_clinvar_vcf_alleles(items: Sequence[ClinVarItem], errors: List[str], fetch_budget: int = ANCHOR_FETCH_BUDGET) -> Dict[str, Tuple[int, str, str]]:
    """
    Fill ref/alt (left-aligned VCF style) on parsed ClinVar records. Indels need the
    reference base preceding the event: fetched from NCBI nuccore in clustered windows
    (bounded by fetch_budget). Returns {variation_id: (vcf_pos, ref, alt)}.
    """
    result: Dict[str, Tuple[int, str, str]] = {}
    pending: List[Tuple[Dict[str, Any], Tuple[str, int, str, str], int]] = []
    for raw, rec in items:
        sp = extract_clinvar_spdi(raw)
        if not sp:
            continue
        pos, ref, alt, need = spdi_to_vcf(sp[1], sp[2], sp[3])
        if need is None and pos and ref and alt:
            rec["ref"], rec["alt"] = ref, alt
            result[rec["variation_id"]] = (pos, ref, alt)
        else:
            pending.append((rec, sp, need))  # type: ignore[arg-type]

    # Cluster anchor positions per RefSeq accession
    by_acc: Dict[str, List[Tuple[Dict[str, Any], Tuple[str, int, str, str], int]]] = {}
    for item in pending:
        by_acc.setdefault(item[1][0], []).append(item)
    fetches = 0
    for acc, group in by_acc.items():
        group.sort(key=lambda x: x[2])
        clusters: List[List[Tuple[Dict[str, Any], Tuple[str, int, str, str], int]]] = []
        for item in group:
            if clusters and item[2] - clusters[-1][0][2] <= ANCHOR_CLUSTER_SPAN:
                clusters[-1].append(item)
            else:
                clusters.append([item])
        for cl in clusters:
            lo, hi = cl[0][2], cl[-1][2]
            seq = ""
            if fetches < fetch_budget and lo >= 1:
                fetches += 1
                try:
                    seq = get_client().fetch_sequence(acc, lo, hi)
                except NcbiError as exc:
                    _add_error(errors, f"Reference sequence (anchor base) fetch failed: {exc}")
            for rec, sp, need in cl:
                base = seq[need - lo] if seq and 0 <= need - lo < len(seq) else None
                pos, ref, alt, _n = spdi_to_vcf(sp[1], sp[2], sp[3], anchor_base=base)
                if pos and ref and alt:
                    rec["ref"], rec["alt"] = ref, alt
                    result[rec["variation_id"]] = (pos, ref, alt)
    return result


def condition_text(rec: Dict[str, Any], gene_symbol: Optional[str] = None) -> str:
    names: List[str] = []
    for c in rec.get("conditions") or []:
        n = str(c.get("name") or "").strip()
        if n and n.lower() not in PLACEHOLDER_CONDITIONS and n not in names:
            names.append(n)
    if names:
        return " | ".join(names[:3])
    if gene_symbol:
        return f"{gene_symbol}-related genetic phenotype"
    return "Not provided (NCBI ClinVar)"


def primary_gene_symbol(rec: Dict[str, Any]) -> Optional[str]:
    genes = rec.get("genes") or []
    for g in genes:
        if g.get("symbol") and not str(g["symbol"]).startswith("LOC"):
            return g["symbol"]
    return genes[0]["symbol"] if genes else None


def _primary_gene_id(rec: Optional[Dict[str, Any]], dbsnp: Optional[Dict[str, Any]]) -> Optional[str]:
    for src in (rec, dbsnp):
        if not src:
            continue
        for g in src.get("genes") or []:
            if g.get("gene_id") and not str(g.get("symbol", "")).startswith("LOC"):
                return str(g["gene_id"])
    return None


# ---------------------------------------------------------------------------
# dbSNP / Gene / PubMed
# ---------------------------------------------------------------------------

def _rs_number(rsid: Optional[str]) -> Optional[str]:
    if not rsid:
        return None
    m = re.fullmatch(r"(?i)\s*(?:rs)?(\d{1,12})\s*", str(rsid))
    return m.group(1) if m else None


def fetch_dbsnp(rsids: Sequence[str], errors: List[str]) -> Dict[str, Tuple[Dict[str, Any], Dict[str, Any]]]:
    nums = [n for n in (_rs_number(r) for r in rsids) if n]
    if not nums:
        return {}
    try:
        raw_map = get_client().esummary("snp", nums)
    except NcbiError as exc:
        _add_error(errors, f"dbSNP esummary: {exc}")
        return {}
    out: Dict[str, Tuple[Dict[str, Any], Dict[str, Any]]] = {}
    for n in nums:
        raw = raw_map.get(n)
        if raw is None:
            continue
        try:
            parsed = parse_dbsnp_summary(raw)
        except Exception as exc:
            _add_error(errors, f"dbSNP rs{n} could not be parsed: {type(exc).__name__}")
            continue
        out[f"rs{n}"] = (raw, parsed)
    return out


def find_dbsnp_at_locus(chrom: str, position: int, errors: List[str]) -> Optional[Dict[str, Any]]:
    try:
        res = get_client().esearch("snp", f"{normalize_chrom(chrom)}[CHR] AND {int(position)}[CPOS]", retmax=10)
    except NcbiError as exc:
        _add_error(errors, f"dbSNP esearch: {exc}")
        return None
    found = fetch_dbsnp([f"rs{i}" for i in res["ids"]], errors)
    candidates = [p for _r, p in found.values() if p.get("chr") == normalize_chrom(chrom) and p.get("position") == int(position)]
    if not candidates:
        candidates = [p for _r, p in found.values() if p.get("chr") == normalize_chrom(chrom)]
    candidates.sort(key=lambda p: int(p["rsid"][2:]))
    return candidates[0] if candidates else None


def compute_risk_frequency(
    dbsnp: Optional[Dict[str, Any]],
    target_spdi: Optional[Tuple[str, int, str, str]],
    display_allele: Optional[str],
    errors: List[str],
    allow_refsnp: bool = True,
) -> Optional[Dict[str, Any]]:
    if not dbsnp or not target_spdi:
        return None
    rf = risk_allele_frequency_from_dbsnp(dbsnp, target_spdi, display_allele)
    if rf is None and allow_refsnp and dbsnp.get("frequencies"):
        num = _rs_number(dbsnp.get("rsid"))
        if num:
            try:
                rf = risk_allele_frequency_from_refsnp(parse_refsnp_frequencies(get_client().refsnp(num)), target_spdi, display_allele)
            except NcbiError as exc:
                _add_error(errors, f"Variation Services refsnp: {exc}")
    return rf


def fetch_gene(gene_id: Optional[str], symbol: Optional[str], errors: List[str]) -> Optional[Dict[str, Any]]:
    client = get_client()
    try:
        if not gene_id and symbol:
            res = client.esearch("gene", f'{gene_symbol_clause(symbol).replace("[gene]", "[sym]")} AND "Homo sapiens"[orgn]', retmax=5)
            if not res["ids"]:
                _add_error(errors, f"NCBI Gene: no human gene with symbol {symbol.upper()}")
                return None
            raws = client.esummary("gene", res["ids"])
            parsed = [parse_gene_summary(raws[i]) for i in res["ids"] if i in raws]
            exact = [g for g in parsed if g["symbol"].upper() == symbol.upper()]
            return (exact or parsed or [None])[0]
        if gene_id:
            raws = client.esummary("gene", [gene_id])
            raw = raws.get(str(gene_id))
            return parse_gene_summary(raw) if raw else None
    except NcbiError as exc:
        _add_error(errors, f"NCBI Gene: {exc}")
    return None


def fetch_publications(variation_ids: Sequence[str], rs_number: Optional[str], errors: List[str], limit: int = 10) -> List[Dict[str, Any]]:
    client = get_client()
    pmids: List[str] = []
    try:
        if variation_ids:
            pmids = client.elink("clinvar", "pubmed", variation_ids[:20], linkname="clinvar_pubmed")
        if not pmids and rs_number:
            pmids = client.elink("snp", "pubmed", [rs_number], linkname="snp_pubmed_cited")
    except NcbiError as exc:
        _add_error(errors, f"PubMed elink: {exc}")
        return []
    pmids = sorted(set(pmids), key=lambda x: int(x) if x.isdigit() else 0, reverse=True)[:limit]
    if not pmids:
        return []
    try:
        raws = client.esummary("pubmed", pmids)
    except NcbiError as exc:
        _add_error(errors, f"PubMed esummary: {exc}")
        return []
    return [parse_pubmed_summary(raws[p]) for p in pmids if p in raws]


# ---------------------------------------------------------------------------
# The user's own genotype at an NCBI locus
# ---------------------------------------------------------------------------

def locus_ref_for_row(row: Sequence[str], target_spdi: Optional[Spdi], target_vcf: Optional[Tuple[int, str, str]]) -> str:
    """gVCF hom-ref blocks carry only the block's first base as REF; report the reference
    allele at the queried locus instead so genotypes read e.g. 'G / G' rather than 'A / A'."""
    if row[4] != ".":
        return row[3]
    if target_vcf and target_vcf[1]:
        return target_vcf[1]
    if target_spdi and target_spdi[2]:
        return target_spdi[2]
    return row[3]


def _empty_genotype(chrom: str, position: int) -> Dict[str, Any]:
    return {
        "found": False,
        "chr": chrom,
        "position": int(position),
        "ref": None,
        "alt": None,
        "genotype": None,
        "alleles_called": None,
        "your_status": None,
        "is_carried": False,
        "dp": None,
        "gq": None,
        "ad": None,
        "filter": None,
        "is_reference_block": False,
        "raw_vcf_line": None,
        "source_file": None,
        "carries_ncbi_alt": None,
        "zygosity_for_ncbi_alt": "unknown",
    }


def build_your_genotype(
    chrom: str,
    position: int,
    target_spdi: Optional[Spdi] = None,
    target_vcf: Optional[Tuple[int, str, str]] = None,
    f_handle=None,
    locate_spdis: Optional[Sequence[Spdi]] = None,
) -> Dict[str, Any]:
    """
    Read the analysis' gVCF at (chrom, position) and decide whether the sample carries the
    NCBI ALT allele specifically. Handles gVCF hom-ref blocks (END=), multi-allelic ALT
    lists, no-calls, and indel representation differences (SPDI vs left-aligned VCF).
    """
    chrom = normalize_chrom(chrom)
    if target_spdi is None and target_vcf and target_vcf[1] and target_vcf[2]:
        target_spdi = (None, target_vcf[0] - 1, target_vcf[1].upper(), target_vcf[2].upper())
    lookup_pos = int(target_vcf[0]) if target_vcf else int(position)
    # Alleles used to locate the matching VCF row: the NCBI ALT, or (when there is no single
    # NCBI ALT, e.g. a multi-allelic rsID without ClinVar record) every dbSNP allele.
    targets: List[Spdi] = [target_spdi] if target_spdi else [sp for sp in (locate_spdis or []) if sp]
    longest = max([max(len(t[2]), len(t[3])) for t in targets] or [1])
    # Always collect every row at the position plus a flank: gVCFs split one site into
    # several rows (e.g. an insertion row followed by the SNV row at the same POS), and a
    # 1-bp deletion's VCF anchor sits one base before the deleted base.
    window = longest + 10

    try:
        best, ctx = vcf_engine.lookup_snp_indel_locus(chrom, lookup_pos, f_handle=f_handle, include_context=window)
    except Exception:
        best, ctx = None, []
    if best is None and lookup_pos != int(position):
        try:
            best, ctx2 = vcf_engine.lookup_snp_indel_locus(chrom, int(position), f_handle=f_handle, include_context=window)
            ctx = ctx + ctx2
        except Exception:
            pass

    candidates: List[List[str]] = []
    for r in list(ctx) + ([best] if best else []):
        if r and len(r) >= 10 and r[4] != "." and r not in candidates:
            candidates.append(r)

    matched_row: Optional[List[str]] = None
    overlapping: List[Tuple[Tuple[int, int], List[str]]] = []
    alt_idx: Optional[int] = None
    undecidable = False
    for _seq, p0, d, i in targets:
        t_lo, t_hi = p0, p0 + max(1, len(d))
        for r in candidates:
            v0 = int(r[1]) - 1
            overlaps = v0 < t_hi and (v0 + len(r[3])) > t_lo
            if overlaps and all(r is not o[1] for o in overlapping):
                # prefer PASS calls starting exactly at the target, then any overlapping call
                overlapping.append(((0 if r[6] == "PASS" else 1, 0 if v0 == t_lo else 1), r))
            for k, a in enumerate(r[4].split(","), start=1):
                eq = alleles_equivalent(p0, d, i, int(r[1]), r[3], a)
                if eq:
                    matched_row, alt_idx = r, k
                    break
                if eq is None and overlaps:
                    undecidable = True
            if matched_row:
                break
        if matched_row:
            break

    # Not the NCBI allele: prefer a variant call overlapping the locus (e.g. a different
    # indel spanning it) over the plain row/block lookup at the position.
    overlapping_row = min(overlapping, key=lambda o: o[0])[1] if overlapping else None
    row = matched_row or overlapping_row or best
    if row is None:
        return _empty_genotype(chrom, position)

    info_map = parse_info_field(row[7])
    sample_map = parse_sample_field(row[8], row[9])
    expected = target_vcf[2] if target_vcf else (target_spdi[3] if target_spdi else None)
    if expected is None and matched_row is not None and alt_idx is not None:
        expected = matched_row[4].split(",")[alt_idx - 1]
    shown_ref = locus_ref_for_row(row, target_spdi, target_vcf)
    fmt = format_vcf_genotype_and_status(shown_ref, row[4], row[5], row[6], info_map, sample_map, "snp-indel", expected_risk_allele=expected)
    gt = sample_map.get("GT", "./.") or "./."
    idx = [x for x in re.split(r"[/|]", gt) if x != ""]
    is_ref_block = row[4] == "."

    carries: Optional[bool]
    if not idx or all(x == "." for x in idx):
        zyg, carries = "no_call", None
    elif target_spdi and target_spdi[2] == target_spdi[3]:
        # ClinVar "=" records describe the reference allele itself: having reference copies
        # is not "carrying an ALT", so these are neither counted nor ranked as carried.
        n = sum(1 for x in idx if x == "0")
        zyg, carries = ("absent", False) if n == 0 else ("reference_allele", None)
    elif is_ref_block:
        zyg, carries = "absent", False
    elif alt_idx is not None:
        n = sum(1 for x in idx if x == str(alt_idx))
        if n == 0:
            zyg, carries = "absent", False
        elif len(idx) == 1:
            # haploid call (male chrX/chrY non-PAR, mitochondrial): one copy
            zyg, carries = ("homoplasmic" if chrom == "MT" else "hemizygous"), True
        elif n == len(idx):
            zyg, carries = "homozygous", True
        else:
            zyg, carries = "heterozygous", True
    elif not targets:
        alts_called = [x for x in idx if x not in ("0", ".")]
        if not alts_called:
            zyg = "absent"
        elif len(alts_called) == len(idx) and len(set(alts_called)) == 1:
            zyg = "homozygous"
        else:
            zyg = "heterozygous"
        carries = None
    elif undecidable:
        zyg, carries = "unknown", None
    else:
        zyg, carries = "absent", False

    if row[6] != "PASS" and zyg != "no_call":
        # FAIL / filtered calls (and FAIL hom-ref blocks) can neither confirm nor rule out
        # the NCBI allele; the matched row is still reported for inspection.
        zyg, carries = "low_quality", None

    status = fmt["your_status"]
    if carries is False and fmt["is_carried"]:
        status = DIFFERENT_ALT_STATUS

    report_pos = int(row[1]) if matched_row is not None else int(position)
    return {
        "found": True,
        "chr": chrom,
        "position": report_pos,
        "ref": shown_ref,
        "alt": row[4],
        "genotype": fmt["gt"],
        "alleles_called": fmt["alleles_called"],
        "your_status": status,
        "is_carried": bool(fmt["is_carried"]),
        "dp": fmt["dp"],
        "gq": fmt["gq"],
        "ad": fmt["ad"] if fmt["ad"] not in (None, ".") else None,
        "filter": row[6],
        "is_reference_block": is_ref_block,
        "raw_vcf_line": "\t".join(row),
        "source_file": SNP_VCF_NAME,
        "carries_ncbi_alt": carries,
        "zygosity_for_ncbi_alt": zyg,
    }


def genotype_for_clinvar(raw: Dict[str, Any], rec: Dict[str, Any], vcf_alleles: Dict[str, Tuple[int, str, str]], f_handle=None) -> Optional[Dict[str, Any]]:
    if not rec.get("chr") or not rec.get("start"):
        return None
    sp = extract_clinvar_spdi(raw)
    gt = build_your_genotype(rec["chr"], rec["start"], target_spdi=sp, target_vcf=vcf_alleles.get(rec["variation_id"]), f_handle=f_handle)
    if is_complex_record(raw) and gt.get("carries_ncbi_alt") is True:
        # Haplotype / diplotype records need every member allele (and phase); a single-locus
        # read can rule them out but cannot confirm them.
        gt["carries_ncbi_alt"] = None
        gt["zygosity_for_ncbi_alt"] = "unknown"
    return gt


def genotypes_for_items(items: Sequence[ClinVarItem], vcf_alleles: Dict[str, Tuple[int, str, str]]) -> List[Optional[Dict[str, Any]]]:
    out: List[Optional[Dict[str, Any]]] = []
    fh = None
    try:
        if vcf_engine.snp_vcf_path.exists():
            fh = open(vcf_engine.snp_vcf_path, "rb")
        for raw, rec in items:
            out.append(genotype_for_clinvar(raw, rec, vcf_alleles, f_handle=fh))
    finally:
        if fh:
            fh.close()
    return out


def _is_reference_allele_record(raw: Dict[str, Any]) -> bool:
    sp = extract_clinvar_spdi(raw)
    return bool(sp) and sp[2] == sp[3]


def _rank_with_genotype(pair: Tuple[ClinVarItem, Optional[Dict[str, Any]]]) -> Tuple:
    item, gt = pair
    carried = 0 if (gt and gt.get("carries_ncbi_alt") is True) else 1
    # ClinVar "=" (reference-allele) records come after allele-changing records
    return (carried, 1 if _is_reference_allele_record(item[0]) else 0) + _item_order(item)


# ---------------------------------------------------------------------------
# Endpoint bodies
# ---------------------------------------------------------------------------

def status_body() -> Dict[str, Any]:
    client = get_client()
    ok, err = client.ping()
    return {
        "reachable": ok,
        "api_key_configured": bool(client.api_key),
        "rate_limit_per_sec": int(client.rate_limit_per_sec),
        "cache_entries": client.cache.count(),
        "cache_ttl_hours": client.cache_ttl_hours,
        "services": {"eutils": EUTILS_BASE, "variation": VARIATION_BASE},
        "errors": [err] if err else [],
        "fetched_at": utc_now_iso(),
    }


def health_info() -> Dict[str, Any]:
    client = get_client()
    return {"api_key_configured": bool(client.api_key), "cache_entries": client.cache.count()}


COMPLEX_OBJ_TYPES = {"haplotype", "diplotype", "compoundheterozygote", "phase unknown", "distinct chromosomes", "variation", "genotype"}


def is_complex_record(raw: Dict[str, Any]) -> bool:
    """Haplotypes / compound records describe several alleles; prefer simple records."""
    return len(raw.get("variation_set") or []) > 1 or str(raw.get("obj_type") or "").strip().lower() in COMPLEX_OBJ_TYPES


def _item_order(it: ClinVarItem) -> Tuple:
    return (1 if is_complex_record(it[0]) else 0,) + clinvar_sort_key(it[1])


def _select_primary(
    items: List[ClinVarItem],
    prefer_alt: Optional[str] = None,
    position: Optional[int] = None,
    prefer_variation_id: Optional[str] = None,
    prefer_rcvs: Optional[Sequence[str]] = None,
    carried_ids: Optional[Sequence[str]] = None,
) -> Tuple[Optional[ClinVarItem], str]:
    """
    Choose the record shown first (and used for genotype / frequency / catalog sync).
    Returns (item, how) where how is one of:
      'rcv' | 'variation_id' | 'allele'  - matches the catalog record the caller described
      'carried'                          - a simple record whose allele the user carries
      'position' | 'fallback'            - nothing specific matched (severity order)
    Haplotype / diplotype records are only used when no simple record exists.
    """
    if not items:
        return (None, "")
    simple = [it for it in items if not is_complex_record(it[0])]
    rcvs = [r for r in (prefer_rcvs or []) if r]
    for r in rcvs:
        for it in simple:
            if r in it[1]["rcv_accessions"]:
                return (it, "rcv")
    if rcvs:
        # The catalog's own RCV belongs to a haplotype / diplotype: show that record (the
        # caller refuses to sync complex records) instead of another allele's record.
        for it in items:
            if rcvs[0] in it[1]["rcv_accessions"]:
                return (it, "rcv")
    if prefer_variation_id:
        for it in simple:
            if it[1]["variation_id"] == str(prefer_variation_id):
                return (it, "variation_id")
    if prefer_alt:
        alts = {a.strip().upper() for a in str(prefer_alt).split(",") if a.strip() and a.strip() not in (".", "*")}
        for it in simple:
            if (it[1].get("alt") or "").upper() in alts:
                return (it, "allele")
    if carried_ids:
        for it in simple:
            if it[1]["variation_id"] in carried_ids:
                return (it, "carried")
    pool = simple or items
    if position:
        exact = [it for it in pool if it[1].get("start") == position]
        if exact:
            return (exact[0], "position")
    return (pool[0], "fallback")


def _carried_variation_ids(items: Sequence[ClinVarItem], vcf_alleles: Dict[str, Tuple[int, str, str]]) -> Tuple[List[str], Dict[str, Optional[Dict[str, Any]]]]:
    simple = [it for it in items if not is_complex_record(it[0])]
    gts = genotypes_for_items(simple, vcf_alleles)
    by_id = {it[1]["variation_id"]: gt for it, gt in zip(simple, gts)}
    return ([vid for vid, gt in by_id.items() if gt and gt.get("carries_ncbi_alt") is True], by_id)


def variant_lookup(
    rsid: Optional[str] = None,
    rcv_id: Optional[str] = None,
    vcv_id: Optional[str] = None,
    chrom: Optional[str] = None,
    position: Optional[int] = None,
    include_publications: bool = True,
    prefer_ref: Optional[str] = None,  # kept for API compatibility; block-derived REFs are unreliable
    prefer_alt: Optional[str] = None,
    prefer_variation_id: Optional[str] = None,
    prefer_rcvs: Optional[Sequence[str]] = None,
    meta: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    with get_client().budget(LOOKUP_BUDGET):
        return _variant_lookup(rsid, rcv_id, vcv_id, chrom, position, include_publications, prefer_alt, prefer_variation_id, prefer_rcvs, meta)


def _search_rsid_ids(rsid_norm: str, errors: List[str]) -> List[str]:
    total, ids, _q = search_clinvar_ids(f"{rsid_norm}[VRID]", RSID_RETMAX, errors)
    if total > len(ids) and ids:
        _t2, more, _q2 = search_clinvar_ids(f"{rsid_norm}[VRID]", min(total, RSID_RETMAX_HARD), errors)
        if len(more) > len(ids):
            ids = more
    if total > len(ids):
        _add_error(errors, f"ClinVar lists {total} records for {rsid_norm}; only the {len(ids)} most recent were checked")
    return ids


def _variant_lookup(
    rsid: Optional[str],
    rcv_id: Optional[str],
    vcv_id: Optional[str],
    chrom: Optional[str],
    position: Optional[int],
    include_publications: bool,
    prefer_alt: Optional[str],
    prefer_variation_id: Optional[str],
    prefer_rcvs: Optional[Sequence[str]],
    meta: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    errors: List[str] = []
    rs_num = _rs_number(rsid)
    rsid_norm = f"rs{rs_num}" if rs_num else None
    chrom_n = normalize_chrom(chrom) if chrom else None
    query = {
        "rsid": rsid,
        "rcv_id": rcv_id,
        "vcv_id": vcv_id,
        "chr": chrom,
        "position": position,
        "include_publications": include_publications,
    }

    # 1. Resolve ClinVar variation ids
    cv_ids: List[str] = []
    rcv_norm = None
    if vcv_id:
        m = re.fullmatch(r"(?i)\s*(?:VCV)?0*(\d{1,12})(?:\.\d+)?\s*", vcv_id)
        if m:
            cv_ids.append(m.group(1))
    if rcv_id:
        m = re.fullmatch(r"(?i)\s*RCV0*(\d{1,12})(?:\.\d+)?\s*", rcv_id)
        if m:
            rcv_norm = f"RCV{int(m.group(1)):09d}"
            cv_ids += search_clinvar_ids(rcv_norm, 20, errors)[1]
    if rsid_norm:
        cv_ids += [i for i in _search_rsid_ids(rsid_norm, errors) if i not in cv_ids]
    if chrom_n and position and not cv_ids and not (vcv_id or rcv_id or rsid_norm):
        cv_ids += search_clinvar_ids(f"{chrom_n}[chr] AND {int(position)}[chrpos38]", 20, errors)[1]

    items = fetch_clinvar(cv_ids, errors)
    vcv_uid = cv_ids[0] if (vcv_id and cv_ids) else None
    if vcv_uid or rcv_norm or rsid_norm:
        # keep records matching any of the identifiers the user supplied
        items = [
            it
            for it in items
            if (vcv_uid and it[1]["variation_id"] == vcv_uid)
            or (rcv_norm and rcv_norm in it[1]["rcv_accessions"])
            or (rsid_norm and (it[1].get("rsid") == rsid_norm or rsid_norm in record_rsids(it[0])))
        ]
    elif chrom_n and position:
        items = [it for it in items if it[1].get("chr") == chrom_n and (it[1].get("start") or 0) <= position <= (it[1].get("stop") or it[1].get("start") or 0)] or items
    items.sort(key=_item_order)

    vcf_alleles = fill_clinvar_vcf_alleles(items, errors)
    primary, how = _select_primary(items, prefer_alt, None, prefer_variation_id, prefer_rcvs)
    simple_gts: Dict[str, Optional[Dict[str, Any]]] = {}
    if how == "fallback":
        # No catalog-specific match: show the record for the allele the user actually
        # carries (e.g. CYP2C19*2 c.681G>A rather than a rarer c.681G>T record).
        carried_ids, simple_gts = _carried_variation_ids(items, vcf_alleles)
        primary, how = _select_primary(items, None, position, None, None, carried_ids)
    if meta is not None:
        meta["how"] = how
        meta["primary_complex"] = bool(primary) and is_complex_record(primary[0])
    primary_rec = primary[1] if primary else None
    primary_spdi = extract_clinvar_spdi(primary[0]) if primary else None

    # 2. dbSNP
    dbsnp: Optional[Dict[str, Any]] = None
    rs_for_dbsnp = rsid_norm or (primary_rec.get("rsid") if primary_rec else None)
    if rs_for_dbsnp:
        found = fetch_dbsnp([rs_for_dbsnp], errors)
        if rs_for_dbsnp in found:
            dbsnp = found[rs_for_dbsnp][1]
    elif chrom_n and position:
        dbsnp = find_dbsnp_at_locus(chrom_n, position, errors)

    # 3. Resolved locus
    r_chr = (primary_rec or {}).get("chr") or (dbsnp or {}).get("chr") or chrom_n
    r_pos = (primary_rec or {}).get("start") or (dbsnp or {}).get("position") or position
    r_rsid = (primary_rec or {}).get("rsid") or (dbsnp or {}).get("rsid") or rsid_norm

    # 4. The user's genotype at the locus
    your_genotype = None
    if r_chr and r_pos:
        if primary:
            cached_gt = simple_gts.get(primary_rec["variation_id"]) if primary_rec else None
            your_genotype = cached_gt or genotype_for_clinvar(primary[0], primary_rec, vcf_alleles) or None
        else:
            sps = dbsnp_spdis({"spdi": dbsnp.get("spdi")}) if dbsnp else []
            target_spdi = sps[0] if len(sps) == 1 else None
            your_genotype = build_your_genotype(r_chr, int(r_pos), target_spdi=target_spdi, locate_spdis=sps)

    # 5. Risk allele frequency (frequency of the ClinVar ALT, not the minor allele blindly)
    risk = None
    if dbsnp:
        if primary_spdi:
            risk = compute_risk_frequency(dbsnp, primary_spdi, primary_rec.get("alt") if primary_rec else None, errors)
        elif len(dbsnp.get("alts") or []) == 1:
            sps = dbsnp_spdis({"spdi": dbsnp.get("spdi")})
            if sps:
                risk = compute_risk_frequency(dbsnp, sps[0], dbsnp["alts"][0], errors)

    # 6. Gene + publications
    gene = None
    gid = _primary_gene_id(primary_rec, dbsnp)
    if gid:
        gene = fetch_gene(gid, None, errors)
    publications: List[Dict[str, Any]] = []
    if include_publications and (items or r_rsid):
        publications = fetch_publications([it[1]["variation_id"] for it in items], _rs_number(r_rsid), errors)

    links = {
        "clinvar": primary_rec["url"] if primary_rec else (f"{NCBI_WEB}/clinvar/?term={rcv_norm}" if rcv_norm else None),
        "dbsnp": dbsnp["url"] if dbsnp else (f"{NCBI_WEB}/snp/{r_rsid}" if r_rsid else None),
        "gene": gene["url"] if gene else None,
        "pubmed": None,
    }
    if publications:
        if items:
            links["pubmed"] = f"{NCBI_WEB}/pubmed?LinkName=clinvar_pubmed&from_uid={primary_rec['variation_id'] if primary_rec else items[0][1]['variation_id']}"
        elif r_rsid:
            links["pubmed"] = f"{NCBI_WEB}/pubmed?LinkName=snp_pubmed_cited&from_uid={_rs_number(r_rsid)}"

    if not items and not dbsnp and not errors:
        _add_error(errors, "No NCBI ClinVar or dbSNP record found for this query")

    ordered = [primary] + [it for it in items if it is not primary] if primary else items
    return {
        "query": query,
        "resolved": {"chr": r_chr, "position": int(r_pos) if r_pos else None, "rsid": r_rsid},
        "clinvar": [it[1] for it in ordered],
        "dbsnp": dbsnp,
        "gene": gene,
        "publications": publications[:10],
        "your_genotype": your_genotype,
        "risk_allele_frequency": risk,
        "links": links,
        "errors": errors,
        "fetched_at": utc_now_iso(),
    }


RISK_KEYWORDS = ("risk factor", "risk allele", "association", "protective", "affects")


def _significance_term(significance: str) -> str:
    return SIGNIFICANCE_FILTERS.get(significance, "")


def _post_filter(items: List[ClinVarItem], significance: str) -> List[ClinVarItem]:
    if significance == "risk":
        return [it for it in items if any(k in it[1]["classification"].lower() for k in RISK_KEYWORDS)]
    return items


def _pool_size(limit: int) -> int:
    """How many ClinVar ids to fetch for client-side ranking (carriers first, then significance/stars)."""
    return max(limit, min(200, max(100, limit * 2)))


USER_LOCI_MAX_POSITIONS = 180
USER_LOCI_CHUNK = 60
USER_LOCI_MAX_IDS = 100
USER_LOCI_MAX_SPAN = 3_000_000
USER_LOCI_PADDING = 2_000


def user_carried_positions(gene: Optional[Dict[str, Any]], max_positions: int = USER_LOCI_MAX_POSITIONS) -> Tuple[Optional[str], List[int]]:
    """GRCh38 positions inside the gene where the user's gVCF has a non-reference call.
    Indels also contribute POS+1 (ClinVar places deletions after the VCF anchor base)."""
    if not gene or not gene.get("chr") or not gene.get("start") or not gene.get("stop"):
        return (None, [])
    chrom = normalize_chrom(str(gene["chr"]))
    start, stop = int(gene["start"]), int(gene["stop"])
    if stop - start > USER_LOCI_MAX_SPAN:
        return (chrom, [])
    rows = vcf_engine.read_vcf_region(chrom, max(1, start - USER_LOCI_PADDING), stop + USER_LOCI_PADDING, variants_only=True, limit=2000)
    positions: List[int] = []
    for r in rows:
        if r.get("source_file") != SNP_VCF_NAME or not r.get("is_carried"):
            continue
        p = int(r["position"])
        ref = str(r.get("ref") or "")
        candidates = [p]
        if any(len(a) != len(ref) for a in str(r.get("alt") or "").split(",") if a and a != "."):
            candidates.append(p + 1)
        for c in candidates:
            if c not in positions:
                positions.append(c)
        if len(positions) >= max_positions:
            break
    return (chrom, positions[:max_positions])


def search_user_loci(term: str, gene: Optional[Dict[str, Any]], errors: List[str], cap: int = USER_LOCI_MAX_IDS) -> List[str]:
    """
    ClinVar records (matching `term`) located where the user carries a non-reference allele.
    Large genes have thousands of ClinVar records and esearch returns the newest first, so
    without this the user's own variants (often older, common alleles) never reach the
    carriers-first ranking.
    """
    try:
        chrom, positions = user_carried_positions(gene)
    except Exception as exc:  # VCF read problems are non-fatal
        _add_error(errors, f"VCF read for gene region failed: {type(exc).__name__}")
        return []
    ids: List[str] = []
    for i in range(0, len(positions), USER_LOCI_CHUNK):
        chunk = positions[i:i + USER_LOCI_CHUNK]
        clause = " OR ".join(f"{p}[chrpos38]" for p in chunk)
        _t, hit, _q = search_clinvar_ids(f"({term}) AND {chrom}[chr] AND ({clause})", cap, errors)
        ids += [x for x in hit if x not in ids]
        if len(ids) >= cap:
            break
    return ids[:cap]


def gene_variants(symbol: str, significance: str = "pathogenic", limit: int = 40) -> Dict[str, Any]:
    with get_client().budget(SCAN_BUDGET):
        return _gene_variants(symbol, significance, limit)


def _gene_variants(symbol: str, significance: str, limit: int) -> Dict[str, Any]:
    errors: List[str] = []
    gene = fetch_gene(None, symbol, errors)
    term = gene_symbol_clause(symbol)
    sig_term = _significance_term(significance)
    if sig_term:
        term = f"{term} AND {sig_term}"
    total, ids, _tq = search_clinvar_pool(term, _pool_size(limit), errors)
    if total > len(ids):
        own = search_user_loci(term, gene, errors)
        ids = own + [i for i in ids if i not in own]
    items = _post_filter(fetch_clinvar(ids, errors), significance)
    vcf_alleles = fill_clinvar_vcf_alleles(items, errors)
    gts = genotypes_for_items(items, vcf_alleles)
    pairs = sorted(zip(items, gts), key=_rank_with_genotype)[:limit]
    variants = [{"clinvar": it[1], "your_genotype": gt} for it, gt in pairs]
    return {
        "gene": gene,
        "significance": significance,
        "clinvar_total": total,
        "variants": variants,
        "carried_count": sum(1 for v in variants if v["your_genotype"] and v["your_genotype"].get("carries_ncbi_alt") is True),
        "errors": errors,
        "fetched_at": utc_now_iso(),
    }


def _run_translated_search(term: str, retmax: int, errors: List[str]) -> Tuple[int, List[str], str, str]:
    queries = translate_search_term(term)
    if not queries:
        return (0, [], "", "")
    total, ids, translated, kind = 0, [], queries[0][1], queries[0][0]
    for kind, q in queries:
        total, ids, translated = search_clinvar_pool(q, retmax, errors)
        if total > 0:
            break
    return (total, ids, translated, kind)


def search(term: str, limit: int = 25) -> Dict[str, Any]:
    with get_client().budget(LOOKUP_BUDGET):
        return _search(term, limit)


def _search(term: str, limit: int) -> Dict[str, Any]:
    errors: List[str] = []
    total, ids, translated, kind = _run_translated_search(term, _pool_size(limit), errors)
    items = fetch_clinvar(ids, errors)
    if kind == "rcv":
        rcv = translate_search_term(term)[0][1]
        items = [it for it in items if rcv in it[1]["rcv_accessions"]] or items
    vcf_alleles = fill_clinvar_vcf_alleles(items, errors)
    gts = genotypes_for_items(items, vcf_alleles)
    pairs = sorted(zip(items, gts), key=_rank_with_genotype)[:limit]
    return {
        "term": term,
        "translated_query": translated,
        "total": total,
        "items": [{"clinvar": it[1], "your_genotype": gt} for it, gt in pairs],
        "errors": errors,
        "fetched_at": utc_now_iso(),
    }


# ---------------------------------------------------------------------------
# Catalog (SQLite) integration
# ---------------------------------------------------------------------------

def _frequency_text(rf: Optional[Dict[str, Any]]) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    if not rf:
        return (None, None, None)
    val, text = format_frequency(rf["frequency"])
    return (val, f"{text} [dbSNP {rf['study']}]", f"NCBI dbSNP {rf['study']}")


def _evidence(review_status: str, filt: Optional[str], dp: Optional[int], gq: Optional[int], vcv: Optional[str]) -> Tuple[str, int]:
    text, stars = map_review_status_to_evidence(review_status, filt or ".", dp, gq)
    return (f"{text} | NCBI ClinVar {vcv or ''}".rstrip(), stars)


def _ncbi_fields_for_row(rec: Dict[str, Any], rf: Optional[Dict[str, Any]], filt: Optional[str], dp: Optional[int], gq: Optional[int], gene_symbol: Optional[str]) -> Dict[str, Any]:
    evidence, stars = _evidence(rec["review_status"], filt, dp, gq, rec.get("vcv_id"))
    fields: Dict[str, Any] = {
        "classification": rec["classification_normalized"],
        "review_status": rec["review_status"] or "no classification provided",
        "evidence": evidence,
        "evidence_stars": stars,
        "condition": condition_text(rec, gene_symbol),
        "all_rcv_ids": ", ".join(rec["rcv_accessions"]),
        "ncbi_variation_id": rec["variation_id"],
        "ncbi_last_evaluated": rec.get("last_evaluated"),
        "ncbi_synced_at": utc_now_iso(),
    }
    val, text, source = _frequency_text(rf)
    if val is not None:
        fields.update({"risk_frequency": val, "risk_frequency_text": text, "frequency_source": source})
    return fields


def _update_row(record_id: int, fields: Dict[str, Any]) -> None:
    indexer.ensure_ncbi_columns()
    cols = list(fields.keys())
    conn = get_db_connection()
    try:
        conn.execute(f"UPDATE variants SET {', '.join(c + ' = ?' for c in cols)} WHERE id = ?", [fields[c] for c in cols] + [record_id])
        conn.commit()
    finally:
        conn.close()


def _get_row(record_id: int) -> Optional[Dict[str, Any]]:
    conn = get_db_connection()
    try:
        row = conn.execute("SELECT * FROM variants WHERE id = ?", (record_id,)).fetchone()
        return dict(row) if row else None
    finally:
        conn.close()


def _find_existing_row(rec: Dict[str, Any], raw: Dict[str, Any], vcf_pos: Optional[int], variant_id: str) -> Optional[int]:
    """
    Catalog row that already represents this ClinVar variation: the row synced/imported
    from it (ncbi_variation_id), else a not-yet-linked row at the same locus whose RCV
    belongs to this variation. Haplotype / diplotype records never take over a
    single-variant row (they share the rsid / position / ALT of one member allele).
    """
    indexer.ensure_ncbi_columns()
    conn = get_db_connection()
    try:
        row = conn.execute("SELECT id FROM variants WHERE ncbi_variation_id = ? ORDER BY id LIMIT 1", (rec["variation_id"],)).fetchone()
        if row:
            return int(row[0])
        # a row already verified against the local ClinVar release (genomescope/clinvar_scan.py)
        row = conn.execute("SELECT id FROM variants WHERE clinvar_variation_id = ? ORDER BY id LIMIT 1", (rec["variation_id"],)).fetchone()
        if row:
            return int(row[0])
        if is_complex_record(raw):
            return None
        rcvs = set(rec.get("rcv_accessions") or [])
        if not rcvs:
            return None
        positions = sorted({p for p in (rec.get("start"), vcf_pos) if p})
        for p in positions:
            rows = conn.execute(
                "SELECT id, rcv_id, all_rcv_ids FROM variants WHERE chr = ? AND position = ? AND source_file = ? AND variant_id = ? "
                "AND (alt_allele = ? OR risk_version = ?) AND (ncbi_variation_id IS NULL OR ncbi_variation_id = '') ORDER BY id",
                (rec["chr"], p, SNP_VCF_NAME, variant_id, rec.get("alt") or "", rec.get("alt") or ""),
            ).fetchall()
            for r in rows:
                row_rcvs = {str(r["rcv_id"] or "").strip()} | {x.strip() for x in str(r["all_rcv_ids"] or "").split(",")}
                if row_rcvs & rcvs:
                    return int(r["id"])
        return None
    finally:
        conn.close()


def _row_id_by_unique_key(unique_key: str) -> Optional[int]:
    conn = get_db_connection()
    try:
        row = conn.execute("SELECT id FROM variants WHERE unique_key = ?", (unique_key,)).fetchone()
        return int(row[0]) if row else None
    finally:
        conn.close()


def _variant_type_for(rec: Dict[str, Any]) -> str:
    if rec.get("ref") and rec.get("alt"):
        return determine_variant_type(rec["ref"], rec["alt"], {}, "snp-indel")
    vt = (rec.get("variant_type") or "").lower()
    if "deletion" in vt:
        return "Deletion"
    if "insertion" in vt or "duplication" in vt:
        return "Insertion"
    if "single nucleotide" in vt:
        return "SNV"
    return "Indel"


def import_variation_ids(variation_ids: Sequence[str], max_refsnp_fallbacks: int = 10) -> Dict[str, Any]:
    # Import writes NCBI data into the catalog: fetch fresh records, not week-old cache
    with get_client().budget(IMPORT_BUDGET, fresh=True):
        return _import_variation_ids(variation_ids, max_refsnp_fallbacks)


def _import_variation_ids(variation_ids: Sequence[str], max_refsnp_fallbacks: int) -> Dict[str, Any]:
    errors: List[str] = []
    skipped: List[Dict[str, str]] = []
    ids: List[str] = []
    for v in variation_ids:
        m = re.fullmatch(r"(?i)\s*(?:VCV)?0*(\d{1,12})(?:\.\d+)?\s*", str(v))
        if not m:
            skipped.append({"variation_id": str(v), "reason": "invalid ClinVar variation id"})
            continue
        if m.group(1) not in ids:
            ids.append(m.group(1))

    items = fetch_clinvar(ids, errors)
    got = {it[1]["variation_id"] for it in items}
    for i in ids:
        if i not in got:
            skipped.append({"variation_id": i, "reason": "not found in NCBI ClinVar" if not errors else "NCBI ClinVar unavailable"})

    located: List[ClinVarItem] = []
    for it in items:
        if not it[1].get("chr") or not it[1].get("start"):
            skipped.append({"variation_id": it[1]["variation_id"], "reason": "no GRCh38 location in ClinVar"})
        else:
            located.append(it)

    vcf_alleles = fill_clinvar_vcf_alleles(located, errors)
    dbsnp_map = fetch_dbsnp([it[1]["rsid"] for it in located if it[1].get("rsid")], errors)
    gts = genotypes_for_items(located, vcf_alleles)

    record_ids: List[int] = []
    refsnp_used = 0
    new_rows: List[Dict[str, Any]] = []
    updates: List[Tuple[int, Dict[str, Any]]] = []
    for (raw, rec), gt in zip(located, gts):
        if not gt or not gt.get("found"):
            skipped.append({"variation_id": rec["variation_id"], "reason": f"no VCF coverage at chr{rec['chr']}:{rec['start']}"})
            continue
        dbsnp = dbsnp_map.get(rec.get("rsid") or "", (None, None))[1]
        allow = refsnp_used < max_refsnp_fallbacks
        rf = compute_risk_frequency(dbsnp, extract_clinvar_spdi(raw), rec.get("alt"), errors, allow_refsnp=allow)
        if rf is None and dbsnp and allow:
            refsnp_used += 1
        gene_symbol = primary_gene_symbol(rec)
        fields = _ncbi_fields_for_row(rec, rf, gt.get("filter"), gt.get("dp"), gt.get("gq"), gene_symbol)
        vcf_pos, _vref, _valt = vcf_alleles.get(rec["variation_id"], (rec["start"], rec.get("ref"), rec.get("alt")))
        variant_id = rec.get("rsid") or rec.get("vcv_id") or f"chr{rec['chr']}:{rec['start']}"
        existing = _find_existing_row(rec, raw, vcf_pos, variant_id)
        # Only a confirmed (PASS) call of this exact allele counts as carried; complex
        # records and undecidable loci (carries_ncbi_alt None) are not.
        carried = bool(gt.get("is_carried")) and gt.get("carries_ncbi_alt") is True

        info_map = parse_info_field(gt["raw_vcf_line"].split("\t")[7]) if gt.get("raw_vcf_line") else {}
        sample_map = parse_sample_field(*gt["raw_vcf_line"].split("\t")[8:10]) if gt.get("raw_vcf_line") else {}
        vparts = gt["raw_vcf_line"].split("\t")
        shown_ref = locus_ref_for_row(vparts, extract_clinvar_spdi(raw), vcf_alleles.get(rec["variation_id"]))
        fmt = format_vcf_genotype_and_status(shown_ref, vparts[4], vparts[5], vparts[6], info_map, sample_map, "snp-indel", expected_risk_allele=rec.get("alt"))
        your_data = fmt["your_data"]

        full = dict(fields)
        full.update(
            {
                "rcv_id": rec["rcv_accessions"][0] if rec["rcv_accessions"] else (rec.get("vcv_id") or f"VCV{rec['variation_id']}"),
                "chr": rec["chr"],
                "position": int(vcf_pos or rec["start"]),
                "end_position": int(rec.get("stop") or rec["start"]),
                "variant_id": variant_id,
                "gene": gene_symbol or "Intergenic",
                "gene_name": "",
                "ref_allele": rec.get("ref") or gt.get("ref"),
                "alt_allele": rec.get("alt") or "",
                "risk_version": rec.get("alt") or "",
                "variant_type": _variant_type_for(rec),
                "your_data": your_data,
                "your_status": gt.get("your_status"),
                "genotype": gt.get("genotype"),
                "alleles_called": gt.get("alleles_called"),
                "hgvs": rec.get("title") or "",
                "molecular_consequence": ", ".join(rec.get("molecular_consequences") or []) or rec.get("variant_type") or "",
                "dp": gt.get("dp"),
                "gq": gt.get("gq"),
                "ad": gt.get("ad"),
                "qual": fmt["qual"],
                "filter_status": gt.get("filter"),
                "is_carried": 1 if carried else 0,
                "source_file": SNP_VCF_NAME,
                "raw_vcf_line": gt.get("raw_vcf_line"),
            }
        )
        if "risk_frequency" not in full:
            full["risk_frequency"], full["risk_frequency_text"], full["frequency_source"] = None, "Rare / Not Reported", None
        if existing is not None:
            # Keep the existing row id; refresh its NCBI-derived annotation and VCF read
            upd = {k: v for k, v in full.items() if k not in ("rcv_id", "variant_id", "chr", "position")}
            if "risk_frequency" not in fields:
                for k in ("risk_frequency", "risk_frequency_text", "frequency_source"):
                    upd.pop(k, None)
            updates.append((existing, upd))
        else:
            full["unique_key"] = f"NCBI:{rec['chr']}:{full['position']}:VCV{int(rec['variation_id']):09d}"
            new_rows.append(full)

    for rid, upd in updates:
        _update_row(rid, upd)
        record_ids.append(rid)
    if new_rows:
        upsert_records(new_rows)
        for r in new_rows:
            rid = _row_id_by_unique_key(r["unique_key"])
            if rid is not None:
                record_ids.append(rid)

    return {
        "imported": len(record_ids),
        "record_ids": record_ids,
        "skipped": skipped,
        "errors": errors,
        "fetched_at": utc_now_iso(),
    }


def _freq_match(a: Optional[float], b: Optional[float]) -> bool:
    if a is None or b is None:
        return False
    diff = abs(float(a) - float(b))
    return diff <= 0.005 or diff <= 0.2 * max(abs(float(a)), abs(float(b)))


def _vcf_alt_of_row(row: Dict[str, Any]) -> Optional[str]:
    parts = str(row.get("raw_vcf_line") or "").split("\t")
    if len(parts) >= 5 and parts[4] not in ("", "."):
        return parts[4]
    return None


def _comparisons(row: Dict[str, Any], body: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Catalog vs NCBI. match is None (shown as n/a) when NCBI has no value to compare."""
    rec = body["clinvar"][0] if body.get("clinvar") else None
    rf = body.get("risk_allele_frequency")
    yg = body.get("your_genotype") or {}
    out: List[Dict[str, Any]] = []

    cat_allele = row.get("risk_version") or row.get("alt_allele")
    vcf_alt = _vcf_alt_of_row(row)
    cat_allele_text = f"{cat_allele or '?'}" + (f" (your VCF ALT {vcf_alt})" if vcf_alt and vcf_alt != cat_allele else "")
    ncbi_allele = f"{rec.get('ref') or '?'}>{rec.get('alt') or '?'}" if rec else None
    if not rec or not rec.get("alt"):
        allele_match: Optional[bool] = None
    elif yg.get("carries_ncbi_alt") is False and yg.get("is_carried"):
        allele_match = False  # the user's call is a different ALT than the ClinVar record
    elif yg.get("carries_ncbi_alt") is True:
        allele_match = True
    else:
        allele_match = (cat_allele or "").upper() == (rec.get("alt") or "").upper()
    out.append({"field": "allele", "catalog": cat_allele_text, "ncbi": ncbi_allele, "match": allele_match})

    cat_class = row.get("classification")
    ncbi_class = rec["classification_normalized"] if rec else None
    out.append({"field": "classification", "catalog": cat_class, "ncbi": ncbi_class, "match": None if not rec else (cat_class or "").lower() == (ncbi_class or "").lower()})
    cat_rev = row.get("review_status")
    ncbi_rev = rec["review_status"] if rec else None
    out.append({"field": "review_status", "catalog": cat_rev, "ncbi": ncbi_rev, "match": None if not rec else (cat_rev or "").strip().lower() == (ncbi_rev or "").strip().lower()})
    cat_cond = row.get("condition") or ""
    ncbi_cond = condition_text(rec, primary_gene_symbol(rec)) if rec else None
    cond_match: Optional[bool] = None
    if rec:
        names = [c["name"].lower() for c in rec.get("conditions") or [] if c.get("name") and c["name"].lower() not in PLACEHOLDER_CONDITIONS]
        cl = cat_cond.lower()
        cond_match = (bool(names) and any(n in cl or (cl and cl in n) for n in names)) or cl == (ncbi_cond or "").lower()
    out.append({"field": "condition", "catalog": cat_cond, "ncbi": ncbi_cond, "match": cond_match})
    cat_rcv = row.get("rcv_id")
    ncbi_rcvs = rec["rcv_accessions"] if rec else []
    out.append({"field": "rcv_id", "catalog": cat_rcv, "ncbi": ncbi_rcvs, "match": None if not rec else (bool(cat_rcv) and cat_rcv in ncbi_rcvs)})
    ncbi_f = rf["frequency"] if rf else None
    out.append({"field": "risk_frequency", "catalog": row.get("risk_frequency"), "ncbi": ncbi_f, "match": None if ncbi_f is None else _freq_match(row.get("risk_frequency"), ncbi_f)})
    return out


STRUCTURAL_RCV_PREFIXES = ("RCV-CNV-", "RCV-SV-", "RCV-LOH-")


def is_structural_row(row: Dict[str, Any]) -> bool:
    """CNV / SV / LOH catalog rows (Canvas / Manta) - not single-locus NCBI lookups."""
    src = str(row.get("source_file") or "")
    alt = str(row.get("alt_allele") or "")
    rcv = str(row.get("rcv_id") or "").upper()
    vt = str(row.get("variant_type") or "").lower()
    return (
        (bool(src) and src != SNP_VCF_NAME)
        or alt.startswith("<")
        or rcv.startswith(STRUCTURAL_RCV_PREFIXES)
        or vt.startswith(("cnv", "sv ", "loh"))
    )


def _structural_body(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "query": {"rsid": None, "rcv_id": None, "vcv_id": None, "chr": row.get("chr"), "position": row.get("position"), "include_publications": False},
        "resolved": {"chr": row.get("chr"), "position": row.get("position"), "rsid": None},
        "clinvar": [],
        "dbsnp": None,
        "gene": None,
        "publications": [],
        "your_genotype": None,
        "risk_allele_frequency": None,
        "links": {"clinvar": None, "dbsnp": None, "gene": None, "pubmed": None},
        "errors": [
            f"NCBI single-locus lookup is not applicable to structural / copy-number variants "
            f"({row.get('variant_type') or 'SV/CNV'} chr{row.get('chr')}:{row.get('position')}-{row.get('end_position')}); "
            "the breakpoint base is not the variant."
        ],
        "fetched_at": utc_now_iso(),
    }


def _sync_refusal(row: Dict[str, Any], body: Dict[str, Any], meta: Dict[str, Any]) -> Optional[str]:
    """Why the NCBI primary record must not be written onto this catalog row (None = OK)."""
    if not body.get("clinvar"):
        return "Sync skipped: no NCBI ClinVar record resolved for this catalog entry"
    rec = body["clinvar"][0]
    if meta.get("primary_complex"):
        return f"Sync skipped: ClinVar {rec.get('vcv_id') or rec['variation_id']} is a haplotype / diplotype record, not this single variant"
    if meta.get("how") not in ("rcv", "variation_id", "allele"):
        return (
            f"Sync skipped: no NCBI ClinVar record matches this catalog entry "
            f"({row.get('rcv_id')}, allele {row.get('risk_version') or row.get('alt_allele') or '?'}); "
            f"the record shown ({rec.get('vcv_id') or rec['variation_id']}) describes another allele"
        )
    yg = body.get("your_genotype") or {}
    if yg.get("carries_ncbi_alt") is False and yg.get("is_carried"):
        return (
            f"Sync skipped: your VCF call at this locus is a different ALT allele ({yg.get('alt')}) "
            f"than ClinVar {rec.get('vcv_id') or rec['variation_id']} ({rec.get('ref')}>{rec.get('alt')})"
        )
    return None


def catalog_record_ncbi(record_id: int, sync: bool = False) -> Optional[Dict[str, Any]]:
    row = _get_row(record_id)
    if row is None:
        return None
    meta: Dict[str, Any] = {}
    if is_structural_row(row):
        body = _structural_body(row)
    else:
        vid = str(row.get("variant_id") or "")
        rsid = vid if re.fullmatch(r"(?i)rs\d+", vid) else None
        rcv = row.get("rcv_id") if re.fullmatch(r"RCV\d{6,12}", str(row.get("rcv_id") or "")) else None
        vcv = row.get("ncbi_variation_id") or None
        chrom = row.get("chr")
        pos = row.get("position")
        # The ClinVar allele this catalog row describes (risk_version), not the gVCF block REF
        prefer_alt = row.get("risk_version") or row.get("alt_allele")
        prefer_rcvs = [rcv] if rcv else []
        prefer_rcvs += [r.strip() for r in str(row.get("all_rcv_ids") or "").split(",") if re.fullmatch(r"RCV\d{6,12}", r.strip()) and r.strip() not in prefer_rcvs]
        # Sync writes into the catalog: bypass the response cache so the data is current
        with get_client().budget(LOOKUP_BUDGET, fresh=bool(sync)):
            if rsid or rcv or vcv:
                body = variant_lookup(
                    rsid=rsid,
                    rcv_id=None if rsid else rcv,
                    vcv_id=vcv if (vcv and not rsid) else None,
                    include_publications=True,
                    prefer_alt=prefer_alt,
                    prefer_variation_id=vcv,
                    prefer_rcvs=prefer_rcvs,
                    meta=meta,
                )
                if vcv and rsid and not any(c["variation_id"] == str(vcv) for c in body.get("clinvar") or []):
                    body["errors"].append(f"Stored ClinVar variation {vcv} is no longer linked to {rsid}")
            else:
                body = variant_lookup(chrom=chrom, position=int(pos) if pos else None, include_publications=True, prefer_alt=prefer_alt, prefer_rcvs=prefer_rcvs, meta=meta)

    synced = False
    if sync:
        refusal = "Sync skipped: NCBI lookup is not applicable to structural / copy-number variants" if is_structural_row(row) else _sync_refusal(row, body, meta)
        if refusal:
            _add_error(body["errors"], refusal)
        else:
            rec = body["clinvar"][0]
            fields = _ncbi_fields_for_row(rec, body.get("risk_allele_frequency"), row.get("filter_status"), row.get("dp"), row.get("gq"), row.get("gene") or primary_gene_symbol(rec))
            if fields["condition"].endswith("-related genetic phenotype") and row.get("condition"):
                fields["condition"] = row["condition"]
            elif row.get("condition") and any(c["field"] == "condition" and c["match"] for c in _comparisons(row, body)):
                # Catalog text already names the NCBI condition(s), often with more traits
                # (other RCVs); keep it rather than narrowing it to the VCV-level trait list.
                fields["condition"] = row["condition"]
            _update_row(record_id, fields)
            row = _get_row(record_id) or row
            synced = True

    body["record_id"] = record_id
    body["catalog_vs_ncbi"] = _comparisons(row, body)
    body["synced"] = synced
    return body


def fallback_import_for_search(
    rcv_id: Optional[str] = None,
    variant_id: Optional[str] = None,
    gene: Optional[str] = None,
    vcv_id: Optional[str] = None,
    chrom: Optional[str] = None,
    position: Optional[str] = None,
    max_import: int = 20,
) -> int:
    """
    Used by GET /api/variants when an *identifier* (RCV, rsID, VCV, gene symbol or
    chr:position) is absent from the catalog and the MyVariant resolver found nothing.
    Never used for free text or condition searches. Single attempt per request and a
    short overall budget so a slow NCBI cannot stall catalog searches.
    """
    term: Optional[str] = None
    if rcv_id and re.fullmatch(r"(?i)\s*RCV\d{1,12}(?:\.\d+)?\s*", rcv_id):
        term = rcv_id.strip()
    elif variant_id and re.fullmatch(r"(?i)\s*rs\d{1,12}\s*", variant_id):
        term = variant_id.strip()
    elif vcv_id and re.fullmatch(r"(?i)\s*VCV\d{1,12}(?:\.\d+)?\s*", vcv_id):
        term = vcv_id.strip()
    elif gene and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9\-.]{0,29}", gene.strip()):
        term = gene.strip().upper()
    elif chrom and position and str(position).replace(",", "").strip().isdigit():
        term = f"{normalize_chrom(chrom)}:{str(position).replace(',', '').strip()}"
    if not term:
        return 0

    with get_client().budget(FALLBACK_BUDGET, attempts=1):
        errors: List[str] = []
        queries = translate_search_term(term)
        if gene and term == gene.strip().upper():
            queries = [("gene", gene_symbol_clause(gene))]
        queries = [q for q in queries if q[0] in ("rcv", "vcv", "rsid", "locus", "gene")]
        ids: List[str] = []
        for _kind, qq in queries:
            total, ids, _t = search_clinvar_ids(qq, 60, errors)
            if total > 0:
                break
        if not ids:
            return 0
        items = fetch_clinvar(ids, errors)
        items.sort(key=_item_order)
        chosen = [it[1]["variation_id"] for it in items if it[1].get("chr") and it[1].get("start")][:max_import]
        if not chosen:
            return 0
        return int(import_variation_ids(chosen, max_refsnp_fallbacks=3)["imported"])
