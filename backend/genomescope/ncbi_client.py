"""
NCBI client for GenomeScope WGS Explorer.

Layers
------
1. Pure parsing / normalisation helpers (no I/O) - unit-testable with fixture JSON:
   parse_esearch, parse_esummary_result, parse_clinvar_summary, parse_dbsnp_summary,
   parse_gene_summary, parse_pubmed_summary, parse_elink_ids, parse_global_maf,
   parse_spdi, spdi_to_vcf, alleles_equivalent, risk_allele_frequency_from_dbsnp,
   risk_allele_frequency_from_refsnp, translate_search_term, ...
2. Infrastructure: a global thread-safe RateLimiter and a persistent SQLite NcbiCache
   (table ncbi_cache in the active analysis' catalog DB).
3. NcbiClient: HTTP access (httpx) to E-utilities and Variation Services with
   rate limiting, caching, retry/backoff. All network traffic goes through
   NcbiClient._transport so tests can monkeypatch a single method.

Environment
-----------
NCBI_API_KEY          optional, raises the rate limit to 10 req/s and is sent as api_key
NCBI_EMAIL            optional, sent as email= when set (never hardcoded)
NCBI_CACHE_TTL_HOURS  cache TTL (default 168 = 7 days)
"""

from __future__ import annotations

import json
import os
import re
import contextlib
import threading
import time
import urllib.parse
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import httpx

from genomescope import indexer
from genomescope.alleles import alleles_equivalent, minimal_alleles  # noqa: F401  (re-exported)
from genomescope.vcf_engine import normalize_chrom
from genomescope.workspace import offline_mode

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
VARIATION_BASE = "https://api.ncbi.nlm.nih.gov/variation/v0/"
NCBI_WEB = "https://www.ncbi.nlm.nih.gov"
PUBMED_WEB = "https://pubmed.ncbi.nlm.nih.gov"
TOOL_NAME = "genomescope_wgs_explorer"

DEFAULT_TIMEOUT = 15.0
CONNECT_TIMEOUT = 3.0
MAX_ATTEMPTS = 3
RATE_SAFETY_FACTOR = 1.1
# Without an API key NCBI allows 3 req/s counted on *arrival*; TLS handshakes on new
# pooled connections can bunch arrivals, so keep a wider margin (~0.42 s spacing).
RATE_SAFETY_FACTOR_NO_KEY = 1.25
RATE_LIMIT_COOLDOWN = 1.0
DEFAULT_CACHE_TTL_HOURS = 168.0
# Empty esearch results ("not found") expire quickly so newly linked records show up.
NEGATIVE_CACHE_TTL_SECONDS = 3600.0
# Circuit breaker: after this many consecutive requests failed at the transport level,
# fail fast (serving stale cache where possible) for BREAKER_COOLDOWN_SECONDS.
BREAKER_THRESHOLD = 2
BREAKER_COOLDOWN_SECONDS = 60.0

# Preferred population-frequency studies (dbSNP esummary study names)
PREFERRED_FREQ_STUDIES = ["ALFA", "GnomAD_genomes", "TOPMED", "1000Genomes", "GnomAD_exomes", "1000Genomes_30X", "ExAC"]
# Variation Services study names -> esummary naming
REFSNP_STUDY_ALIASES = {"dbGaP_PopFreq": "ALFA"}

# RefSeq chromosome accessions (GRCh38) -> chromosome names
_NC_TO_CHROM = {f"NC_{i:06d}": str(i) for i in range(1, 23)}
_NC_TO_CHROM.update({"NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT"})


class NcbiError(Exception):
    """Raised for non-fatal upstream NCBI problems (callers convert to `errors` entries)."""


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# Pure helpers: E-utilities envelopes
# ---------------------------------------------------------------------------

def parse_esearch(payload: Dict[str, Any]) -> Dict[str, Any]:
    res = payload.get("esearchresult") if isinstance(payload, dict) else None
    if not isinstance(res, dict):
        raise NcbiError("Malformed esearch response")
    if res.get("ERROR"):
        raise NcbiError(f"esearch error: {res.get('ERROR')}")
    try:
        count = int(res.get("count") or 0)
    except (TypeError, ValueError):
        count = 0
    return {
        "count": count,
        "ids": [str(i) for i in (res.get("idlist") or [])],
        "translated": str(res.get("querytranslation") or ""),
    }


def parse_esummary_result(payload: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Return {uid: record} for records without an 'error' key."""
    res = payload.get("result") if isinstance(payload, dict) else None
    if not isinstance(res, dict):
        if isinstance(payload, dict) and payload.get("error"):
            raise NcbiError(f"esummary error: {payload.get('error')}")
        raise NcbiError("Malformed esummary response")
    out: Dict[str, Dict[str, Any]] = {}
    for uid in res.get("uids") or []:
        rec = res.get(str(uid))
        if isinstance(rec, dict) and not rec.get("error"):
            out[str(uid)] = rec
    return out


def parse_elink_ids(payload: Dict[str, Any], linkname: Optional[str] = None) -> List[str]:
    ids: List[str] = []
    for ls in (payload or {}).get("linksets") or []:
        for db in ls.get("linksetdbs") or []:
            if linkname and db.get("linkname") != linkname:
                continue
            for i in db.get("links") or []:
                s = str(i)
                if s not in ids:
                    ids.append(s)
    return ids


def parse_fasta_sequence(text: str) -> str:
    return "".join(line.strip() for line in (text or "").splitlines() if line and not line.startswith(">")).upper()


# ---------------------------------------------------------------------------
# Pure helpers: alleles, SPDI, VCF normalisation
# ---------------------------------------------------------------------------

def chrom_from_refseq(seq_id: str) -> Optional[str]:
    base = (seq_id or "").split(".")[0]
    return _NC_TO_CHROM.get(base)


def parse_spdi(spdi: Optional[str]) -> Optional[Tuple[str, int, str, str]]:
    """'NC_000012.12:111803961:G:A' -> ('NC_000012.12', 111803961, 'G', 'A'). Position is 0-based."""
    if not spdi or not isinstance(spdi, str):
        return None
    parts = spdi.strip().split(":")
    if len(parts) != 4:
        return None
    seq, pos, dele, ins = parts
    try:
        p = int(pos)
    except ValueError:
        return None
    if not re.fullmatch(r"[ACGTNacgtn]*", dele) or not re.fullmatch(r"[ACGTNacgtn]*", ins):
        return None
    return (seq, p, dele.upper(), ins.upper())


def spdi_to_vcf(
    pos0: int, deleted: str, inserted: str, anchor_base: Optional[str] = None
) -> Tuple[Optional[int], Optional[str], Optional[str], Optional[int]]:
    """
    Convert an SPDI allele to a left-aligned VCF (POS, REF, ALT).

    Returns (pos, ref, alt, anchor_needed_pos1). When the minimal representation is an
    indel (one side empty), VCF requires the reference base preceding the event; if
    `anchor_base` is not provided, anchor_needed_pos1 is the 1-based coordinate of the
    base that must be fetched and a best-effort fallback is returned:
    the un-normalised SPDI alleles when they share a leading base (valid VCF), else None.
    """
    if deleted == inserted:
        return (pos0 + 1, deleted or None, inserted or None, None)
    mpos, md, mi = minimal_alleles(pos0, deleted, inserted)
    if md and mi:
        return (mpos + 1, md, mi, None)
    anchor_pos1 = mpos  # 1-based coordinate of base at 0-based (mpos - 1)
    if anchor_base:
        a = anchor_base.upper()[:1]
        return (anchor_pos1, a + md, a + mi, None)
    if deleted and inserted and deleted[0] == inserted[0]:
        return (pos0 + 1, deleted, inserted, anchor_pos1)
    return (None, None, None, anchor_pos1)


def allele_names(pos0: int, deleted: str, inserted: str) -> Tuple[set, set]:
    """Names under which dbSNP may report the ref and alt allele of an SPDI (incl. '-' for empty)."""
    _p, md, mi = minimal_alleles(pos0, deleted, inserted)
    alt = {inserted or "-", mi or "-"}
    ref = {deleted or "-", md or "-"}
    return (ref - alt, alt - ref)


# ---------------------------------------------------------------------------
# Pure helpers: ClinVar
# ---------------------------------------------------------------------------

def review_status_to_stars(review_status: Optional[str]) -> int:
    rev = (review_status or "").lower()
    if "practice guideline" in rev:
        return 4
    if "expert panel" in rev:
        return 3
    if "multiple submitters" in rev and "no conflict" in rev:
        return 2
    if "conflicting" in rev or "single submitter" in rev:
        return 1
    return 0


def normalize_ncbi_date(raw: Optional[str]) -> Optional[str]:
    """'2021/03/24 00:00' -> '2021-03-24'; ClinVar's '1/01/01 00:00' placeholder -> None."""
    if not raw:
        return None
    m = re.match(r"^\s*(\d{4})[/-](\d{1,2})[/-](\d{1,2})", str(raw))
    if not m:
        return None
    y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
    if y < 1900:
        return None
    return f"{y:04d}-{mo:02d}-{d:02d}"


def classification_rank(normalized: str) -> int:
    s = (normalized or "").lower()
    if s.startswith("pathogenic"):
        return 0
    if s == "likely pathogenic":
        return 1
    if s in ("risk factor", "drug response") or "risk allele" in s:
        return 2
    if s in ("association", "protective"):
        return 3
    if s.startswith("conflicting"):
        return 4
    if s.startswith("uncertain"):
        return 5
    if s == "likely benign":
        return 6
    if s.startswith("benign"):
        return 7
    return 8


def normalize_clinvar_classification(raw_class: str) -> str:
    """ClinVar germline classification -> catalog label. An empty germline description
    (somatic-only or unclassified records) is 'Not Classified', never VUS."""
    if not str(raw_class or "").strip():
        return "Not Classified"
    return indexer.normalize_classification(raw_class)


def _classification_description(raw: Dict[str, Any], key: str) -> str:
    block = raw.get(key)
    return str(block.get("description") or "").strip() if isinstance(block, dict) else ""


def _grch38_location(vset: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    locs = vset.get("variation_loc") or []
    best = None
    for loc in locs:
        if str(loc.get("assembly_name", "")).upper().startswith("GRCH38"):
            if best is None or loc.get("status") == "current":
                best = loc
    return best


def _to_int(v: Any) -> Optional[int]:
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return None


def record_rsids(raw: Dict[str, Any]) -> List[str]:
    """Every dbSNP rsID referenced by any member allele of a ClinVar record (haplotypes and
    diplotypes list several)."""
    out: List[str] = []
    for vs in raw.get("variation_set") or []:
        for x in vs.get("variation_xrefs") or []:
            if x.get("db_source") == "dbSNP" and str(x.get("db_id") or "").isdigit():
                rs = f"rs{x['db_id']}"
                if rs not in out:
                    out.append(rs)
    return out


def extract_clinvar_spdi(raw: Dict[str, Any]) -> Optional[Tuple[str, int, str, str]]:
    for vset in raw.get("variation_set") or []:
        sp = parse_spdi(vset.get("canonical_spdi"))
        if sp:
            return sp
    return None


def parse_clinvar_summary(raw: Dict[str, Any]) -> Dict[str, Any]:
    """ClinVar esummary record -> NcbiClinVarRecord (ref/alt filled from SNV/MNV SPDI only;
    indels need an anchor base, see ncbi_service.fill_clinvar_vcf_alleles)."""
    uid = str(raw.get("uid") or "")
    gc = raw.get("germline_classification") or {}
    if not isinstance(gc, dict):
        gc = {}
    raw_class = str(gc.get("description") or "")
    if not raw_class and isinstance(raw.get("clinical_significance"), dict):
        raw_class = str(raw["clinical_significance"].get("description") or "")
    review_status = str(gc.get("review_status") or "")

    conditions: List[Dict[str, Any]] = []
    for t in gc.get("trait_set") or []:
        name = str(t.get("trait_name") or "").strip()
        if not name:
            continue
        medgen = None
        omim: List[str] = []
        for x in t.get("trait_xrefs") or []:
            src = str(x.get("db_source") or "")
            did = str(x.get("db_id") or "")
            if src == "MedGen" and did and medgen is None:
                medgen = did
            elif src == "OMIM" and did and did not in omim:
                omim.append(did)
        conditions.append(
            {
                "name": name,
                "medgen_id": medgen,
                "omim_ids": omim,
                "url": f"{NCBI_WEB}/medgen/{medgen}" if medgen else None,
            }
        )

    genes = [
        {"symbol": str(g.get("symbol") or ""), "gene_id": str(g.get("geneid") or "")}
        for g in raw.get("genes") or []
        if g.get("symbol")
    ]

    vsets = raw.get("variation_set") or []
    vset = vsets[0] if vsets else {}
    loc = None
    for vs in vsets:
        loc = _grch38_location(vs)
        if loc:
            vset = vs
            break
    chrom = normalize_chrom(loc.get("chr")) if loc and loc.get("chr") else None
    start = _to_int(loc.get("start")) if loc else None
    stop = _to_int(loc.get("stop")) if loc else None

    rsid = None
    for vs in vsets:
        for x in vs.get("variation_xrefs") or []:
            if x.get("db_source") == "dbSNP" and str(x.get("db_id") or "").isdigit():
                rsid = f"rs{x['db_id']}"
                break
        if rsid:
            break

    ref = alt = None
    sp = extract_clinvar_spdi(raw)
    if sp:
        pos, r, a, _need = spdi_to_vcf(sp[1], sp[2], sp[3])
        if _need is None and r and a:
            ref, alt = r, a
        if chrom is None:
            chrom = chrom_from_refseq(sp[0])
        if start is None and pos is not None:
            start = sp[1] + 1
            stop = sp[1] + max(1, len(sp[2]))

    accession = str(raw.get("accession") or "")
    protein = str(raw.get("protein_change") or "").strip()
    somatic_parts = []
    onco = _classification_description(raw, "oncogenicity_classification")
    impact = _classification_description(raw, "clinical_impact_classification")
    if onco:
        somatic_parts.append(f"Oncogenicity: {onco}")
    if impact:
        somatic_parts.append(f"Somatic clinical impact: {impact}")
    return {
        "variation_id": uid,
        "vcv_id": accession if accession.upper().startswith("VCV") else None,
        "title": str(raw.get("title") or ""),
        "variant_type": str(vset.get("variant_type") or raw.get("obj_type") or ""),
        "rcv_accessions": [str(r) for r in ((raw.get("supporting_submissions") or {}).get("rcv") or [])],
        "classification": raw_class,
        "classification_normalized": normalize_clinvar_classification(raw_class),
        "somatic_classification": " | ".join(somatic_parts) or None,
        "review_status": review_status,
        "review_stars": review_status_to_stars(review_status),
        "last_evaluated": normalize_ncbi_date(gc.get("last_evaluated")),
        "conditions": conditions,
        "genes": genes,
        "chr": chrom,
        "start": start,
        "stop": stop,
        "assembly": "GRCh38",
        "ref": ref,
        "alt": alt,
        "rsid": rsid,
        "molecular_consequences": [str(m) for m in raw.get("molecular_consequence_list") or []],
        "protein_change": protein or None,
        "url": f"{NCBI_WEB}/clinvar/variation/{uid}/",
    }


def clinvar_sort_key(rec: Dict[str, Any]) -> Tuple[int, int, int]:
    """Most severe classification first, then review stars, then the longest-established record (lowest id)."""
    return (classification_rank(rec.get("classification_normalized", "")), -int(rec.get("review_stars") or 0), _to_int(rec.get("variation_id")) or 0)


# ---------------------------------------------------------------------------
# Pure helpers: dbSNP
# ---------------------------------------------------------------------------

def parse_global_maf(entry: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """{"study":"GnomAD_genomes","freq":"A=0.0079314/1183"} -> {"study","allele","frequency","count"}"""
    study = str(entry.get("study") or "")
    raw = str(entry.get("freq") or "")
    m = re.match(r"^\s*([^=]*)=([0-9.eE+-]+)/(\d+)\s*$", raw)
    if not study or not m:
        return None
    try:
        freq = float(m.group(2))
    except ValueError:
        return None
    allele = m.group(1).strip() or "-"
    return {"study": study, "allele": allele, "frequency": freq, "count": int(m.group(3))}


def _freq_has_data(f: Dict[str, Any]) -> bool:
    return not (float(f.get("frequency") or 0) == 0.0 and int(f.get("count") or 0) == 0)


def pick_preferred_frequency(freqs: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    by_study = {f["study"]: f for f in freqs}
    for s in PREFERRED_FREQ_STUDIES:
        f = by_study.get(s)
        if f and _freq_has_data(f):
            return {"study": f["study"], "allele": f["allele"], "frequency": f["frequency"]}
    return None


def dbsnp_spdis(raw: Dict[str, Any]) -> List[Tuple[str, int, str, str]]:
    out = []
    for s in str(raw.get("spdi") or "").split(","):
        sp = parse_spdi(s)
        if sp:
            out.append(sp)
    return out


def parse_dbsnp_summary(raw: Dict[str, Any]) -> Dict[str, Any]:
    uid = str(raw.get("snp_id") or raw.get("uid") or "")
    sps = dbsnp_spdis(raw)
    chrom = normalize_chrom(raw.get("chr")) if raw.get("chr") else None
    position = None
    chrpos = str(raw.get("chrpos") or "")
    if ":" in chrpos:
        c, p = chrpos.split(":", 1)
        chrom = chrom or normalize_chrom(c)
        position = _to_int(p)
    ref = None
    alts: List[str] = []
    if sps:
        position = sps[0][1] + 1
        ref = sps[0][2] or "-"
        for sp in sps:
            a = sp[3] or "-"
            if a not in alts:
                alts.append(a)
        if chrom is None:
            chrom = chrom_from_refseq(sps[0][0])
    freqs = [f for f in (parse_global_maf(e) for e in raw.get("global_mafs") or []) if f]
    genes = [{"symbol": str(g.get("name") or ""), "gene_id": str(g.get("gene_id") or "")} for g in raw.get("genes") or [] if g.get("name")]
    return {
        "rsid": f"rs{uid}",
        "chr": chrom,
        "position": position,
        "ref": ref,
        "alts": alts,
        "spdi": str(raw.get("spdi")) if raw.get("spdi") else None,
        "variant_class": str(raw.get("snp_class") or "") or None,
        "function_class": [s for s in str(raw.get("fxn_class") or "").split(",") if s],
        "clinical_significance": [s for s in str(raw.get("clinical_significance") or "").split(",") if s],
        "genes": genes,
        "frequencies": freqs,
        "preferred_frequency": pick_preferred_frequency(freqs),
        "url": f"{NCBI_WEB}/snp/rs{uid}",
    }


def _allele_kind(deleted: str, inserted: str) -> str:
    _p, md, mi = minimal_alleles(0, deleted, inserted)
    if md and not mi:
        return "del"
    if mi and not md:
        return "ins"
    return "sub"


def risk_allele_frequency_from_dbsnp(
    dbsnp: Dict[str, Any],
    target_spdi: Optional[Tuple[str, int, str, str]],
    display_allele: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Frequency of the ClinVar ALT allele from dbSNP esummary global_mafs.
    global_mafs reports *a* minor allele per study, which may be the REF or another ALT:
    - reported allele == ClinVar ALT       -> use it (an explicit 0 is a real 0)
    - reported allele == REF, biallelic rs  -> ALT = 1 - f
    - otherwise                            -> skip that study (unknown)
    Indels: esummary names deletions '-' and insertions by their minimal bases, which is
    ambiguous when the rs has several deletion / insertion alleles; those short names are
    only trusted when the rs has exactly one allele of that kind (else use refsnp counts).
    """
    if not dbsnp or not target_spdi:
        return None
    _seq, p0, d, i = target_spdi
    ref_names, alt_names = allele_names(p0, d, i)
    # The ClinVar allele must be one of the rs alleles
    target_min = minimal_alleles(p0, d, i)[1:]
    rs_spdis = dbsnp_spdis({"spdi": dbsnp.get("spdi")})
    rs_alleles = [minimal_alleles(sp[1], sp[2], sp[3])[1:] for sp in rs_spdis]
    if rs_alleles and target_min not in rs_alleles:
        return None
    kind = _allele_kind(d, i)
    if kind in ("del", "ins"):
        same_kind = sum(1 for sp in rs_spdis if _allele_kind(sp[2], sp[3]) == kind)
        if same_kind != 1:
            # Only the full inserted sequence identifies the allele unambiguously
            alt_names = {i} if i else set()
            ref_names = set()
    biallelic = len(dbsnp.get("alts") or []) == 1
    disp = display_allele or (i or "-")
    by_study = {f["study"]: f for f in dbsnp.get("frequencies") or []}
    for study in PREFERRED_FREQ_STUDIES:
        f = by_study.get(study)
        if not f:
            continue
        allele = f["allele"]
        if allele in alt_names:
            return {"allele": disp, "frequency": float(f["frequency"]), "study": study}
        if not _freq_has_data(f):
            continue
        if allele in ref_names and biallelic:
            return {"allele": disp, "frequency": round(1.0 - float(f["frequency"]), 7), "study": study}
    return None


def parse_refsnp_frequencies(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Variation Services refsnp JSON -> [{study, deleted, inserted, allele_count, total_count}] for ALT observations."""
    out: List[Dict[str, Any]] = []
    ps = (payload or {}).get("primary_snapshot_data") or {}
    for aa in ps.get("allele_annotations") or []:
        for fr in aa.get("frequency") or []:
            obs = fr.get("observation") or {}
            d = str(obs.get("deleted_sequence") or "").upper()
            i = str(obs.get("inserted_sequence") or "").upper()
            if d == i:
                continue
            study = str(fr.get("study_name") or "")
            out.append(
                {
                    "study": REFSNP_STUDY_ALIASES.get(study, study),
                    "deleted": d,
                    "inserted": i,
                    "allele_count": int(fr.get("allele_count") or 0),
                    "total_count": int(fr.get("total_count") or 0),
                }
            )
    return out


def risk_allele_frequency_from_refsnp(
    refsnp_freqs: Sequence[Dict[str, Any]],
    target_spdi: Optional[Tuple[str, int, str, str]],
    display_allele: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """Exact per-allele counts (handles multi-allelic sites and indel naming)."""
    if not target_spdi:
        return None
    _seq, p0, d, i = target_spdi
    target = minimal_alleles(p0, d, i)[1:]
    for study in PREFERRED_FREQ_STUDIES:
        for fr in refsnp_freqs:
            if fr["study"] != study or fr["total_count"] <= 0:
                continue
            if minimal_alleles(0, fr["deleted"], fr["inserted"])[1:] == target:
                return {
                    "allele": display_allele or (i or "-"),
                    "frequency": round(fr["allele_count"] / fr["total_count"], 7),
                    "study": study,
                }
    return None


# ---------------------------------------------------------------------------
# Pure helpers: Gene / PubMed
# ---------------------------------------------------------------------------

def parse_gene_summary(raw: Dict[str, Any]) -> Dict[str, Any]:
    uid = str(raw.get("uid") or "")
    start = stop = None
    gi = raw.get("genomicinfo") or []
    if gi:
        a, b = _to_int(gi[0].get("chrstart")), _to_int(gi[0].get("chrstop"))
        if a is not None and b is not None:
            # esummary coordinates are 0-based; minus-strand genes have start > stop
            start, stop = min(a, b) + 1, max(a, b) + 1
    aliases = [s.strip() for s in str(raw.get("otheraliases") or "").split(",") if s.strip()]
    chrom = str(raw.get("chromosome") or "") or None
    return {
        "gene_id": uid,
        "symbol": str(raw.get("name") or ""),
        "name": str(raw.get("description") or ""),
        "summary": str(raw.get("summary") or ""),
        "aliases": aliases,
        "chr": normalize_chrom(chrom) if chrom else None,
        "start": start,
        "stop": stop,
        "map_location": str(raw.get("maplocation") or "") or None,
        "mim_ids": [str(m) for m in raw.get("mim") or []],
        "url": f"{NCBI_WEB}/gene/{uid}",
    }


def parse_pubmed_summary(raw: Dict[str, Any]) -> Dict[str, Any]:
    uid = str(raw.get("uid") or "")
    names = [str(a.get("name")) for a in raw.get("authors") or [] if a.get("name")]
    authors = ", ".join(names[:3]) + (" et al." if len(names) > 3 else "")
    pubdate = str(raw.get("pubdate") or raw.get("epubdate") or raw.get("sortpubdate") or "")
    m = re.search(r"(\d{4})", pubdate)
    return {
        "pmid": uid,
        "title": str(raw.get("title") or ""),
        "journal": str(raw.get("fulljournalname") or raw.get("source") or ""),
        "year": m.group(1) if m else "",
        "authors": authors,
        "url": f"{PUBMED_WEB}/{uid}/",
    }


# ---------------------------------------------------------------------------
# Pure helpers: search term translation / sanitisation
# ---------------------------------------------------------------------------

_SAFE_TEXT = re.compile(r"[^A-Za-z0-9 \-_.,'/+()>:]")
SIGNIFICANCE_FILTERS = {
    "pathogenic": '("clinsig pathogenic"[prop] OR "clinsig likely pathogenic"[prop])',
    "risk": '("clinsig risk factor"[prop] OR "clinsig other"[prop] OR "clinsig likely risk allele"[prop] OR "clinsig established risk allele"[prop])',
    "drug": '"clinsig drug response"[prop]',
    "all": "",
}


def sanitize_phrase(text: str) -> str:
    """Strip E-utilities syntax characters (quotes, brackets, #, &, ...) and collapse whitespace."""
    s = _SAFE_TEXT.sub(" ", str(text or ""))
    s = s.replace("(", " ").replace(")", " ")
    return re.sub(r"\s+", " ", s).strip()[:200]


def is_gene_symbol(term: str) -> bool:
    t = term.strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9\-]{1,14}", t):
        return False
    # Gene symbols are written upper-case or contain digits (BRCA1, HFE, CYP2C19, C9orf72)
    return t.upper() == t or any(c.isdigit() for c in t)


def gene_symbol_clause(symbol: str) -> str:
    return f'"{sanitize_phrase(symbol).upper()}"[gene]'


def translate_search_term(term: str) -> List[Tuple[str, str]]:
    """
    Map a user term to ClinVar esearch queries (ordered; later ones are fallbacks if
    earlier ones return nothing). Returns [(kind, query)].
    """
    t = (term or "").strip()
    m = re.fullmatch(r"(?i)RCV0*(\d{1,12})(?:\.\d+)?", t)
    if m:
        return [("rcv", f"RCV{int(m.group(1)):09d}")]
    m = re.fullmatch(r"(?i)VCV0*(\d{1,12})(?:\.\d+)?", t)
    if m:
        return [("vcv", f"VCV{int(m.group(1)):09d}[VACC]")]
    m = re.fullmatch(r"(?i)rs(\d{1,12})", t)
    if m:
        return [("rsid", f"rs{m.group(1)}[VRID]")]
    m = re.fullmatch(r"(?i)(?:chr)?([0-9]{1,2}|X|Y|MT|M):([0-9,]{1,12})", t)
    if m:
        c = normalize_chrom(m.group(1))
        p = int(m.group(2).replace(",", ""))
        return [("locus", f"{c}[chr] AND {p}[chrpos38]")]
    safe = sanitize_phrase(t)
    if not safe:
        return []
    words = [w for w in safe.split(" ") if w]
    quoted_words = " AND ".join(f'"{w}"[All Fields]' for w in words)
    if len(words) == 1 and is_gene_symbol(words[0]):
        return [("gene", gene_symbol_clause(words[0])), ("text", quoted_words)]
    return [("condition", f'"{safe}"[dis]'), ("text", quoted_words)]


# ---------------------------------------------------------------------------
# Infrastructure: rate limiter + persistent cache
# ---------------------------------------------------------------------------

class RateLimiter:
    """Global, thread-safe minimum-interval limiter (N requests per second)."""

    def __init__(self, rate_per_sec: float, safety_factor: float = RATE_SAFETY_FACTOR):
        self.rate_per_sec = max(0.1, float(rate_per_sec))
        # Safety margin: NCBI counts arrivals, and network jitter can bunch
        # requests that left exactly 1/rate apart into the same one-second window.
        self.interval = (1.0 / self.rate_per_sec) * max(1.0, float(safety_factor))
        self._lock = threading.Lock()
        self._next_slot = 0.0

    def acquire(self) -> float:
        with self._lock:
            now = time.monotonic()
            slot = max(now, self._next_slot)
            self._next_slot = slot + self.interval
        wait = slot - now
        if wait > 0:
            time.sleep(wait)
        return wait

    def cooldown(self, seconds: float) -> None:
        """Push every thread's next slot back (used after an HTTP 429 from NCBI)."""
        with self._lock:
            self._next_slot = max(self._next_slot, time.monotonic() + max(0.0, seconds))


class NcbiCache:
    """Persistent cache in table ncbi_cache(key TEXT PRIMARY KEY, payload TEXT, fetched_at REAL)."""

    def __init__(self, ttl_hours: float):
        self.ttl_seconds = float(ttl_hours) * 3600.0
        self._ready_paths: set = set()
        self._lock = threading.Lock()

    def _conn(self):
        conn = indexer.get_db_connection()
        path = str(indexer.db_path())
        if path not in self._ready_paths:
            with self._lock:
                conn.execute("CREATE TABLE IF NOT EXISTS ncbi_cache (key TEXT PRIMARY KEY, payload TEXT, fetched_at REAL)")
                conn.commit()
                self._ready_paths.add(path)
        return conn

    def get_entry(self, key: str) -> Optional[Tuple[str, float]]:
        """(payload, fetched_at epoch seconds) regardless of age, or None."""
        try:
            conn = self._conn()
            try:
                row = conn.execute("SELECT payload, fetched_at FROM ncbi_cache WHERE key = ?", (key,)).fetchone()
            finally:
                conn.close()
        except Exception:
            return None
        if not row:
            return None
        return (row[0], float(row[1] or 0))

    def get(self, key: str, allow_stale: bool = False) -> Optional[str]:
        entry = self.get_entry(key)
        if not entry:
            return None
        if not allow_stale and time.time() - entry[1] > self.ttl_seconds:
            return None
        return entry[0]

    def set(self, key: str, payload: str) -> None:
        try:
            conn = self._conn()
            try:
                conn.execute("INSERT OR REPLACE INTO ncbi_cache (key, payload, fetched_at) VALUES (?, ?, ?)", (key, payload, time.time()))
                conn.commit()
            finally:
                conn.close()
        except Exception:
            pass

    def count(self) -> int:
        try:
            conn = self._conn()
            try:
                return int(conn.execute("SELECT COUNT(*) FROM ncbi_cache").fetchone()[0])
            finally:
                conn.close()
        except Exception:
            return 0

    def clear(self) -> None:
        try:
            conn = self._conn()
            try:
                conn.execute("DELETE FROM ncbi_cache")
                conn.commit()
            finally:
                conn.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# HTTP client
# ---------------------------------------------------------------------------

def _is_negative_esearch(url: str, payload: str) -> bool:
    if not url.endswith("esearch.fcgi"):
        return False
    try:
        return int(((json.loads(payload) or {}).get("esearchresult") or {}).get("count") or 0) == 0
    except (ValueError, TypeError, AttributeError):
        return False


class NcbiClient:
    def __init__(self) -> None:
        self.api_key = (os.environ.get("NCBI_API_KEY") or "").strip() or None
        self.email = (os.environ.get("NCBI_EMAIL") or "").strip() or None
        try:
            ttl = float(os.environ.get("NCBI_CACHE_TTL_HOURS") or DEFAULT_CACHE_TTL_HOURS)
        except ValueError:
            ttl = DEFAULT_CACHE_TTL_HOURS
        self.cache_ttl_hours = ttl
        self.rate_limit_per_sec = 10 if self.api_key else 3
        self.limiter = RateLimiter(self.rate_limit_per_sec, RATE_SAFETY_FACTOR if self.api_key else RATE_SAFETY_FACTOR_NO_KEY)
        self.cache = NcbiCache(ttl)
        self.timeout = DEFAULT_TIMEOUT
        self.connect_timeout = CONNECT_TIMEOUT
        self.backoff_base = 0.5
        self.breaker_threshold = BREAKER_THRESHOLD
        self.breaker_cooldown = BREAKER_COOLDOWN_SECONDS
        self._breaker_lock = threading.Lock()
        self._consecutive_failures = 0
        self._breaker_open_until = 0.0
        self._ctx = threading.local()
        self._http: Optional[httpx.Client] = None
        self._http_lock = threading.Lock()

    # -- per-call context: deadline budget, attempts cap, cache bypass -----------
    @contextlib.contextmanager
    def budget(self, seconds: Optional[float] = None, attempts: Optional[int] = None, fresh: Optional[bool] = None):
        """
        Bound every NCBI request made by this thread inside the block:
        - seconds: overall wall-clock deadline (nested budgets never extend an outer one)
        - attempts: maximum attempts per request (e.g. 1 for interactive fallbacks)
        - fresh: bypass cached responses (still refreshes the cache), e.g. for catalog sync
        """
        prev = (getattr(self._ctx, "deadline", None), getattr(self._ctx, "attempts", None), getattr(self._ctx, "fresh", False))
        deadline = prev[0]
        if seconds is not None:
            d = time.monotonic() + max(0.0, float(seconds))
            deadline = d if deadline is None else min(deadline, d)
        self._ctx.deadline = deadline
        if attempts is not None:
            self._ctx.attempts = attempts if prev[1] is None else min(prev[1], attempts)
        if fresh is not None:
            self._ctx.fresh = bool(fresh) or prev[2]
        try:
            yield self
        finally:
            self._ctx.deadline, self._ctx.attempts, self._ctx.fresh = prev

    def _remaining(self) -> Optional[float]:
        d = getattr(self._ctx, "deadline", None)
        return None if d is None else d - time.monotonic()

    def breaker_open(self) -> bool:
        with self._breaker_lock:
            return time.monotonic() < self._breaker_open_until

    def _record_result(self, transport_ok: bool) -> None:
        with self._breaker_lock:
            if transport_ok:
                self._consecutive_failures = 0
                self._breaker_open_until = 0.0
            else:
                self._consecutive_failures += 1
                if self._consecutive_failures >= self.breaker_threshold:
                    self._breaker_open_until = time.monotonic() + self.breaker_cooldown

    # -- transport (monkeypatch target in tests) --------------------------------
    def _transport(self, url: str, params: Dict[str, Any], timeout: float) -> Tuple[int, str, Dict[str, str]]:
        if offline_mode():
            raise httpx.ConnectError("offline mode (GENOMESCOPE_OFFLINE=1): NCBI is not contacted")
        with self._http_lock:
            if self._http is None:
                self._http = httpx.Client(
                    timeout=httpx.Timeout(DEFAULT_TIMEOUT, connect=CONNECT_TIMEOUT),
                    headers={"User-Agent": f"{TOOL_NAME}/1.0 (+FastAPI; NCBI E-utilities client)"},
                    follow_redirects=True,
                )
            client = self._http
        # Short connect timeout (fail fast when NCBI is unreachable), longer read timeout
        resp = client.get(url, params=params, timeout=httpx.Timeout(timeout, connect=min(timeout, self.connect_timeout)))
        return (resp.status_code, resp.text, dict(resp.headers))

    @staticmethod
    def cache_key(url: str, params: Dict[str, Any]) -> str:
        clean = sorted((k, str(v)) for k, v in params.items() if k not in ("api_key", "email", "tool"))
        return url + "?" + urllib.parse.urlencode(clean)

    def _with_identity(self, url: str, params: Dict[str, Any]) -> Dict[str, Any]:
        p = dict(params)
        if url.startswith(EUTILS_BASE):
            p["tool"] = TOOL_NAME
            if self.email:
                p["email"] = self.email
            if self.api_key:
                p["api_key"] = self.api_key
        return p

    @staticmethod
    def _decode(text: str, kind: str) -> Any:
        return json.loads(text) if kind == "json" else text

    def _stale_or_raise(self, key: str, kind: str, use_cache: bool, message: str) -> Any:
        """Network trouble: serve an expired cache entry when one exists (better than nothing)."""
        if use_cache:
            stale = self.cache.get(key, allow_stale=True)
            if stale is not None:
                try:
                    return self._decode(stale, kind)
                except ValueError:
                    pass
        raise NcbiError(message)

    def request(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        kind: str = "json",
        use_cache: bool = True,
        attempts: int = MAX_ATTEMPTS,
        timeout: Optional[float] = None,
        probe: bool = False,
    ) -> Any:
        params = params or {}
        key = self.cache_key(url, params)
        short = url.replace(EUTILS_BASE, "eutils/")
        fresh = bool(getattr(self._ctx, "fresh", False))
        if use_cache and not fresh:
            entry = self.cache.get_entry(key)
            if entry is not None:
                ttl = self.cache.ttl_seconds
                if _is_negative_esearch(url, entry[0]):
                    ttl = min(ttl, NEGATIVE_CACHE_TTL_SECONDS)
                if time.time() - entry[1] <= ttl:
                    try:
                        return self._decode(entry[0], kind)
                    except ValueError:
                        pass
        if not probe and self.breaker_open():
            return self._stale_or_raise(key, kind, use_cache, f"NCBI unreachable (recent connection failures; retrying shortly) for {short}")
        ctx_attempts = getattr(self._ctx, "attempts", None)
        if ctx_attempts is not None:
            attempts = min(attempts, ctx_attempts)
        send = self._with_identity(url, params)
        last_err = "unknown error"
        transport_ok = False
        attempted = 0
        for attempt in range(max(1, attempts)):
            remaining = self._remaining()
            if remaining is not None and remaining <= 0.05:
                last_err = "time budget exhausted"
                break
            self.limiter.acquire()
            req_timeout = timeout or self.timeout
            remaining = self._remaining()
            if remaining is not None:
                if remaining <= 0.05:
                    last_err = "time budget exhausted"
                    break
                req_timeout = min(req_timeout, remaining)
            attempted += 1
            try:
                status, text, headers = self._transport(url, send, req_timeout)
                transport_ok = True
            except httpx.HTTPError as exc:
                last_err = f"{type(exc).__name__}: {exc}"
                status, text, headers = 0, "", {}
            except Exception as exc:  # pragma: no cover - defensive
                last_err = f"{type(exc).__name__}: {exc}"
                status, text, headers = 0, "", {}
            headers = headers or {}
            quota_left = headers.get("x-ratelimit-remaining") or headers.get("X-RateLimit-Remaining")
            if quota_left is not None and str(quota_left).strip() == "0":
                # NCBI says this second's quota is used up: hold every thread back briefly
                self.limiter.cooldown(RATE_LIMIT_COOLDOWN)
            if status == 200:
                self._record_result(True)
                if kind == "json":
                    try:
                        data = json.loads(text)
                    except ValueError:
                        last_err = "invalid JSON from NCBI"
                        break
                    if isinstance(data, dict) and isinstance(data.get("error"), str) and len(data) <= 2:
                        raise NcbiError(f"NCBI error: {data['error'][:200]}")
                    if use_cache:
                        self.cache.set(key, text)
                    return data
                if use_cache:
                    self.cache.set(key, text)
                return text
            if status and status not in (429,) and status < 500:
                self._record_result(True)
                raise NcbiError(f"NCBI HTTP {status} for {short}")
            if status:
                last_err = f"HTTP {status}"
            if status == 429:
                self.limiter.cooldown(RATE_LIMIT_COOLDOWN)
            if attempt + 1 < attempts:
                delay = self.backoff_base * (2 ** attempt)
                ra = headers.get("retry-after") or headers.get("Retry-After")
                if ra and str(ra).isdigit():
                    delay = max(delay, min(10.0, float(ra)))
                remaining = self._remaining()
                if remaining is not None:
                    delay = min(delay, max(0.0, remaining))
                time.sleep(delay)
        if attempted and not transport_ok:
            self._record_result(False)  # only real connection failures trip the breaker
        return self._stale_or_raise(key, kind, use_cache, f"NCBI request failed after {attempts} attempt(s) ({last_err}) for {short}")

    # -- E-utilities ------------------------------------------------------------
    def esearch(self, db: str, term: str, retmax: int = 20, retstart: int = 0, sort: Optional[str] = None) -> Dict[str, Any]:
        params: Dict[str, Any] = {"db": db, "term": term, "retmode": "json", "retmax": int(retmax), "retstart": int(retstart)}
        if sort:
            params["sort"] = sort
        return parse_esearch(self.request(EUTILS_BASE + "esearch.fcgi", params))

    def esummary(self, db: str, ids: Iterable[str], chunk: int = 100) -> Dict[str, Dict[str, Any]]:
        uniq: List[str] = []
        for i in ids:
            s = str(i).strip()
            if s and s not in uniq:
                uniq.append(s)
        out: Dict[str, Dict[str, Any]] = {}
        for k in range(0, len(uniq), chunk):
            part = uniq[k:k + chunk]
            payload = self.request(EUTILS_BASE + "esummary.fcgi", {"db": db, "id": ",".join(part), "retmode": "json"})
            out.update(parse_esummary_result(payload))
        return out

    def elink(self, dbfrom: str, db: str, ids: Iterable[str], linkname: Optional[str] = None) -> List[str]:
        params: Dict[str, Any] = {"dbfrom": dbfrom, "db": db, "id": ",".join(str(i) for i in ids), "retmode": "json"}
        if linkname:
            params["linkname"] = linkname
        return parse_elink_ids(self.request(EUTILS_BASE + "elink.fcgi", params), linkname)

    def fetch_sequence(self, accession: str, start1: int, stop1: int) -> str:
        text = self.request(
            EUTILS_BASE + "efetch.fcgi",
            {"db": "nuccore", "id": accession, "seq_start": int(start1), "seq_stop": int(stop1), "rettype": "fasta", "retmode": "text"},
            kind="text",
        )
        return parse_fasta_sequence(text)

    def refsnp(self, rs_number: str) -> Dict[str, Any]:
        return self.request(VARIATION_BASE + f"refsnp/{int(rs_number)}")

    def ping(self) -> Tuple[bool, Optional[str]]:
        try:
            self.request(
                EUTILS_BASE + "esearch.fcgi",
                {"db": "clinvar", "term": "17865[VID]", "retmode": "json", "retmax": 0},
                use_cache=False,
                attempts=1,
                timeout=8.0,
                probe=True,  # bypasses (and, on success, closes) the circuit breaker
            )
            return (True, None)
        except NcbiError as exc:
            return (False, str(exc))


_client_lock = threading.Lock()
_client: Optional[NcbiClient] = None


def get_client() -> NcbiClient:
    global _client
    with _client_lock:
        if _client is None:
            _client = NcbiClient()
        return _client


def set_client(client: Optional[NcbiClient]) -> None:
    """Replace the global client (tests)."""
    global _client
    with _client_lock:
        _client = client
