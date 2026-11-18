import json
import sqlite3
import ssl
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from genomescope import workspace
from genomescope.alleles import DIFFERENT_ALT_STATUS, UNVERIFIED_ALLELE_STATUS, is_plain_allele
from genomescope.vcf_engine import (
    CNV_VCF_NAME,
    SNP_VCF_NAME,
    SV_VCF_NAME,
    determine_variant_type,
    format_vcf_genotype_and_status,
    normalize_chrom,
    parse_info_field,
    parse_sample_field,
    vcf_engine,
)

def _tls_context() -> ssl.SSLContext:
    """Verified TLS for MyVariant.info / mygene.info (certifi bundle when available)."""
    try:
        import certifi  # shipped with httpx

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:  # pragma: no cover - certifi missing
        return ssl.create_default_context()


SSL_CTX = _tls_context()
HTTP_HEADERS = {"User-Agent": "GenomeScope/2.0 (variant annotation client)"}


def db_path() -> Path:
    """SQLite catalog of the active analysis workspace."""
    return workspace.current().db_path


def get_db_connection() -> sqlite3.Connection:
    path = db_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path), check_same_thread=False, timeout=30)
    conn.row_factory = sqlite3.Row
    return conn


def _http_get_json(url: str, timeout: float) -> Any:
    req = urllib.request.Request(url, headers=HTTP_HEADERS)
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as resp:
        return json.loads(resp.read().decode("utf-8", errors="replace"))


def init_db() -> None:
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            unique_key TEXT UNIQUE,
            rcv_id TEXT,
            all_rcv_ids TEXT,
            chr TEXT,
            position INTEGER,
            end_position INTEGER,
            variant_id TEXT,
            gene TEXT,
            gene_name TEXT,
            ref_allele TEXT,
            alt_allele TEXT,
            risk_version TEXT,
            risk_frequency REAL,
            risk_frequency_text TEXT,
            variant_type TEXT,
            your_data TEXT,
            your_status TEXT,
            genotype TEXT,
            alleles_called TEXT,
            condition TEXT,
            classification TEXT,
            evidence TEXT,
            evidence_stars INTEGER,
            review_status TEXT,
            hgvs TEXT,
            molecular_consequence TEXT,
            dp INTEGER,
            gq INTEGER,
            ad TEXT,
            qual REAL,
            filter_status TEXT,
            is_carried INTEGER,
            source_file TEXT,
            raw_vcf_line TEXT
        )
        """
    )
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_rcv ON variants(rcv_id)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_chr_pos ON variants(chr, position)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_vid ON variants(variant_id)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_type ON variants(variant_type)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_status ON variants(your_status)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_class ON variants(classification)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_carried ON variants(is_carried)")
    conn.commit()
    conn.close()
    ensure_ncbi_columns(force=True)


# Columns added for NCBI (ClinVar / dbSNP) provenance. Added via guarded ALTER TABLE so
# existing databases are migrated in place.
NCBI_EXTRA_COLUMNS: List[Tuple[str, str]] = [
    ("ncbi_variation_id", "TEXT"),
    ("ncbi_last_evaluated", "TEXT"),
    ("ncbi_synced_at", "TEXT"),
    ("frequency_source", "TEXT"),
    # Local ClinVar release scan (genomescope/clinvar_scan.py): which ClinVar variation a row was
    # verified against, the release it was verified with, and how the allele was matched
    # (exact | equivalent | adjudicated_carried | adjudicated_not_carried | ... | unverified).
    ("clinvar_variation_id", "TEXT"),
    ("clinvar_release", "TEXT"),
    ("clinvar_match", "TEXT"),
]
_migrated_db_paths: set = set()


def ensure_ncbi_columns(force: bool = False) -> None:
    path = str(db_path())
    if path in _migrated_db_paths and not force:
        return
    conn = get_db_connection()
    try:
        cur = conn.cursor()
        cur.execute("PRAGMA table_info(variants)")
        existing = {r[1] for r in cur.fetchall()}
        if not existing:
            return  # table not created yet; init_db() will call us again
        for name, ctype in NCBI_EXTRA_COLUMNS:
            if name not in existing:
                cur.execute(f"ALTER TABLE variants ADD COLUMN {name} {ctype}")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_ncbi_vid ON variants(ncbi_variation_id)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_variants_clinvar_vid ON variants(clinvar_variation_id)")
        conn.commit()
        _migrated_db_paths.add(path)
    finally:
        conn.close()


def normalize_classification(raw_sig: str) -> str:
    s = (raw_sig or "").lower().strip()

    def suffixed(base: str) -> str:
        # Keep the secondary ClinVar term (e.g. "...; risk factor") visible on the label
        if "drug response" in s:
            return f"{base} / Drug Response"
        if "risk factor" in s:
            return f"{base} / Risk Factor"
        return base

    # Aggregate "Pathogenic/Likely pathogenic" is at least as severe as Likely Pathogenic:
    # handle it before the generic "likely pathogenic" branch so it is not downgraded.
    if "pathogenic/likely pathogenic" in s:
        return suffixed("Pathogenic/Likely Pathogenic")
    if "pathogenic" in s and "likely pathogenic" not in s and "conflicting" not in s:
        return suffixed("Pathogenic")
    # Conflicts must stay visible even when a secondary "risk factor" / "drug response" term exists
    if "conflicting" in s:
        return suffixed("Conflicting Interpretations")
    if "likely pathogenic" in s:
        return "Likely Pathogenic"
    if "risk factor" in s:
        return "Risk Factor"
    if "drug response" in s:
        return "Drug Response"
    if "protective" in s:
        return "Protective"
    if "association" in s or "affects" in s:
        return "Association"
    if "uncertain" in s:
        return "Uncertain Significance"
    if "benign/likely benign" in s:
        return "Benign/Likely Benign"
    if "likely benign" in s:
        return "Likely Benign"
    if "benign" in s:
        return "Benign"
    return raw_sig.title() if raw_sig else "Uncertain Significance"


def map_review_status_to_evidence(review_status: str, filt: str, dp: Optional[int], gq: Optional[int]) -> Tuple[str, int]:
    rev = (review_status or "").lower()
    if "practice guideline" in rev:
        stars = 4
        label = "Practice Guideline (4 Stars)"
    elif "expert panel" in rev:
        stars = 3
        label = "Expert Panel Reviewed (3 Stars)"
    elif "multiple submitters" in rev and "no conflicts" in rev:
        stars = 2
        label = "Strong - Multiple Submitters (2 Stars)"
    elif "conflicting" in rev:
        stars = 1
        label = "Conflicting Interpretations (1 Star)"
    elif "single submitter" in rev:
        stars = 1
        label = "Moderate - Single Submitter (1 Star)"
    else:
        stars = 0
        label = "Limited / No Assertion Criteria"

    read_parts = [f"VCF {filt}"]
    if dp is not None:
        read_parts.append(f"DP {dp}")
    if gq is not None and gq > 0:
        read_parts.append(f"GQ {gq}")
    return (f"{label} | {', '.join(read_parts)}", stars)


def format_frequency(af: Optional[float]) -> Tuple[Optional[float], str]:
    if af is None:
        return (None, "Rare / Not Reported")
    try:
        val = float(af)
        pct = val * 100.0
        if pct >= 1.0:
            return (val, f"{pct:.2f}% ({val:.4f})")
        elif pct >= 0.01:
            return (val, f"{pct:.3f}% ({val:.5f})")
        else:
            return (val, f"{pct:.4f}% ({val:.6f})")
    except Exception:
        return (None, "Rare / Not Reported")


def extract_clinvar_records_from_hit(hit: Dict[str, Any]) -> List[Dict[str, Any]]:
    cv = hit.get("clinvar")
    if not isinstance(cv, dict):
        return []
    chrom = normalize_chrom(str(cv.get("chrom", "")))
    hg38 = cv.get("hg38")
    if not isinstance(hg38, dict) or not hg38.get("start"):
        return []
    pos = int(hg38["start"])
    end_pos = int(hg38.get("end", pos))

    rsid_val = cv.get("rsid") or hit.get("dbsnp", {}).get("rsid")
    if isinstance(rsid_val, list):
        rsid_val = rsid_val[0]
    rsid = str(rsid_val) if rsid_val else f"chr{chrom}:{pos}"
    if rsid.isdigit():
        rsid = f"rs{rsid}"

    gene_obj = cv.get("gene", {})
    gene_symbol = "Intergenic"
    if isinstance(gene_obj, dict):
        gene_symbol = gene_obj.get("symbol") or "Intergenic"
    elif isinstance(gene_obj, list) and gene_obj:
        gene_symbol = gene_obj[0].get("symbol") or "Intergenic"

    ref = str(cv.get("ref", "N"))
    alt = str(cv.get("alt", "N"))

    # Extract allele frequency
    af = None
    gnomad = hit.get("gnomad_genome", {})
    if isinstance(gnomad, dict):
        af_obj = gnomad.get("af", {})
        if isinstance(af_obj, dict):
            af = af_obj.get("af")
        elif isinstance(af_obj, (int, float)):
            af = af_obj
    if af is None:
        exac = hit.get("exac", {})
        if isinstance(exac, dict):
            af = exac.get("af")

    hgvs_obj = cv.get("hgvs", {})
    hgvs_str = ""
    if isinstance(hgvs_obj, dict):
        coding = hgvs_obj.get("coding")
        protein = hgvs_obj.get("protein")
        c_str = coding[0] if isinstance(coding, list) and coding else (coding if isinstance(coding, str) else "")
        p_str = protein[0] if isinstance(protein, list) and protein else (protein if isinstance(protein, str) else "")
        hgvs_str = f"{c_str} ({p_str})" if (c_str and p_str) else (c_str or p_str)

    rcvs = cv.get("rcv", [])
    if isinstance(rcvs, dict):
        rcvs = [rcvs]
    elif not isinstance(rcvs, list):
        rcvs = []

    all_rcv_list: List[str] = []
    valid_rcvs: List[Tuple[str, str, str, str]] = []
    for r in rcvs:
        if not isinstance(r, dict):
            continue
        acc = r.get("accession", "")
        if acc:
            all_rcv_list.append(acc)
        sig = r.get("clinical_significance", "Uncertain significance")
        rev = r.get("review_status", "")
        cond_obj = r.get("conditions", {})
        cond_names: List[str] = []
        if isinstance(cond_obj, dict):
            n = cond_obj.get("name")
            if n:
                cond_names.append(str(n))
        elif isinstance(cond_obj, list):
            for co in cond_obj:
                if isinstance(co, dict) and co.get("name"):
                    cond_names.append(str(co["name"]))
        cond_str = "; ".join(cond_names) if cond_names else "Genetic Condition"
        valid_rcvs.append((acc, sig, rev, cond_str))

    if not valid_rcvs:
        return []

    def rev_score(rev_str: str) -> int:
        r_low = (rev_str or "").lower()
        if "practice guideline" in r_low:
            return 4
        if "expert panel" in r_low:
            return 3
        if "multiple submitters" in r_low and "no conflicts" in r_low:
            return 2
        if "single submitter" in r_low or "conflicting" in r_low:
            return 1
        return 0

    def rcv_priority(item: Tuple[str, str, str, str]) -> Tuple[int, int, int]:
        _acc, sig, rev, cond = item
        c_low = cond.lower()
        has_named_cond = 0 if c_low in ("not provided", "not specified", "none provided", "see cases", "genetic condition") else 1
        s_low = sig.lower()
        sig_score = 0
        if "pathogenic" in s_low:
            sig_score = 5
        elif "risk factor" in s_low or "drug response" in s_low:
            sig_score = 4
        elif "association" in s_low or "protective" in s_low:
            sig_score = 3
        elif "uncertain" in s_low:
            sig_score = 2
        else:
            sig_score = 1
        return (has_named_cond, sig_score, rev_score(rev))

    valid_rcvs.sort(key=rcv_priority, reverse=True)
    best_acc, best_sig, best_rev, best_cond = valid_rcvs[0]
    # Use the highest-ranked review_status across the variant's RCV submissions
    best_rev_overall = max((r[2] for r in valid_rcvs), key=rev_score, default=best_rev)
    if rev_score(best_rev_overall) > rev_score(best_rev):
        best_rev = best_rev_overall

    # Combine unique named conditions across RCVs if available
    named_conditions: List[str] = []
    for _acc, _sig, _rev, cond in valid_rcvs:
        if cond.lower() not in ("not provided", "not specified", "none provided", "see cases", "genetic condition"):
            for part in cond.split("; "):
                p_clean = part.strip()
                if p_clean and p_clean.lower() not in ("not provided", "not specified", "none provided", "see cases") and p_clean not in named_conditions:
                    named_conditions.append(p_clean)
    final_condition = " | ".join(named_conditions[:3]) if named_conditions else best_cond
    if final_condition.lower() in ("not provided", "not specified"):
        final_condition = f"{gene_symbol}-related genetic phenotype"

    return [
        {
            "rcv_id": best_acc or "RCV000000000",
            "all_rcv_ids": ", ".join(dict.fromkeys(all_rcv_list)),
            "chr": chrom,
            "position": pos,
            "end_position": end_pos,
            "variant_id": rsid,
            "gene": gene_symbol,
            "ref_allele": ref,
            "alt_allele": alt,
            "af": af,
            "classification": normalize_classification(best_sig),
            "review_status": best_rev,
            "condition": final_condition,
            "hgvs": hgvs_str or str(cv.get("rcv", [{}])[0].get("preferred_name", "")) if isinstance(cv.get("rcv"), list) else hgvs_str,
            "clinvar_type": str(cv.get("type", "single nucleotide variant")),
        }
    ]


def legacy_allele_verdict(cv_rec: Dict[str, Any], vcf_parts: List[str]) -> Tuple[str, Optional[int]]:
    """
    Is the sample's call at this position the ClinVar allele of `cv_rec`?
    Returns (verdict, alt_index): "ref_block" | "match" | "different" | "unverified".

    The MyVariant-derived record joins to the gVCF by position only, so the alleles must be
    compared before anything is called "carried". Only a substitution whose ClinVar allele is
    one of the row's ALTs is confirmed here; indels and unusable alleles stay "unverified"
    (authoritative, representation-aware matching is done by genomescope.clinvar_scan).
    """
    alt = vcf_parts[4]
    if alt == ".":
        return ("ref_block", None)
    alts = [a.upper() for a in alt.split(",")]
    cref, calt = str(cv_rec.get("ref_allele", "")).upper(), str(cv_rec.get("alt_allele", "")).upper()
    if not (is_plain_allele(cref) and is_plain_allele(calt)) or cref == calt:
        return ("unverified", None)
    if len(cref) == 1 and len(calt) == 1 and int(vcf_parts[1]) == int(cv_rec["position"]):
        if calt in alts:
            return ("match", alts.index(calt) + 1)
        return ("different", None)
    return ("unverified", None)


def build_record_from_clinvar_and_vcf(
    cv_rec: Dict[str, Any], vcf_parts: Optional[List[str]]
) -> Optional[Dict[str, Any]]:
    if not vcf_parts or len(vcf_parts) < 10:
        return None

    vcf_ref = vcf_parts[3]
    vcf_alt = vcf_parts[4]
    vcf_qual = vcf_parts[5]
    vcf_filt = vcf_parts[6]
    info_map = parse_info_field(vcf_parts[7])
    sample_map = parse_sample_field(vcf_parts[8], vcf_parts[9])

    fmt = format_vcf_genotype_and_status(
        vcf_ref,
        vcf_alt,
        vcf_qual,
        vcf_filt,
        info_map,
        sample_map,
        source_type="snp-indel",
        expected_risk_allele=cv_rec["alt_allele"],
    )

    vtype = determine_variant_type(
        cv_rec["ref_allele"],
        cv_rec["alt_allele"],
        info_map,
        source_type="snp-indel",
    )

    freq_val, freq_text = format_frequency(cv_rec["af"])
    evidence_str, stars = map_review_status_to_evidence(
        cv_rec["review_status"], vcf_filt, fmt["dp"], fmt["gq"]
    )

    vid = cv_rec["variant_id"]
    if vcf_parts[2] != "." and vcf_parts[2].startswith("rs"):
        vid = vcf_parts[2]

    # Never call a variant "carried" unless the ClinVar allele itself was found in the call
    verdict, alt_k = legacy_allele_verdict(cv_rec, vcf_parts)
    is_carried = bool(fmt["is_carried"])
    your_status = fmt["your_status"]
    if verdict == "match":
        gt_idx = [x for x in sample_map.get("GT", "./.").replace("|", "/").split("/")]
        if str(alt_k) not in gt_idx:
            verdict = "different"
    if verdict == "different":
        is_carried, your_status = False, DIFFERENT_ALT_STATUS
    elif verdict == "unverified" and is_carried:
        is_carried, your_status = False, UNVERIFIED_ALLELE_STATUS

    unique_key = f"SNP:{cv_rec['chr']}:{cv_rec['position']}:{vid}:{cv_rec['rcv_id']}"

    return {
        "unique_key": unique_key,
        "rcv_id": cv_rec["rcv_id"],
        "all_rcv_ids": cv_rec["all_rcv_ids"],
        "chr": cv_rec["chr"],
        "position": cv_rec["position"],
        "end_position": int(info_map.get("END", cv_rec["end_position"])),
        "variant_id": vid,
        "gene": cv_rec["gene"],
        "gene_name": cv_rec.get("gene_name", ""),
        "ref_allele": vcf_ref,
        "alt_allele": vcf_alt if vcf_alt != "." else cv_rec["alt_allele"],
        "risk_version": cv_rec["alt_allele"],
        "risk_frequency": freq_val,
        "risk_frequency_text": freq_text,
        "variant_type": vtype,
        "your_data": fmt["your_data"],
        "your_status": your_status,
        "genotype": fmt["gt"],
        "alleles_called": fmt["alleles_called"],
        "condition": cv_rec["condition"],
        "classification": cv_rec["classification"],
        "evidence": evidence_str,
        "evidence_stars": stars,
        "review_status": cv_rec["review_status"],
        "hgvs": cv_rec["hgvs"],
        "molecular_consequence": cv_rec.get("clinvar_type", ""),
        "dp": fmt["dp"],
        "gq": fmt["gq"],
        "ad": fmt["ad"],
        "qual": fmt["qual"],
        "filter_status": vcf_filt,
        "is_carried": 1 if is_carried else 0,
        "clinvar_match": "legacy_" + verdict,
        "source_file": SNP_VCF_NAME,
        "raw_vcf_line": "\t".join(vcf_parts),
    }


UPSERT_BASE_COLUMNS: List[str] = [
    "unique_key",
    "rcv_id",
    "all_rcv_ids",
    "chr",
    "position",
    "end_position",
    "variant_id",
    "gene",
    "gene_name",
    "ref_allele",
    "alt_allele",
    "risk_version",
    "risk_frequency",
    "risk_frequency_text",
    "variant_type",
    "your_data",
    "your_status",
    "genotype",
    "alleles_called",
    "condition",
    "classification",
    "evidence",
    "evidence_stars",
    "review_status",
    "hgvs",
    "molecular_consequence",
    "dp",
    "gq",
    "ad",
    "qual",
    "filter_status",
    "is_carried",
    "source_file",
    "raw_vcf_line",
]
# Fields that "Sync from NCBI" / NCBI import write. Once a row has been synced
# (ncbi_synced_at set), a later non-NCBI upsert (MyVariant resolver, catalog rebuild)
# must not silently revert them.
NCBI_SYNCED_FIELDS = {
    "all_rcv_ids",
    "condition",
    "classification",
    "evidence",
    "evidence_stars",
    "review_status",
    "risk_frequency",
    "risk_frequency_text",
}


def _upsert_sql() -> str:
    ncbi_cols = [name for name, _t in NCBI_EXTRA_COLUMNS]
    cols = UPSERT_BASE_COLUMNS + ncbi_cols
    sets: List[str] = []
    for c in cols:
        if c == "unique_key":
            continue
        if c in ncbi_cols:
            # NCBI provenance: only an NCBI-derived record (non-NULL value) may change it
            sets.append(f"{c} = COALESCE(excluded.{c}, variants.{c})")
        elif c in NCBI_SYNCED_FIELDS:
            sets.append(
                f"{c} = CASE WHEN variants.ncbi_synced_at IS NOT NULL AND excluded.ncbi_synced_at IS NULL "
                f"THEN variants.{c} ELSE excluded.{c} END"
            )
        else:
            sets.append(f"{c} = excluded.{c}")
    placeholders = ", ".join(["?"] * len(cols))
    # ON CONFLICT ... DO UPDATE keeps the row id stable (INSERT OR REPLACE deleted the row
    # and re-inserted it under a new AUTOINCREMENT id, breaking open drawers / links).
    return (
        f"INSERT INTO variants ({', '.join(cols)}) VALUES ({placeholders}) "
        f"ON CONFLICT(unique_key) DO UPDATE SET {', '.join(sets)}"
    )


def upsert_records(records: List[Dict[str, Any]]) -> int:
    if not records:
        return 0
    ensure_ncbi_columns()
    cols = UPSERT_BASE_COLUMNS + [name for name, _t in NCBI_EXTRA_COLUMNS]
    conn = get_db_connection()
    try:
        rows = [tuple(r.get(c) for c in cols) for r in records]
        conn.executemany(_upsert_sql(), rows)
        conn.commit()
    finally:
        conn.close()
    return len(rows)


def query_myvariant_batch(q: str, size: int = 1000) -> List[Dict[str, Any]]:
    if workspace.offline_mode():
        return []
    url = (
        "https://myvariant.info/v1/query?q="
        + urllib.parse.quote(q)
        + f"&assembly=hg38&size={size}&fields=clinvar,dbsnp.rsid,gnomad_genome.af.af,exac.af"
    )
    try:
        data = _http_get_json(url, timeout=20)
        return data.get("hits", [])
    except Exception:
        return []


def annotate_genomic_interval_genes(chrom: str, start_pos: int, end_pos: int) -> Tuple[str, str]:
    if workspace.offline_mode():
        return (f"CHR{chrom}_LOCUS", "Genomic structural interval")
    span = max(0, end_pos - start_pos)
    query_end = end_pos if span <= 2_500_000 else start_pos + 2_500_000
    url = f"https://mygene.info/v3/query?q=chr{chrom}:{start_pos}-{query_end}&species=human&size=5&fields=symbol,name"
    try:
        data = _http_get_json(url, timeout=8)
        hits = data.get("hits", [])
        symbols = []
        names = []
        for h in hits:
            sym = h.get("symbol")
            if sym and not sym.startswith("LOC") and sym not in symbols:
                symbols.append(sym)
                if h.get("name"):
                    names.append(h["name"])
        if not symbols and hits:
            for h in hits:
                if h.get("symbol") and h["symbol"] not in symbols:
                    symbols.append(h["symbol"])
        if symbols:
            return (", ".join(symbols[:3]), "; ".join(names[:2]))
    except Exception:
        pass
    return (f"CHR{chrom}_LOCUS", "Genomic structural interval")


def index_cnv_and_sv_files() -> int:
    records: List[Dict[str, Any]] = []

    # 1. Index non-reference CNVs from cnv.vcf
    if vcf_engine.cnv_vcf_path.exists():
        with open(vcf_engine.cnv_vcf_path, "r", encoding="utf-8", errors="replace") as f:
            for idx, line in enumerate(f):
                if line.startswith("#"):
                    continue
                raw_line = line.rstrip("\r\n")
                parts = raw_line.split("\t")
                if len(parts) < 10:
                    continue
                chrom = normalize_chrom(parts[0])
                if "_" in chrom or chrom.startswith("Un"):
                    continue
                pos = int(parts[1])
                vid = parts[2]
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                filt = parts[6]
                if alt == ".":
                    continue  # skip multi-megabase hom-ref blocks in CNV catalog
                info_map = parse_info_field(parts[7])
                sample_map = parse_sample_field(parts[8], parts[9])
                end_pos = int(info_map.get("END", pos))
                cnv_len = int(info_map.get("CNVLEN", max(1, end_pos - pos)))

                fmt = format_vcf_genotype_and_status(ref, alt, qual, filt, info_map, sample_map, "cnv")
                vtype = determine_variant_type(ref, alt, info_map, "cnv")

                if filt == "PASS":
                    gene_sym, gene_desc = annotate_genomic_interval_genes(chrom, pos, end_pos)
                else:
                    gene_sym, gene_desc = (f"CHR{chrom}_CNV", "Copy number interval")

                svtype = info_map.get("SVTYPE", "CNV")
                if svtype == "LOH":
                    classification = "Risk Factor"
                    condition = f"Region of Loss of Heterozygosity ({cnv_len:,} bp) spanning {gene_sym}"
                    evidence = f"WGS Canvas LOH Call | {filt} (QS {qual})"
                    stars = 2 if filt == "PASS" else 0
                elif "LOSS" in vid or "<CN0>" in alt:
                    classification = "Likely Pathogenic" if (filt == "PASS" and cnv_len > 30000 and gene_sym and not gene_sym.startswith("CHR")) else "Uncertain Significance"
                    condition = f"Copy Number Deletion ({cnv_len:,} bp) affecting {gene_sym}" + (f" ({gene_desc})" if gene_desc else "")
                    evidence = f"WGS Canvas CNV Read | {filt} (QS {qual}, Bins {sample_map.get('BC', '.')})"
                    stars = 2 if filt == "PASS" else 0
                else:
                    classification = "Uncertain Significance"
                    condition = f"Copy Number Duplication/Gain ({cnv_len:,} bp) affecting {gene_sym}"
                    evidence = f"WGS Canvas CNV Read | {filt} (QS {qual}, Bins {sample_map.get('BC', '.')})"
                    stars = 2 if filt == "PASS" else 0

                rcv_syn = f"RCV-CNV-{chrom}-{pos}"
                records.append(
                    {
                        "unique_key": f"CNV:{chrom}:{pos}:{vid}",
                        "rcv_id": rcv_syn,
                        "all_rcv_ids": rcv_syn,
                        "chr": chrom,
                        "position": pos,
                        "end_position": end_pos,
                        "variant_id": vid,
                        "gene": gene_sym,
                        "gene_name": gene_desc,
                        "ref_allele": ref,
                        "alt_allele": alt,
                        "risk_version": alt,
                        "risk_frequency": 0.005 if filt == "PASS" else None,
                        "risk_frequency_text": "0.500% (Structural CNV)" if filt == "PASS" else "Rare / Structural",
                        "variant_type": vtype,
                        "your_data": fmt["your_data"],
                        "your_status": fmt["your_status"],
                        "genotype": fmt["gt"],
                        "alleles_called": fmt["alleles_called"],
                        "condition": condition,
                        "classification": classification,
                        "evidence": evidence,
                        "evidence_stars": stars,
                        "review_status": "WGS Canvas Copy Number Caller",
                        "hgvs": f"NC_0000{chrom.zfill(2) if chrom.isdigit() else chrom}:g.{pos}_{end_pos}{svtype.lower()}",
                        "molecular_consequence": vtype,
                        "dp": fmt["dp"],
                        "gq": fmt["gq"],
                        "ad": fmt["ad"],
                        "qual": fmt["qual"],
                        "filter_status": filt,
                        "is_carried": 1 if fmt["is_carried"] else 0,
                        "source_file": CNV_VCF_NAME,
                        "raw_vcf_line": raw_line,
                    }
                )

    # 2. Index high-quality PASS Structural Variants from sv.vcf across chromosomes
    if vcf_engine.sv_vcf_path.exists():
        sv_candidates: List[List[str]] = []
        with open(vcf_engine.sv_vcf_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 10:
                    continue
                chrom = normalize_chrom(parts[0])
                if "_" in chrom or chrom.startswith("Un"):
                    continue
                if parts[6] != "PASS":
                    continue
                info_map = parse_info_field(parts[7])
                svtype = info_map.get("SVTYPE", "")
                if svtype in ("DEL", "INS", "DUP"):
                    try:
                        qual_f = float(parts[5]) if parts[5] != "." else 0.0
                    except ValueError:
                        qual_f = 0.0
                    if qual_f >= 400:
                        sv_candidates.append(parts)

        # Pick top ~60 distributed across chromosomes so all SV types are represented with gene annotations
        by_chrom: Dict[str, List[List[str]]] = {}
        for p in sv_candidates:
            by_chrom.setdefault(p[0], []).append(p)

        selected_svs: List[List[str]] = []
        for ch, items in by_chrom.items():
            selected_svs.extend(items[:3])

        for parts in selected_svs[:65]:
            chrom = normalize_chrom(parts[0])
            pos = int(parts[1])
            vid = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filt = parts[6]
            info_map = parse_info_field(parts[7])
            sample_map = parse_sample_field(parts[8], parts[9])
            end_pos = int(info_map.get("END", pos))
            svlen = info_map.get("SVLEN", str(end_pos - pos))

            fmt = format_vcf_genotype_and_status(ref, alt, qual, filt, info_map, sample_map, "sv")
            vtype = determine_variant_type(ref, alt, info_map, "sv")
            gene_sym, gene_desc = annotate_genomic_interval_genes(chrom, pos, max(pos + 1, end_pos))

            rcv_syn = f"RCV-SV-{chrom}-{pos}"
            records.append(
                {
                    "unique_key": f"SV:{chrom}:{pos}:{vid}",
                    "rcv_id": rcv_syn,
                    "all_rcv_ids": rcv_syn,
                    "chr": chrom,
                    "position": pos,
                    "end_position": end_pos,
                    "variant_id": vid,
                    "gene": gene_sym,
                    "gene_name": gene_desc,
                    "ref_allele": ref[:30],
                    "alt_allele": alt[:30],
                    "risk_version": fmt["risk_version"],
                    "risk_frequency": 0.002,
                    "risk_frequency_text": "0.200% (Structural Variant)",
                    "variant_type": vtype,
                    "your_data": fmt["your_data"],
                    "your_status": fmt["your_status"],
                    "genotype": fmt["gt"],
                    "alleles_called": fmt["alleles_called"],
                    "condition": f"Structural {vtype} ({svlen} bp) at {gene_sym}" + (f" ({gene_desc})" if gene_desc else ""),
                    "classification": "Uncertain Significance" if gene_sym.startswith("CHR") else "Risk Factor",
                    "evidence": f"WGS Manta SV Caller | PASS (QUAL {qual}, GQ {sample_map.get('GQ', '.')})",
                    "evidence_stars": 2,
                    "review_status": "WGS Manta Paired/Split-Read Support",
                    "hgvs": f"chr{chrom}:g.{pos}_{end_pos}({svlen}bp)",
                    "molecular_consequence": vtype,
                    "dp": fmt["dp"],
                    "gq": fmt["gq"],
                    "ad": fmt["ad"],
                    "qual": fmt["qual"],
                    "filter_status": filt,
                    "is_carried": 1 if fmt["is_carried"] else 0,
                    "source_file": SV_VCF_NAME,
                    "raw_vcf_line": "\t".join(parts),
                }
            )

    return upsert_records(records)


def build_initial_database(force: bool = False) -> Dict[str, Any]:
    init_db()
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM variants")
    existing_count = cur.fetchone()[0]
    conn.close()

    if existing_count > 500 and not force:
        return {"status": "cached", "total_records": existing_count}

    t0 = time.time()
    vcf_engine.ensure_loaded()

    queries = [
        'clinvar.rcv.clinical_significance:("Risk factor" OR "Drug response" OR Association OR Protective) AND _exists_:dbsnp.rsid',
        'clinvar.rcv.clinical_significance:(Pathogenic OR "Likely pathogenic") AND gnomad_genome.af.af:>0.0003 AND _exists_:dbsnp.rsid',
        'clinvar.gene.symbol:(BRCA1 OR BRCA2 OR APOE OR MTHFR OR CFTR OR HFE OR ALDH2 OR CYP2C19 OR CYP2D6 OR CYP2C9 OR VKORC1 OR SLCO1B1 OR F5 OR F2 OR SERPINA1 OR GJB2 OR PAH OR LDLR OR PCSK9 OR TP53 OR COMT OR BDNF OR VDR OR LCT OR MCM6 OR ACTN3 OR HERC2 OR OCA2 OR MC1R OR UGT1A1 OR DPYD OR TPMT OR NUDT15 OR CYP3A5 OR NAT2 OR BCHE) AND gnomad_genome.af.af:>0.001 AND _exists_:dbsnp.rsid',
    ]

    all_hits: List[Dict[str, Any]] = []
    for q in queries:
        all_hits.extend(query_myvariant_batch(q, size=1000))

    clinvar_candidates: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for h in all_hits:
        extracted = extract_clinvar_records_from_hit(h)
        for rec in extracted:
            key = (rec["chr"], rec["position"])
            if key not in clinvar_candidates:
                clinvar_candidates[key] = rec

    records_to_insert: List[Dict[str, Any]] = []
    if vcf_engine.snp_vcf_path.exists():
        with open(vcf_engine.snp_vcf_path, "rb") as f_handle:
            for (ch, pos), cv_rec in clinvar_candidates.items():
                vcf_parts, _ctx = vcf_engine.lookup_snp_indel_locus(ch, pos, f_handle=f_handle, include_context=0)
                if vcf_parts:
                    row = build_record_from_clinvar_and_vcf(cv_rec, vcf_parts)
                    if row:
                        records_to_insert.append(row)

    snp_inserted = upsert_records(records_to_insert)
    cnv_sv_inserted = index_cnv_and_sv_files()

    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM variants")
    total = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM variants WHERE is_carried = 1")
    carried = cur.fetchone()[0]
    conn.close()

    return {
        "status": "built",
        "snp_indel_indexed": snp_inserted,
        "cnv_sv_indexed": cnv_sv_inserted,
        "total_records": total,
        "carried_variants": carried,
        "elapsed_seconds": round(time.time() - t0, 2),
    }


def resolve_live_query(
    rcv_id: Optional[str] = None,
    variant_id: Optional[str] = None,
    gene: Optional[str] = None,
    chr_val: Optional[str] = None,
    pos_val: Optional[str] = None,
    condition: Optional[str] = None,
) -> int:
    """
    On-demand live resolver:
    1. If chr + position is provided, directly reads the VCF at that coordinate/range and enriches via MyVariant.info.
    2. If rcv_id, variant_id (rsID), gene, or condition is provided and not yet well-represented in SQLite,
       queries MyVariant.info (GRCh38), reads the corresponding coordinates from the analysis' VCFs, and caches in SQLite.
    """
    new_records: List[Dict[str, Any]] = []

    # Case A: Direct MyVariant.info query for RCV ID, Variant ID, Gene, or Condition
    mv_clauses: List[str] = []
    if rcv_id and rcv_id.upper().startswith("RCV") and not rcv_id.upper().startswith("RCV-"):
        mv_clauses.append(f'clinvar.rcv.accession:"{rcv_id.strip().upper()}"')
    if variant_id and variant_id.lower().startswith("rs"):
        mv_clauses.append(f'dbsnp.rsid:"{variant_id.strip().lower()}"')
    if gene and len(gene.strip()) >= 2:
        mv_clauses.append(f'clinvar.gene.symbol:"{gene.strip().upper()}"')
    if condition and len(condition.strip()) >= 3 and not mv_clauses:
        mv_clauses.append(f'clinvar.rcv.conditions.name:"{condition.strip()}"')

    if mv_clauses:
        q = " AND ".join(mv_clauses)
        hits = query_myvariant_batch(q, size=100)
        if hits and vcf_engine.snp_vcf_path.exists():
            with open(vcf_engine.snp_vcf_path, "rb") as f_handle:
                for h in hits:
                    for cv_rec in extract_clinvar_records_from_hit(h):
                        vcf_parts, _ctx = vcf_engine.lookup_snp_indel_locus(
                            cv_rec["chr"], cv_rec["position"], f_handle=f_handle
                        )
                        if vcf_parts:
                            row = build_record_from_clinvar_and_vcf(cv_rec, vcf_parts)
                            if row:
                                new_records.append(row)

    # Case B: Direct VCF coordinate read if chr and pos_val are specified
    if chr_val and pos_val:
        chrom = normalize_chrom(chr_val)
        pos_clean = pos_val.replace(",", "").strip()
        start_p: Optional[int] = None
        end_p: Optional[int] = None
        if "-" in pos_clean:
            p_parts = pos_clean.split("-", 1)
            if p_parts[0].strip().isdigit() and p_parts[1].strip().isdigit():
                start_p = int(p_parts[0].strip())
                end_p = int(p_parts[1].strip())
        elif pos_clean.isdigit():
            start_p = int(pos_clean)
            end_p = start_p

        if start_p is not None and end_p is not None:
            # Query MyVariant for any ClinVar records in this region first
            mv_region_q = f"chr{chrom}:{start_p}-{min(end_p, start_p + 100000)} AND _exists_:clinvar.rcv"
            hits = query_myvariant_batch(mv_region_q, size=50)
            if hits and vcf_engine.snp_vcf_path.exists():
                with open(vcf_engine.snp_vcf_path, "rb") as f_handle:
                    for h in hits:
                        for cv_rec in extract_clinvar_records_from_hit(h):
                            vcf_parts, _ctx = vcf_engine.lookup_snp_indel_locus(
                                cv_rec["chr"], cv_rec["position"], f_handle=f_handle
                            )
                            if vcf_parts:
                                row = build_record_from_clinvar_and_vcf(cv_rec, vcf_parts)
                                if row:
                                    new_records.append(row)

            # Also read raw VCF rows directly at [start_p, end_p] so even non-ClinVar loci return real VCF reads
            raw_reads = vcf_engine.read_vcf_region(chrom, start_p, end_p, variants_only=False, limit=30)
            gene_sym, gene_desc = annotate_genomic_interval_genes(chrom, start_p, end_p)
            for r in raw_reads:
                vid = r["variant_id"]
                unique_key = f"RAW:{r['source_file']}:{r['chr']}:{r['position']}:{vid}"
                new_records.append(
                    {
                        "unique_key": unique_key,
                        "rcv_id": "VCF-DIRECT-READ",
                        "all_rcv_ids": "VCF-DIRECT-READ",
                        "chr": r["chr"],
                        "position": r["position"],
                        "end_position": r["end_position"],
                        "variant_id": vid,
                        "gene": gene_sym,
                        "gene_name": gene_desc,
                        "ref_allele": r["ref"],
                        "alt_allele": r["alt"],
                        "risk_version": r["risk_version"],
                        "risk_frequency": None,
                        "risk_frequency_text": "Direct WGS Read",
                        "variant_type": r["variant_type"],
                        "your_data": r["your_data"],
                        "your_status": r["your_status"],
                        "genotype": r["your_data"].split("(")[1].split(")")[0] if "(" in r["your_data"] else "./.",
                        "alleles_called": r["your_data"].split("(")[0].strip(),
                        "condition": f"Direct VCF Read at {gene_sym} (chr{r['chr']}:{r['position']:,})",
                        "classification": "Variant Call" if r["is_carried"] else "Reference Block",
                        "evidence": f"Direct WGS VCF Read | {r['filter']} (DP {r['dp'] or '.'}, GQ {r['gq'] or '.'})",
                        "evidence_stars": 2 if r["filter"] == "PASS" else 0,
                        "review_status": "Direct WGS VCF Read",
                        "hgvs": f"chr{r['chr']}:g.{r['position']}{r['ref']}>{r['alt']}",
                        "molecular_consequence": r["variant_type"],
                        "dp": r["dp"],
                        "gq": r["gq"],
                        "ad": r["ad"],
                        "qual": r["qual"],
                        "filter_status": r["filter"],
                        "is_carried": 1 if r["is_carried"] else 0,
                        "source_file": r["source_file"],
                        "raw_vcf_line": r["raw_vcf_line"],
                    }
                )

    return upsert_records(new_records)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="(Re)build the variant catalog of one analysis workspace.")
    ap.add_argument("--workspace", type=Path, required=True, help="analysis job folder (contains inputs/ and work/)")
    args = ap.parse_args()
    with workspace.activate(workspace.Workspace(args.workspace)):
        res = build_initial_database(force=True)
    print("Database build complete:", json.dumps(res, indent=2))
