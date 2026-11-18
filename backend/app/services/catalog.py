"""
Catalog queries for the active analysis (summary, filter options, multi-field search,
live-resolution helpers). Moved unchanged from the original API module; the hard-coded
sample id / source file name are now taken from the active workspace.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from genomescope.indexer import get_db_connection
from genomescope.ncbi_client import is_gene_symbol
from genomescope.vcf_engine import SNP_VCF_NAME, engine_for, normalize_chrom
from genomescope.workspace import Workspace


def parse_position_filter(pos_str: str) -> Tuple[Optional[int], Optional[int]]:
    cleaned = pos_str.replace(",", "").strip()
    if ":" in cleaned:
        cleaned = cleaned.split(":", 1)[1].strip()
    if "-" in cleaned:
        parts = cleaned.split("-", 1)
        if parts[0].strip().isdigit() and parts[1].strip().isdigit():
            return (int(parts[0].strip()), int(parts[1].strip()))
    if cleaned.isdigit():
        val = int(cleaned)
        return (val, val)
    return (None, None)


def apply_risk_frequency_clause(
    rf_str: str, where_clauses: List[str], params: List[Any]
) -> None:
    s = rf_str.strip().lower().replace("%", "")
    if not s or s == "all":
        return
    if s in ("rare", "<0.01", "< 0.01", "<1"):
        where_clauses.append("(risk_frequency IS NOT NULL AND risk_frequency < 0.01)")
        return
    if s in ("uncommon", "<0.05", "< 0.05", "0.01-0.05"):
        where_clauses.append("(risk_frequency IS NOT NULL AND risk_frequency < 0.05)")
        return
    if s in ("common", ">0.05", "> 0.05", ">=0.05"):
        where_clauses.append("(risk_frequency IS NOT NULL AND risk_frequency >= 0.05)")
        return

    # Range like 0.01-0.10
    if "-" in s:
        p1, p2 = s.split("-", 1)
        try:
            v1, v2 = float(p1.strip()), float(p2.strip())
            if v1 > 1.0 or v2 > 1.0:
                v1, v2 = v1 / 100.0, v2 / 100.0
            where_clauses.append("(risk_frequency IS NOT NULL AND risk_frequency >= ? AND risk_frequency <= ?)")
            params.extend([v1, v2])
            return
        except ValueError:
            pass

    m = re.match(r"^(<=|>=|<|>|=)?\s*([0-9]*\.?[0-9]+)$", s)
    if m:
        op = m.group(1) or "<="
        val = float(m.group(2))
        if val > 1.0:
            val = val / 100.0
        if op == "=":
            where_clauses.append(" LOWER(risk_frequency_text) LIKE ? ")
            params.append(f"%{s}%")
        else:
            where_clauses.append(f"(risk_frequency IS NOT NULL AND risk_frequency {op} ?)")
            params.append(val)
    else:
        where_clauses.append("LOWER(risk_frequency_text) LIKE ?")
        params.append(f"%{s}%")


def build_summary(ws: Workspace) -> Dict[str, Any]:
    conn = get_db_connection()
    cur = conn.cursor()

    cur.execute("SELECT COUNT(*) FROM variants")
    total_records = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM variants WHERE is_carried = 1")
    carried_total = cur.fetchone()[0]

    cur.execute(
        "SELECT COUNT(*) FROM variants WHERE is_carried = 1 AND (classification LIKE 'Pathogenic%' OR classification = 'Likely Pathogenic')"
    )
    pathogenic_carried = cur.fetchone()[0]

    cur.execute(
        "SELECT COUNT(*) FROM variants WHERE is_carried = 1 AND classification IN ('Risk Factor', 'Association')"
    )
    risk_factor_carried = cur.fetchone()[0]

    cur.execute(
        "SELECT COUNT(*) FROM variants WHERE is_carried = 1 AND classification LIKE '%Drug Response%'"
    )
    drug_response_carried = cur.fetchone()[0]

    cur.execute(
        "SELECT COUNT(*) FROM variants WHERE is_carried = 1 AND source_file != ?", (SNP_VCF_NAME,)
    )
    structural_cnv_carried = cur.fetchone()[0]

    cur.execute(
        "SELECT classification, COUNT(*) as total, SUM(is_carried) as carried FROM variants GROUP BY classification ORDER BY carried DESC, total DESC"
    )
    by_classification = [dict(r) for r in cur.fetchall()]

    cur.execute(
        "SELECT your_status, COUNT(*) as count FROM variants GROUP BY your_status ORDER BY count DESC"
    )
    by_status = [dict(r) for r in cur.fetchall()]

    cur.execute(
        "SELECT variant_type, COUNT(*) as total, SUM(is_carried) as carried FROM variants GROUP BY variant_type ORDER BY total DESC"
    )
    by_variant_type = [dict(r) for r in cur.fetchall()]

    cur.execute(
        "SELECT gene, COUNT(*) as carried_count FROM variants WHERE is_carried = 1 AND gene != 'Intergenic' AND gene NOT LIKE 'CHR%' GROUP BY gene ORDER BY carried_count DESC LIMIT 18"
    )
    top_carried_genes = [dict(r) for r in cur.fetchall()]

    conn.close()

    return {
        "sample_id": ws.sample_label,
        "reference_genome": ws.reference,
        "total_indexed_loci": total_records,
        "carried_variants": carried_total,
        "pathogenic_or_likely_carried": pathogenic_carried,
        "risk_factor_carried": risk_factor_carried,
        "drug_response_carried": drug_response_carried,
        "structural_cnv_carried": structural_cnv_carried,
        "by_classification": by_classification,
        "by_status": by_status,
        "by_variant_type": by_variant_type,
        "top_carried_genes": top_carried_genes,
        "vcf_files": engine_for(ws).get_vcf_files_metadata(),
    }


def filter_options() -> Dict[str, Any]:
    conn = get_db_connection()
    cur = conn.cursor()

    cur.execute("SELECT DISTINCT variant_type FROM variants WHERE variant_type IS NOT NULL ORDER BY variant_type")
    variant_types = [r[0] for r in cur.fetchall()]

    cur.execute("SELECT DISTINCT your_status FROM variants WHERE your_status IS NOT NULL ORDER BY your_status")
    statuses = [r[0] for r in cur.fetchall()]

    cur.execute("SELECT DISTINCT classification FROM variants WHERE classification IS NOT NULL ORDER BY classification")
    classifications = [r[0] for r in cur.fetchall()]

    cur.execute("SELECT DISTINCT source_file FROM variants WHERE source_file IS NOT NULL ORDER BY source_file")
    source_files = [r[0] for r in cur.fetchall()]

    cur.execute(
        "SELECT gene, COUNT(*) as c FROM variants WHERE gene != 'Intergenic' AND gene NOT LIKE 'CHR%' GROUP BY gene ORDER BY c DESC LIMIT 60"
    )
    popular_genes = [r[0] for r in cur.fetchall()]

    conn.close()

    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    evidence_options = [
        "Practice Guideline (4 Stars)",
        "Expert Panel Reviewed (3 Stars)",
        "Strong - Multiple Submitters (2 Stars)",
        "Moderate - Single Submitter (1 Star)",
        "Conflicting Interpretations (1 Star)",
        "VCF PASS",
    ]
    risk_frequency_presets = [
        {"label": "All Frequencies", "value": ""},
        {"label": "Very Rare (< 0.1% / < 0.001)", "value": "<0.001"},
        {"label": "Rare (< 1% / < 0.01)", "value": "<0.01"},
        {"label": "Low Frequency (< 5% / < 0.05)", "value": "<0.05"},
        {"label": "Moderate (1% - 10% / 0.01-0.10)", "value": "0.01-0.10"},
        {"label": "Common (>= 5% / >= 0.05)", "value": ">=0.05"},
    ]

    return {
        "chromosomes": chromosomes,
        "variant_types": variant_types,
        "statuses": statuses,
        "classifications": classifications,
        "evidence_options": evidence_options,
        "source_files": source_files,
        "popular_genes": popular_genes,
        "risk_frequency_presets": risk_frequency_presets,
    }


def execute_variant_search(
    q: Optional[str],
    rcv_id: Optional[str],
    chr_val: Optional[str],
    position: Optional[str],
    variant_id: Optional[str],
    gene: Optional[str],
    risk_frequency: Optional[str],
    variant_type: Optional[str],
    your_status: Optional[str],
    condition: Optional[str],
    classification: Optional[str],
    evidence: Optional[str],
    source_file: Optional[str],
    carried_only: bool,
    sort_by: str,
    sort_dir: str,
    page: int,
    page_size: int,
) -> Tuple[int, List[Dict[str, Any]]]:
    where_clauses: List[str] = []
    params: List[Any] = []

    if carried_only:
        where_clauses.append("is_carried = 1")

    if q and q.strip():
        q_clean = q.strip()
        # Check if q is chr:pos format like "12:111803962" or "chr12:111803962"
        m_locus = re.match(r"^(?:chr)?([0-9XYMThm]+):([0-9,]+(?:\-[0-9,]+)?)$", q_clean, re.IGNORECASE)
        if m_locus:
            c_norm = normalize_chrom(m_locus.group(1))
            p_start, p_end = parse_position_filter(m_locus.group(2))
            where_clauses.append("chr = ?")
            params.append(c_norm)
            if p_start is not None and p_end is not None:
                where_clauses.append("(position <= ? AND end_position >= ?)")
                params.extend([p_end, p_start])
        else:
            like_q = f"%{q_clean.lower()}%"
            where_clauses.append(
                """(
                    LOWER(rcv_id) LIKE ? OR
                    LOWER(all_rcv_ids) LIKE ? OR
                    LOWER(variant_id) LIKE ? OR
                    LOWER(gene) LIKE ? OR
                    LOWER(condition) LIKE ? OR
                    LOWER(classification) LIKE ? OR
                    LOWER(your_status) LIKE ? OR
                    LOWER(variant_type) LIKE ? OR
                    LOWER(evidence) LIKE ? OR
                    CAST(position AS TEXT) LIKE ?
                )"""
            )
            params.extend([like_q] * 10)

    if rcv_id and rcv_id.strip():
        r_like = f"%{rcv_id.strip().lower()}%"
        where_clauses.append("(LOWER(rcv_id) LIKE ? OR LOWER(all_rcv_ids) LIKE ?)")
        params.extend([r_like, r_like])

    if chr_val and chr_val.strip() and chr_val.strip().lower() != "all":
        c_norm = normalize_chrom(chr_val)
        where_clauses.append("UPPER(chr) = ?")
        params.append(c_norm.upper())

    if position and position.strip():
        p_start, p_end = parse_position_filter(position)
        if p_start is not None and p_end is not None:
            if p_start == p_end:
                where_clauses.append("((position = ?) OR (position <= ? AND end_position >= ?))")
                params.extend([p_start, p_start, p_start])
            else:
                where_clauses.append("(position <= ? AND end_position >= ?)")
                params.extend([p_end, p_start])
        else:
            where_clauses.append("CAST(position AS TEXT) LIKE ?")
            params.append(f"%{position.strip()}%")

    if variant_id and variant_id.strip():
        where_clauses.append("LOWER(variant_id) LIKE ?")
        params.append(f"%{variant_id.strip().lower()}%")

    if gene and gene.strip():
        g_like = f"%{gene.strip().lower()}%"
        where_clauses.append("(LOWER(gene) LIKE ? OR LOWER(gene_name) LIKE ?)")
        params.extend([g_like, g_like])

    if risk_frequency and risk_frequency.strip():
        apply_risk_frequency_clause(risk_frequency, where_clauses, params)

    if variant_type and variant_type.strip() and variant_type.strip().lower() != "all":
        where_clauses.append("LOWER(variant_type) LIKE ?")
        params.append(f"%{variant_type.strip().lower()}%")

    if your_status and your_status.strip() and your_status.strip().lower() != "all":
        where_clauses.append("LOWER(your_status) LIKE ?")
        params.append(f"%{your_status.strip().lower()}%")

    if condition and condition.strip():
        where_clauses.append("LOWER(condition) LIKE ?")
        params.append(f"%{condition.strip().lower()}%")

    if classification and classification.strip() and classification.strip().lower() != "all":
        where_clauses.append("LOWER(classification) LIKE ?")
        params.append(f"%{classification.strip().lower()}%")

    if evidence and evidence.strip() and evidence.strip().lower() != "all":
        where_clauses.append("LOWER(evidence) LIKE ?")
        params.append(f"%{evidence.strip().lower()}%")

    if source_file and source_file.strip() and source_file.strip().lower() != "all":
        where_clauses.append("LOWER(source_file) LIKE ?")
        params.append(f"%{source_file.strip().lower()}%")

    where_sql = f"WHERE {' AND '.join(where_clauses)}" if where_clauses else ""

    allowed_sorts = {
        "priority": "is_carried DESC, evidence_stars DESC, CASE WHEN classification LIKE 'Pathogenic%' THEN 1 WHEN classification = 'Likely Pathogenic' THEN 2 WHEN classification LIKE '%Drug Response%' THEN 3 WHEN classification = 'Risk Factor' THEN 4 ELSE 5 END ASC, position ASC",
        "position": "CAST(chr AS INTEGER) ASC, position ASC",
        "gene": "gene ASC, position ASC",
        "risk_frequency": "CASE WHEN risk_frequency IS NULL THEN 1 ELSE 0 END ASC, risk_frequency DESC",
        "evidence": "evidence_stars DESC, is_carried DESC",
        "classification": "classification ASC, is_carried DESC",
        "your_status": "is_carried DESC, your_status ASC",
        "variant_id": "variant_id ASC",
        "rcv_id": "rcv_id ASC",
    }

    if sort_by in allowed_sorts and sort_by != "priority":
        direction = "DESC" if sort_dir.lower() == "desc" else "ASC"
        if sort_by in ("position", "gene", "classification", "variant_id", "rcv_id"):
            order_sql = f"ORDER BY {sort_by} {direction}, id ASC"
        else:
            order_sql = f"ORDER BY {allowed_sorts[sort_by]}"
    else:
        order_sql = f"ORDER BY {allowed_sorts['priority']}"

    conn = get_db_connection()
    cur = conn.cursor()

    cur.execute(f"SELECT COUNT(*) FROM variants {where_sql}", params)
    total = cur.fetchone()[0]

    offset = (max(1, page) - 1) * page_size
    cur.execute(
        f"SELECT * FROM variants {where_sql} {order_sql} LIMIT ? OFFSET ?",
        params + [page_size, offset],
    )
    rows = [dict(r) for r in cur.fetchall()]
    conn.close()
    return (total, rows)


_LIVE_RCV_RE = re.compile(r"^RCV\d{1,12}(?:\.\d+)?$", re.IGNORECASE)
_LIVE_RS_RE = re.compile(r"^rs\d{1,12}$", re.IGNORECASE)
_LIVE_VCV_RE = re.compile(r"^VCV\d{1,12}(?:\.\d+)?$", re.IGNORECASE)
_LIVE_LOCUS_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M):([0-9,]{1,15}(?:-[0-9,]{1,15})?)$", re.IGNORECASE)


def live_identifiers(
    q: Optional[str],
    rcv_id: Optional[str],
    variant_id: Optional[str],
    gene: Optional[str],
    chr_val: Optional[str],
    position: Optional[str],
) -> Dict[str, Optional[str]]:
    """Identifiers a zero-result catalog search may resolve live (never free text)."""
    live: Dict[str, Optional[str]] = {
        "rcv_id": (rcv_id or "").strip() or None,
        "variant_id": (variant_id or "").strip() or None,
        "gene": (gene or "").strip() or None,
        "vcv_id": None,
        "chr": (chr_val or "").strip() or None if (position or "").strip() else None,
        "position": (position or "").strip() or None if (chr_val or "").strip() else None,
    }
    if q and q.strip() and not any([rcv_id, variant_id, gene, position]):
        q_s = q.strip()
        m = _LIVE_LOCUS_RE.match(q_s)
        if _LIVE_RCV_RE.match(q_s):
            live["rcv_id"] = q_s
        elif _LIVE_RS_RE.match(q_s):
            live["variant_id"] = q_s
        elif _LIVE_VCV_RE.match(q_s):
            live["vcv_id"] = q_s
        elif m:
            live["chr"], live["position"] = m.group(1), m.group(2)
        elif is_gene_symbol(q_s):
            live["gene"] = q_s
    return live


def identifiers_in_catalog(live: Dict[str, Optional[str]], condition: Optional[str]) -> bool:
    """True when the identifier (ignoring every other filter) already exists in the catalog."""
    if live.get("vcv_id"):
        m = re.search(r"(\d+)", live["vcv_id"] or "")
        conn = get_db_connection()
        try:
            row = conn.execute("SELECT 1 FROM variants WHERE ncbi_variation_id = ? LIMIT 1", (str(int(m.group(1))) if m else "",)).fetchone()
        except Exception:
            row = None
        finally:
            conn.close()
        return row is not None
    total, _rows = execute_variant_search(
        q=None,
        rcv_id=live.get("rcv_id"),
        chr_val=live.get("chr"),
        position=live.get("position"),
        variant_id=live.get("variant_id"),
        gene=live.get("gene"),
        risk_frequency=None,
        variant_type=None,
        your_status=None,
        condition=condition,
        classification=None,
        evidence=None,
        source_file=None,
        carried_only=False,
        sort_by="priority",
        sort_dir="desc",
        page=1,
        page_size=1,
    )
    return total > 0
