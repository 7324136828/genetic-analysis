"""
Variant explorer API (catalog search, direct VCF reads, NCBI live endpoints).

Every route works on the analysis selected by the ``X-GenomeScope-Job`` header (see
app.api.deps). Behaviour is unchanged from the original application.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Body, Depends, HTTPException, Query
from pydantic import BaseModel

from app.api.deps import optional_analysis, require_app_header, use_analysis
from app.services import catalog
from genomescope import clinvar_scan, ncbi_service
from genomescope.indexer import get_db_connection, resolve_live_query
from genomescope.vcf_engine import normalize_chrom, vcf_engine
from genomescope.workspace import Workspace

# Routes that need a selected analysis
router = APIRouter(prefix="/api", tags=["explorer"], dependencies=[Depends(use_analysis)])
# Routes that also answer without one (status pills)
status_router = APIRouter(prefix="/api", tags=["status"])


@status_router.get("/clinvar/status")
def clinvar_status(ws: Optional[Workspace] = Depends(optional_analysis)) -> Dict[str, Any]:
    """Local ClinVar release + (when an analysis is selected) its verification state."""
    return clinvar_scan.status()


@status_router.get("/ncbi/status")
def ncbi_status(ws: Optional[Workspace] = Depends(optional_analysis)) -> Dict[str, Any]:
    return ncbi_service.status_body()


@router.get("/summary")
def get_summary(ws: Workspace = Depends(use_analysis)) -> Dict[str, Any]:
    return catalog.build_summary(ws)


@router.get("/filters")
def get_filter_options() -> Dict[str, Any]:
    return catalog.filter_options()


@router.get("/variants")
def search_variants(
    q: Optional[str] = Query(None, description="Global search query"),
    rcv_id: Optional[str] = Query(None, description="ClinVar RCV accession ID"),
    chr: Optional[str] = Query(None, description="Chromosome (1-22, X, Y, MT)"),
    position: Optional[str] = Query(None, description="Genomic position or range (e.g. 111803962 or 111800000-111900000)"),
    variant_id: Optional[str] = Query(None, description="Variant ID (rsID, Manta SV ID, Canvas CNV ID)"),
    gene: Optional[str] = Query(None, description="Gene symbol (e.g. ALDH2, BRCA1, CYP2C19, MTHFR)"),
    risk_frequency: Optional[str] = Query(None, description="Risk allele frequency filter (e.g. <0.01, <0.05, 0.01-0.10)"),
    variant_type: Optional[str] = Query(None, description="Variant type (SNV, Deletion, Insertion, CNV Loss, CNV Gain, LOH, SV)"),
    your_status: Optional[str] = Query(None, description="Personal genotype status from VCF"),
    condition: Optional[str] = Query(None, description="Associated clinical condition or phenotype"),
    classification: Optional[str] = Query(None, description="Clinical significance classification"),
    evidence: Optional[str] = Query(None, description="ClinVar review evidence or VCF quality filter"),
    source_file: Optional[str] = Query(None, description="Source VCF filename"),
    carried_only: bool = Query(False, description="Only return variants carried by the sample (0/1, 1/1, PASS SV/CNV)"),
    sort_by: str = Query("priority", description="Sort field"),
    sort_dir: str = Query("desc", description="Sort direction (asc or desc)"),
    page: int = Query(1, ge=1),
    page_size: int = Query(25, ge=1, le=200),
) -> Dict[str, Any]:
    search_args = dict(
        q=q, rcv_id=rcv_id, chr_val=chr, position=position, variant_id=variant_id, gene=gene,
        risk_frequency=risk_frequency, variant_type=variant_type, your_status=your_status,
        condition=condition, classification=classification, evidence=evidence, source_file=source_file,
        carried_only=carried_only, sort_by=sort_by, sort_dir=sort_dir, page=page, page_size=page_size,
    )
    total, items = catalog.execute_variant_search(**search_args)

    # Live resolution (MyVariant + VCF, then NCBI ClinVar import) only when the queried
    # identifier itself is absent from the catalog. Other filters (carried_only,
    # classification, ...) returning nothing must not trigger writes to the catalog.
    if total == 0:
        live = catalog.live_identifiers(q, rcv_id, variant_id, gene, chr, position)
        use_condition = bool(condition and condition.strip()) and not any(live.values())
        if (any(live.values()) or use_condition) and not catalog.identifiers_in_catalog(live, condition if use_condition else None):
            added = resolve_live_query(
                rcv_id=live["rcv_id"],
                variant_id=live["variant_id"],
                gene=live["gene"],
                chr_val=live["chr"],
                pos_val=live["position"],
                condition=condition if use_condition else None,
            )
            if added == 0 and any(live.values()):
                # Final fallback: live NCBI ClinVar import for identifier lookups only
                try:
                    added = ncbi_service.fallback_import_for_search(
                        rcv_id=live["rcv_id"],
                        variant_id=live["variant_id"],
                        gene=live["gene"],
                        vcv_id=live["vcv_id"],
                        chrom=live["chr"],
                        position=live["position"],
                    )
                except Exception:
                    added = 0
            if added > 0:
                total, items = catalog.execute_variant_search(**search_args)

    return {
        "total": total,
        "page": page,
        "page_size": page_size,
        "total_pages": max(1, (total + page_size - 1) // page_size),
        "items": items,
    }


@router.get("/variants/{record_id}")
def get_variant_detail(record_id: int) -> Dict[str, Any]:
    conn = get_db_connection()
    try:
        row = conn.execute("SELECT * FROM variants WHERE id = ?", (record_id,)).fetchone()
    finally:
        conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="Variant record not found")

    record = dict(row)
    pos = record["position"]
    # Neighbouring VCF reads around this position (+/- 2,500 bp)
    record["context_vcf_reads"] = vcf_engine.read_vcf_region(
        chrom=record["chr"],
        start_pos=max(1, pos - 2500),
        end_pos=pos + 2500,
        variants_only=False,
        limit=20,
    )
    return record


@router.get("/vcf/read")
def read_vcf_directly(
    chr: str = Query(..., description="Chromosome (1-22, X, Y, MT)"),
    position: str = Query(..., description="Position or start-end range on GRCh38 (e.g. 111803962 or 111800000-111810000)"),
    variants_only: bool = Query(False, description="Skip homozygous reference blocks (ALT = '.')"),
    limit: int = Query(50, ge=1, le=200),
) -> Dict[str, Any]:
    chrom = normalize_chrom(chr)
    p_start, p_end = catalog.parse_position_filter(position)
    if p_start is None or p_end is None:
        raise HTTPException(status_code=400, detail="Invalid position format. Use an integer (e.g. 111803962) or range (e.g. 111800000-111810000).")

    reads = vcf_engine.read_vcf_region(chrom=chrom, start_pos=p_start, end_pos=p_end, variants_only=variants_only, limit=limit)
    return {
        "chr": chrom,
        "start_position": p_start,
        "end_position": p_end,
        "variants_only": variants_only,
        "count": len(reads),
        "reads": reads,
    }


# ---------------------------------------------------------------------------
# NCBI (ClinVar, dbSNP, Gene, PubMed) endpoints
# ---------------------------------------------------------------------------

_RSID_RE = re.compile(r"^(?:rs)?\d{1,12}$", re.IGNORECASE)
_RCV_RE = re.compile(r"^RCV\d{1,12}(?:\.\d+)?$", re.IGNORECASE)
_VCV_RE = re.compile(r"^(?:VCV)?\d{1,12}(?:\.\d+)?$", re.IGNORECASE)
_CHR_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M)$", re.IGNORECASE)
_GENE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9\-.]{0,29}$")
_SIGNIFICANCES = ("pathogenic", "risk", "drug", "all")


def _clean(v: Optional[str]) -> Optional[str]:
    if v is None:
        return None
    v = v.strip()
    return v or None


@router.get("/ncbi/variant")
def ncbi_variant(
    rsid: Optional[str] = Query(None, description="dbSNP rsID (e.g. rs671)"),
    rcv_id: Optional[str] = Query(None, description="ClinVar RCV accession (e.g. RCV000019449)"),
    vcv_id: Optional[str] = Query(None, description="ClinVar VCV accession or variation id (e.g. VCV000017865)"),
    chr: Optional[str] = Query(None, description="Chromosome (GRCh38)"),
    position: Optional[str] = Query(None, description="GRCh38 position (1-based)"),
    include_publications: bool = Query(True, description="Include linked PubMed publications"),
) -> Dict[str, Any]:
    rsid, rcv_id, vcv_id, chr, position = map(_clean, (rsid, rcv_id, vcv_id, chr, position))
    if rsid and not _RSID_RE.match(rsid):
        raise HTTPException(status_code=400, detail="Invalid rsid (expected e.g. rs671)")
    if rcv_id and not _RCV_RE.match(rcv_id):
        raise HTTPException(status_code=400, detail="Invalid rcv_id (expected e.g. RCV000019449)")
    if vcv_id and not _VCV_RE.match(vcv_id):
        raise HTTPException(status_code=400, detail="Invalid vcv_id (expected e.g. VCV000017865)")
    pos_int: Optional[int] = None
    if position:
        p = position.replace(",", "")
        if ":" in p:
            p = p.split(":", 1)[1]
        if not p.isdigit():
            raise HTTPException(status_code=400, detail="Invalid position (expected a GRCh38 integer coordinate)")
        pos_int = int(p)
    if chr and not _CHR_RE.match(chr):
        raise HTTPException(status_code=400, detail="Invalid chr (expected 1-22, X, Y or MT)")
    has_locus = bool(chr and pos_int)
    if not (rsid or rcv_id or vcv_id or has_locus):
        raise HTTPException(status_code=400, detail="Provide at least one of rsid, rcv_id, vcv_id, or chr + position")
    if rsid and not rsid.lower().startswith("rs"):
        rsid = f"rs{rsid}"
    return ncbi_service.variant_lookup(
        rsid=rsid,
        rcv_id=rcv_id,
        vcv_id=vcv_id,
        chrom=chr if has_locus else None,
        position=pos_int if has_locus else None,
        include_publications=include_publications,
    )


@router.get("/ncbi/gene/{symbol}")
def ncbi_gene(
    symbol: str,
    significance: str = Query("pathogenic", description="pathogenic | risk | drug | all"),
    limit: int = Query(40, ge=1, le=100),
) -> Dict[str, Any]:
    symbol = symbol.strip()
    if not _GENE_RE.match(symbol):
        raise HTTPException(status_code=400, detail="Invalid gene symbol")
    sig = (significance or "pathogenic").strip().lower()
    if sig not in _SIGNIFICANCES:
        raise HTTPException(status_code=400, detail=f"significance must be one of {', '.join(_SIGNIFICANCES)}")
    return ncbi_service.gene_variants(symbol, sig, limit)


@router.get("/ncbi/search")
def ncbi_search(
    term: str = Query(..., min_length=1, max_length=200, description="RCV, VCV, rsID, gene symbol, condition, or free text"),
    limit: int = Query(25, ge=1, le=100),
) -> Dict[str, Any]:
    t = term.strip()
    if not t:
        raise HTTPException(status_code=400, detail="term must not be empty")
    return ncbi_service.search(t, limit)


class NcbiImportRequest(BaseModel):
    variation_ids: List[str]


@router.post("/ncbi/import", dependencies=[Depends(require_app_header)])
def ncbi_import(payload: NcbiImportRequest = Body(...)) -> Dict[str, Any]:
    ids = [str(v).strip() for v in payload.variation_ids if str(v).strip()]
    if not 1 <= len(ids) <= 50:
        raise HTTPException(status_code=400, detail="variation_ids must contain between 1 and 50 ClinVar variation ids")
    return ncbi_service.import_variation_ids(ids)


@router.get("/variants/{record_id}/ncbi")
def variant_ncbi(record_id: int, sync: bool = Query(False, description="Deprecated: use POST /api/variants/{id}/ncbi/sync")) -> Dict[str, Any]:
    if sync:
        # GET must stay side-effect free (an <img src> on any page could trigger it)
        raise HTTPException(status_code=405, detail="Catalog sync moved to POST /api/variants/{id}/ncbi/sync")
    body = ncbi_service.catalog_record_ncbi(record_id, sync=False)
    if body is None:
        raise HTTPException(status_code=404, detail="Variant record not found")
    return body


@router.post("/variants/{record_id}/ncbi/sync", dependencies=[Depends(require_app_header)])
def variant_ncbi_sync(record_id: int) -> Dict[str, Any]:
    body = ncbi_service.catalog_record_ncbi(record_id, sync=True)
    if body is None:
        raise HTTPException(status_code=404, detail="Variant record not found")
    return body
