"""
Pipeline stages of one analysis. Each stage wraps original genomescope logic and works
entirely inside the job's temp folder (the active workspace).

    index    sparse byte-offset index of the gVCF        (genomescope.vcf_engine)
    catalog  annotated variant catalog in SQLite         (genomescope.indexer)
    clinvar  allele-specific ClinVar release scan        (genomescope.clinvar_scan)
    package  CSV / JSON exports + results ZIP            (this module)
"""

from __future__ import annotations

import csv
import json
import os
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from app.config import APP_VERSION
from app.services.catalog import build_summary
from genomescope import clinvar_release as cr
from genomescope import clinvar_scan, indexer
from genomescope.vcf_engine import engine_for
from genomescope.workspace import Workspace, offline_mode

Say = Callable[[str], None]

CLINVAR_MODES = ("off", "cached", "update")
RESULTS_ZIP_NAME = "genomescope_results.zip"

# Columns exported to CSV (raw_vcf_line is included last: it is the user's own data and
# lets them trace every row back to the VCF).
EXPORT_COLUMNS = [
    "id", "variant_id", "rcv_id", "all_rcv_ids", "chr", "position", "end_position", "gene", "gene_name",
    "ref_allele", "alt_allele", "risk_version", "risk_frequency", "risk_frequency_text", "frequency_source",
    "variant_type", "your_status", "genotype", "alleles_called", "your_data", "is_carried", "condition",
    "classification", "evidence", "evidence_stars", "review_status", "hgvs", "molecular_consequence",
    "dp", "gq", "ad", "qual", "filter_status", "source_file", "clinvar_variation_id", "clinvar_release",
    "clinvar_match", "ncbi_variation_id", "ncbi_last_evaluated", "ncbi_synced_at", "raw_vcf_line",
]


def stage_index(ws: Workspace, say: Say) -> Dict[str, Any]:
    t0 = time.time()
    eng = engine_for(ws)
    eng.ensure_loaded()
    if not eng.snp_vcf_path.exists():
        raise RuntimeError("the SNP/indel gVCF is missing from the job folder")
    if not eng.chrom_ranges:
        raise RuntimeError("no variant records were found in the gVCF")
    say(f"Indexed gVCF: {len(eng.blocks):,} seek blocks across {len(eng.chrom_ranges)} contigs")
    return {"seek_blocks": len(eng.blocks), "contigs": len(eng.chrom_ranges), "elapsed_seconds": round(time.time() - t0, 2)}


def stage_catalog(ws: Workspace, say: Say) -> Dict[str, Any]:
    if offline_mode():
        say("Offline mode: MyVariant.info / mygene.info are not contacted; cataloging structural / copy-number calls only.")
    else:
        say("Querying MyVariant.info / mygene.info for ClinVar annotations and joining them with the gVCF ...")
    indexer.init_db()
    res = indexer.build_initial_database(force=True)
    say(f"Catalog built: {res.get('total_records', 0):,} records, {res.get('carried_variants', 0):,} carried "
        f"({res.get('snp_indel_indexed', 0):,} SNV/indel, {res.get('cnv_sv_indexed', 0):,} CNV/SV) in {res.get('elapsed_seconds')} s")
    if not res.get("total_records"):
        say("WARNING: the catalog is empty. MyVariant.info may be unreachable; the explorer still reads the VCF directly.")
    return res


def stage_clinvar(ws: Workspace, mode: str, say: Say) -> Dict[str, Any]:
    if mode not in CLINVAR_MODES:
        raise ValueError(f"unknown ClinVar mode {mode!r}")
    if mode == "update" and offline_mode():
        say("Offline mode: the ClinVar release is not downloaded; using the cached release if present.")
        mode = "cached"
    if mode == "off":
        say("ClinVar release verification skipped (not requested).")
        return {"status": "skipped", "reason": "not requested"}
    if mode == "cached" and not cr.ClinVarIndex.available():
        say("ClinVar release verification skipped: no local ClinVar release is cached yet. "
            "Choose 'Download / refresh ClinVar release' to fetch it (~460 MB, once).")
        return {"status": "skipped", "reason": "no cached ClinVar release"}
    with _clinvar_cache_lock(say):
        summary_path = clinvar_scan.prepare_index(update=(mode == "update"), progress=say)
    index = cr.ClinVarIndex()
    try:
        summary = clinvar_scan.run_scan(ws.snp_vcf_path, index, summary_gz=summary_path, progress=say)
    finally:
        index.close()
    violations = summary.get("postcondition_violations") or []
    if violations:
        say(f"WARNING: {len(violations)} post-condition violation(s); see clinvar_scan_report.json")
    return {
        "status": "done",
        "clinvar_release": summary.get("clinvar_release"),
        "catalog": summary.get("catalog"),
        "postcondition_violations": len(violations),
    }


class _clinvar_cache_lock:
    """Cross-process lock around the shared ClinVar cache (download / index build)."""

    STALE_SECONDS = 3 * 3600

    def __init__(self, say: Say, timeout: float = 3 * 3600):
        self.path = cr.CLINVAR_DIR / ".genomescope.lock"
        self.say = say
        self.timeout = timeout
        self.fd: Optional[int] = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        deadline = time.time() + self.timeout
        announced = False
        while True:
            try:
                self.fd = os.open(str(self.path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(self.fd, str(os.getpid()).encode())
                return self
            except FileExistsError:
                try:
                    if time.time() - self.path.stat().st_mtime > self.STALE_SECONDS:
                        self.path.unlink(missing_ok=True)
                        continue
                except OSError:
                    continue
                if not announced:
                    self.say("Waiting for another analysis that is updating the shared ClinVar cache ...")
                    announced = True
                if time.time() > deadline:
                    raise TimeoutError("timed out waiting for the ClinVar cache lock")
                time.sleep(2)

    def __exit__(self, *exc):
        if self.fd is not None:
            os.close(self.fd)
        self.path.unlink(missing_ok=True)
        return False


def _export_csv(path: Path, where: str = "") -> int:
    conn = indexer.get_db_connection()
    try:
        cols_present = {r[1] for r in conn.execute("PRAGMA table_info(variants)")}
        cols = [c for c in EXPORT_COLUMNS if c in cols_present]
        rows = conn.execute(
            f"SELECT {', '.join(cols)} FROM variants {where} ORDER BY is_carried DESC, evidence_stars DESC, chr, position"
        )
        n = 0
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for r in rows:
                w.writerow(["" if v is None else v for v in r])
                n += 1
        return n
    finally:
        conn.close()


def stage_package(ws: Workspace, say: Say, job_meta: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Write exports into outputs/ and zip them (plus the pipeline log) into archive/."""
    ws.outputs_dir.mkdir(parents=True, exist_ok=True)
    ws.archive_dir.mkdir(parents=True, exist_ok=True)
    all_rows = _export_csv(ws.outputs_dir / "variants_catalog.csv")
    carried_rows = _export_csv(ws.outputs_dir / "carried_variants.csv", "WHERE is_carried = 1")
    summary = build_summary(ws)
    summary_doc = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": f"GenomeScope {APP_VERSION}",
        "analysis": {k: v for k, v in (job_meta or {}).items() if k in ("id", "label", "created_at", "options")},
        "clinvar_status": clinvar_scan.status(),
        "summary": summary,
        "disclaimer": "Research-grade WGS interpretation; not a clinical diagnosis. Confirm clinically relevant findings with an accredited laboratory.",
    }
    (ws.outputs_dir / "analysis_summary.json").write_text(json.dumps(summary_doc, indent=2), encoding="utf-8")

    zip_path = ws.archive_dir / RESULTS_ZIP_NAME
    tmp = zip_path.with_suffix(".zip.part")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(ws.outputs_dir.rglob("*")):
            if f.is_file():
                zf.write(f, arcname=f.relative_to(ws.outputs_dir).as_posix())
        if ws.log_path.exists():
            zf.write(ws.log_path, arcname="pipeline.log")
    tmp.replace(zip_path)
    say(f"Packaged {all_rows:,} catalog rows ({carried_rows:,} carried) into {zip_path.name}")
    return {
        "zip_path": str(zip_path),
        "catalog_rows": all_rows,
        "carried_rows": carried_rows,
        "pathogenic_or_likely_carried": summary.get("pathogenic_or_likely_carried"),
        "drug_response_carried": summary.get("drug_response_carried"),
        "risk_factor_carried": summary.get("risk_factor_carried"),
    }


def zip_is_stale(ws: Workspace) -> bool:
    """True when the catalog changed (NCBI import / sync in the explorer) after packaging."""
    zip_path = ws.archive_dir / RESULTS_ZIP_NAME
    if not zip_path.exists():
        return True
    try:
        return ws.db_path.stat().st_mtime > zip_path.stat().st_mtime
    except OSError:
        return True
