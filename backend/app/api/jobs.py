"""
Analysis jobs: upload VCFs, follow progress, history, discard / delete, results ZIP.

    POST   /api/jobs                    multipart upload -> new job (202)
    GET    /api/jobs                    history (newest first)
    GET    /api/jobs/{id}               one job
    GET    /api/jobs/{id}/logs          tail of the pipeline log
    POST   /api/jobs/{id}/discard       abort a queued / running job and purge its temp folder
    DELETE /api/jobs/{id}               remove a job and all its data
    GET    /api/jobs/{id}/download-zip  results archive
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, File, Form, HTTPException, Query, UploadFile, status
from fastapi.responses import FileResponse, Response

from app.api.deps import require_app_header
from app.config import get_settings
from app.schemas.job import ClinvarMode, Job, to_public
from app.services import analysis, job_store, staging
from app.services.pipeline import initial_stages, manager
from app.utils.temp_manager import create_job_workspace, is_valid_job_id, new_job_id, purge_job_dir
from genomescope import clinvar_release as cr
from genomescope import workspace

router = APIRouter(prefix="/api/jobs", tags=["jobs"])

MAX_LABEL = 80


def _get_job_or_404(job_id: str) -> Dict[str, Any]:
    if not is_valid_job_id(job_id):
        raise HTTPException(status_code=400, detail="Invalid analysis id")
    job = job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Analysis not found")
    return job


def _clean_label(label: Optional[str], job_id: str) -> str:
    text = re.sub(r"[\x00-\x1f\x7f]", "", (label or "")).strip()[:MAX_LABEL]
    return text or f"Analysis {job_id[:8]}"


@router.post("", status_code=status.HTTP_202_ACCEPTED, response_model=Job, dependencies=[Depends(require_app_header)])
def create_job(
    snp_indel: UploadFile = File(..., description="SNP / indel gVCF (.vcf or .vcf.gz), required"),
    sv: Optional[UploadFile] = File(None, description="Structural variants VCF (optional)"),
    cnv: Optional[UploadFile] = File(None, description="Copy number variants VCF (optional)"),
    label: Optional[str] = Form(None, description="Display name for this analysis"),
    clinvar_mode: ClinvarMode = Form("cached", description="off | cached (use local release if present) | update (download / refresh, ~460 MB once)"),
) -> Job:
    settings = get_settings()
    uploads = {"snp_indel": snp_indel, "sv": sv, "cnv": cnv}
    uploads = {k: v for k, v in uploads.items() if v is not None and (v.filename or v.size)}
    if "snp_indel" not in uploads:
        raise HTTPException(status_code=400, detail="The SNP / indel gVCF is required")
    total = sum(int(u.size or 0) for u in uploads.values())
    if total > settings.max_upload_bytes:
        raise HTTPException(status_code=413, detail=f"Upload too large ({total / 1e6:,.0f} MB > {settings.max_upload_mb:,} MB limit)")

    job_id = new_job_id()
    clean_label = _clean_label(label, job_id)
    ws = create_job_workspace(job_id, clean_label)
    try:
        files: List[Dict[str, Any]] = []
        for slot, upload in uploads.items():
            upload.file.seek(0)
            files.append(staging.stage_file(ws, slot, upload.file, upload.filename))
        staging.write_manifest(ws, files)
    except staging.StagingError as e:
        purge_job_dir(ws.root)
        raise HTTPException(status_code=400, detail=str(e))
    except Exception:
        purge_job_dir(ws.root)
        raise
    finally:
        for u in uploads.values():
            u.file.close()

    job = job_store.create(
        job_id=job_id,
        label=clean_label,
        temp_dir=str(ws.root),
        files=files,
        file_size=sum(f["size_bytes"] for f in files),
        options={"clinvar_mode": clinvar_mode},
        stages=initial_stages(),
    )
    manager.submit(job_id)
    return to_public(job)


@router.get("", response_model=List[Job])
def list_jobs() -> List[Job]:
    return [to_public(j) for j in job_store.list_all()]


@router.get("/{job_id}", response_model=Job)
def get_job(job_id: str) -> Job:
    return to_public(_get_job_or_404(job_id))


@router.get("/{job_id}/logs")
def get_job_logs(job_id: str, tail: int = Query(200, ge=1, le=5000)) -> Dict[str, Any]:
    job = _get_job_or_404(job_id)
    log_path = workspace.Workspace(Path(job["temp_dir"])).log_path
    lines: List[str] = []
    if log_path.exists():
        try:
            with open(log_path, "r", encoding="utf-8", errors="replace") as f:
                lines = f.read().splitlines()[-tail:]
        except OSError:
            lines = []
    return {"job_id": job_id, "status": job["status"], "lines": lines}


@router.post("/{job_id}/discard", response_model=Job, dependencies=[Depends(require_app_header)])
def discard_job(job_id: str) -> Job:
    _get_job_or_404(job_id)
    try:
        return to_public(manager.discard(job_id))
    except ValueError as e:
        raise HTTPException(status_code=409, detail=str(e))


@router.delete("/{job_id}", status_code=status.HTTP_204_NO_CONTENT, response_class=Response, dependencies=[Depends(require_app_header)])
def delete_job(job_id: str) -> Response:
    _get_job_or_404(job_id)
    try:
        manager.delete(job_id)
    except OSError as e:
        raise HTTPException(status_code=409, detail=str(e))
    return Response(status_code=status.HTTP_204_NO_CONTENT)


@router.get("/{job_id}/download-zip")
def download_zip(job_id: str) -> FileResponse:
    job = _get_job_or_404(job_id)
    if job["status"] != job_store.COMPLETED:
        raise HTTPException(status_code=400, detail="Analysis is not completed")
    ws = workspace.Workspace(Path(job["temp_dir"]), job["label"])
    if not ws.db_path.exists():
        raise HTTPException(status_code=410, detail="The data of this analysis has been purged")
    if analysis.zip_is_stale(ws):
        # The catalog changed after packaging (NCBI import / sync in the explorer): rebuild.
        with workspace.activate(ws):
            analysis.stage_package(ws, lambda _m: None, {"id": job["id"], "label": job["label"], "created_at": job["created_at"], "options": job.get("options")})
    zip_path = ws.archive_dir / analysis.RESULTS_ZIP_NAME
    if not zip_path.is_file():
        raise HTTPException(status_code=404, detail="Results archive not found")
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", job["label"]).strip("_") or "analysis"
    return FileResponse(path=zip_path, filename=f"{safe}_genomescope_results.zip", media_type="application/zip")


meta_router = APIRouter(prefix="/api", tags=["meta"])


@meta_router.get("/config")
def get_config() -> Dict[str, Any]:
    """Limits and shared-cache state the upload screen shows."""
    s = get_settings()
    index_meta: Optional[Dict[str, Any]] = None
    if cr.ClinVarIndex.available():
        try:
            ix = cr.ClinVarIndex()
            try:
                m = ix.meta()
                index_meta = {"release_date": m.get("release_date"), "rows_indexed": m.get("rows_indexed")}
            finally:
                ix.close()
        except cr.ClinVarError:
            index_meta = None
    return {
        "max_upload_mb": s.max_upload_mb,
        "job_ttl_hours": s.job_ttl_hours,
        "max_concurrent_jobs": s.max_concurrent_jobs,
        "clinvar_cache": index_meta,
        "upload_slots": [{"slot": k, "label": v["label"], "required": k in staging.REQUIRED_SLOTS} for k, v in staging.SLOTS.items()],
    }
