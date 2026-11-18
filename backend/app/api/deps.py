"""
Request dependencies: analysis selection (workspace activation) and CSRF protection.

The frontend sends the selected analysis id in the ``X-GenomeScope-Job`` header (or
``?job=`` for plain links). The dependency resolves it to the job's temp folder and
activates that workspace for the rest of the request, so every genomescope module reads
the right genome and catalog.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from fastapi import Header, HTTPException, Query

from app.services import job_store
from app.utils.temp_manager import is_valid_job_id
from genomescope import workspace
from genomescope.workspace import Workspace

# State-changing routes require this header: cross-site pages cannot send it without a
# CORS preflight (which is refused), and plain <form>/<img> requests cannot set it.
CSRF_HEADER_VALUE = "GenomeScope"
JOB_HEADER = "X-GenomeScope-Job"


def require_app_header(x_requested_with: Optional[str] = Header(None)) -> None:
    if (x_requested_with or "").strip() != CSRF_HEADER_VALUE:
        raise HTTPException(status_code=403, detail=f"Missing X-Requested-With: {CSRF_HEADER_VALUE} header")


def _resolve(job_id: Optional[str], required: bool) -> Optional[Workspace]:
    job_id = (job_id or "").strip()
    if not job_id:
        if required:
            raise HTTPException(status_code=409, detail="No analysis selected. Open an analysis from the History screen or start a new one.")
        return None
    if not is_valid_job_id(job_id):
        raise HTTPException(status_code=400, detail="Invalid analysis id")
    job = job_store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Analysis not found (it may have been deleted or expired)")
    if job["status"] != job_store.COMPLETED:
        raise HTTPException(status_code=409, detail=f"Analysis is not ready yet (status: {job['status']})")
    ws = Workspace(Path(job["temp_dir"]), job["label"])
    if not ws.db_path.exists():
        raise HTTPException(status_code=410, detail="The data of this analysis has been purged from the temp folder")
    return ws


# Async on purpose: FastAPI awaits it in the request's own task, so the context variable
# set here is inherited by the (sync) endpoint that runs next in the thread pool.
async def use_analysis(
    x_genomescope_job: Optional[str] = Header(None),
    job: Optional[str] = Query(None, include_in_schema=False),
) -> Workspace:
    ws = _resolve(x_genomescope_job or job, required=True)
    workspace.set_current(ws)
    return ws  # type: ignore[return-value]


async def optional_analysis(
    x_genomescope_job: Optional[str] = Header(None),
    job: Optional[str] = Query(None, include_in_schema=False),
) -> Optional[Workspace]:
    ws = _resolve(x_genomescope_job or job, required=False)
    workspace.set_current(ws)
    return ws
