"""Public representation of analysis jobs (internal temp-folder paths are never exposed)."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field

JobStatus = Literal["queued", "in_progress", "completed", "failed", "discarded"]
ClinvarMode = Literal["off", "cached", "update"]


class JobFile(BaseModel):
    slot: str
    label: str
    original_filename: str
    size_bytes: int
    compressed_upload: bool = False
    reference: Optional[str] = None


class JobStage(BaseModel):
    key: str
    label: str
    status: str
    message: str = ""


class Job(BaseModel):
    id: str
    label: str
    status: JobStatus
    stage: Optional[str] = None
    progress: int = 0
    message: Optional[str] = None
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    updated_at: str
    elapsed_seconds: Optional[float] = None
    error_message: Optional[str] = None
    file_size: int = 0
    files: List[JobFile] = Field(default_factory=list)
    options: Dict[str, Any] = Field(default_factory=dict)
    stages: List[JobStage] = Field(default_factory=list)
    result: Dict[str, Any] = Field(default_factory=dict)
    can_open: bool = False
    can_download: bool = False


def _parse(ts: Optional[str]) -> Optional[datetime]:
    if not ts:
        return None
    try:
        return datetime.fromisoformat(str(ts).replace("Z", "+00:00"))
    except ValueError:
        return None


def to_public(job: Dict[str, Any]) -> Job:
    start = _parse(job.get("started_at")) or _parse(job.get("created_at"))
    end = _parse(job.get("completed_at")) or (datetime.now(timezone.utc) if job.get("status") in ("queued", "in_progress") else None)
    elapsed = round((end - start).total_seconds(), 1) if (start and end) else None
    result = dict(job.get("result") or {})
    pkg = dict(result.get("package") or {})
    pkg.pop("zip_path", None)  # internal path
    if pkg:
        result["package"] = pkg
    completed = job.get("status") == "completed"
    return Job(
        id=job["id"],
        label=job["label"],
        status=job["status"],
        stage=job.get("stage"),
        progress=int(job.get("progress") or 0),
        message=job.get("message"),
        created_at=job["created_at"],
        started_at=job.get("started_at"),
        completed_at=job.get("completed_at"),
        updated_at=job["updated_at"],
        elapsed_seconds=elapsed,
        error_message=job.get("error_message"),
        file_size=int(job.get("file_size") or 0),
        files=[JobFile(**{k: f.get(k) for k in JobFile.model_fields if k in f}) for f in (job.get("files") or [])],
        options=job.get("options") or {},
        stages=[JobStage(**{k: s.get(k, "") for k in JobStage.model_fields}) for s in (job.get("stages") or [])],
        result=result,
        can_open=completed,
        can_download=completed,
    )
