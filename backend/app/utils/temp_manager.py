"""
System temp-folder lifecycle for analysis jobs.

Layout (see genomescope.workspace)::

    <system temp>/genomescope_jobs/          (GENOMESCOPE_JOBS_DIR), mode 0700
      jobs.db                                job registry
      <job uuid>/{inputs,work,outputs,archive}

Uploaded genomes, intermediate files and outputs never touch the repository. Every
delete is guarded so it can only remove a direct child folder of the jobs root.
"""

from __future__ import annotations

import os
import shutil
import stat
import time
import uuid
from pathlib import Path
from typing import Optional

from app.config import get_settings
from genomescope.workspace import Workspace

_PRIVATE_DIR_MODE = 0o700


def _make_private_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    if os.name != "nt":  # Windows temp dirs are already per-user
        try:
            os.chmod(path, _PRIVATE_DIR_MODE)
        except OSError:
            pass


def jobs_root() -> Path:
    root = get_settings().jobs_dir
    _make_private_dir(root)
    return root


def new_job_id() -> str:
    return str(uuid.uuid4())


def is_valid_job_id(job_id: str) -> bool:
    try:
        return str(uuid.UUID(str(job_id))) == str(job_id)
    except (ValueError, TypeError):
        return False


def job_dir(job_id: str) -> Path:
    if not is_valid_job_id(job_id):
        raise ValueError("invalid job id")
    return jobs_root() / job_id


def create_job_workspace(job_id: str, sample_label: str) -> Workspace:
    """Create <jobs root>/<job_id>/{inputs,work,outputs,archive} (private to the user)."""
    root = job_dir(job_id)
    _make_private_dir(root)
    ws = Workspace(root, sample_label).ensure_dirs()
    for d in (ws.inputs_dir, ws.work_dir, ws.outputs_dir, ws.archive_dir):
        _make_private_dir(d)
    return ws


def _inside_jobs_root(path: Path) -> bool:
    try:
        resolved = Path(path).resolve()
        root = jobs_root().resolve()
    except OSError:
        return False
    return resolved.parent == root and is_valid_job_id(resolved.name)


def _on_rm_error(func, path, _exc):  # pragma: no cover - Windows read-only files
    try:
        os.chmod(path, stat.S_IWRITE)
        func(path)
    except OSError:
        pass


def purge_job_dir(path: Optional[str | Path], attempts: int = 5) -> bool:
    """Delete one job folder (with retries: on Windows a just-killed worker may still hold
    a file handle for a moment). Returns True when the folder is gone."""
    if not path:
        return True
    p = Path(path)
    if not p.exists():
        return True
    if not _inside_jobs_root(p):
        raise ValueError(f"refusing to delete a folder outside the jobs root: {p}")
    from genomescope.vcf_engine import forget_engine

    forget_engine(p)
    for i in range(attempts):
        shutil.rmtree(p, onerror=_on_rm_error)
        if not p.exists():
            return True
        time.sleep(0.3 * (i + 1))
    return not p.exists()


def folder_size(path: Path) -> int:
    total = 0
    for f in Path(path).rglob("*"):
        try:
            if f.is_file():
                total += f.stat().st_size
        except OSError:
            pass
    return total
