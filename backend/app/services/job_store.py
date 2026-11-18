"""
Persistent job registry (SQLite ``jobs.db`` in the jobs root).

States: queued -> in_progress -> completed | failed | discarded
"""

from __future__ import annotations

import json
import sqlite3
import threading
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from app.config import get_settings
from app.utils.temp_manager import jobs_root

QUEUED, IN_PROGRESS, COMPLETED, FAILED, DISCARDED = "queued", "in_progress", "completed", "failed", "discarded"
ACTIVE_STATES = (QUEUED, IN_PROGRESS)
FINAL_STATES = (COMPLETED, FAILED, DISCARDED)

_JSON_COLUMNS = ("files", "options", "stages", "result")
_lock = threading.RLock()

_SCHEMA = """
CREATE TABLE IF NOT EXISTS jobs (
    id TEXT PRIMARY KEY,
    label TEXT NOT NULL,
    status TEXT NOT NULL,
    stage TEXT,
    progress INTEGER NOT NULL DEFAULT 0,
    message TEXT,
    created_at TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT,
    updated_at TEXT NOT NULL,
    error_message TEXT,
    temp_dir TEXT NOT NULL,
    zip_path TEXT,
    file_size INTEGER NOT NULL DEFAULT 0,
    files TEXT,
    options TEXT,
    stages TEXT,
    result TEXT
)
"""


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _connect() -> sqlite3.Connection:
    jobs_root()
    conn = sqlite3.connect(str(get_settings().jobs_db_path), timeout=30, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute(_SCHEMA)
    return conn


def _row_to_job(row: sqlite3.Row) -> Dict[str, Any]:
    job = dict(row)
    for col in _JSON_COLUMNS:
        try:
            job[col] = json.loads(job[col]) if job.get(col) else ([] if col in ("files", "stages") else {})
        except ValueError:
            job[col] = [] if col in ("files", "stages") else {}
    return job


def create(job_id: str, label: str, temp_dir: str, files: List[Dict[str, Any]], file_size: int,
           options: Dict[str, Any], stages: List[Dict[str, Any]]) -> Dict[str, Any]:
    now = utcnow()
    with _lock, _connect() as conn:
        conn.execute(
            "INSERT INTO jobs (id, label, status, stage, progress, message, created_at, updated_at, temp_dir, file_size, files, options, stages) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (job_id, label, QUEUED, "queued", 0, "Waiting for a free worker", now, now, temp_dir, int(file_size),
             json.dumps(files), json.dumps(options), json.dumps(stages)),
        )
    return get(job_id)  # type: ignore[return-value]


def get(job_id: str) -> Optional[Dict[str, Any]]:
    with _lock, _connect() as conn:
        row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    return _row_to_job(row) if row else None


def list_all() -> List[Dict[str, Any]]:
    with _lock, _connect() as conn:
        rows = conn.execute("SELECT * FROM jobs ORDER BY created_at DESC").fetchall()
    return [_row_to_job(r) for r in rows]


def update(job_id: str, **fields: Any) -> Optional[Dict[str, Any]]:
    if not fields:
        return get(job_id)
    fields["updated_at"] = utcnow()
    cols, vals = [], []
    for k, v in fields.items():
        cols.append(f"{k} = ?")
        vals.append(json.dumps(v) if k in _JSON_COLUMNS and v is not None else v)
    with _lock, _connect() as conn:
        conn.execute(f"UPDATE jobs SET {', '.join(cols)} WHERE id = ?", vals + [job_id])
    return get(job_id)


def update_if_status(job_id: str, allowed: tuple, **fields: Any) -> bool:
    """Atomic compare-and-set on status (so a discard can never be overwritten by a finishing worker)."""
    fields["updated_at"] = utcnow()
    cols, vals = [], []
    for k, v in fields.items():
        cols.append(f"{k} = ?")
        vals.append(json.dumps(v) if k in _JSON_COLUMNS and v is not None else v)
    marks = ", ".join("?" for _ in allowed)
    with _lock, _connect() as conn:
        cur = conn.execute(f"UPDATE jobs SET {', '.join(cols)} WHERE id = ? AND status IN ({marks})", vals + [job_id, *allowed])
        return cur.rowcount == 1


def delete(job_id: str) -> None:
    with _lock, _connect() as conn:
        conn.execute("DELETE FROM jobs WHERE id = ?", (job_id,))
