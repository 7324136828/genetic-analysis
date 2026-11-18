"""
Job orchestration: queue, worker subprocesses, progress tracking, discard / delete,
TTL cleanup and crash recovery.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from app.config import BACKEND_DIR, get_settings
from app.services import job_store
from app.services.job_store import ACTIVE_STATES, COMPLETED, DISCARDED, FAILED, FINAL_STATES, IN_PROGRESS, QUEUED, utcnow
from app.utils.temp_manager import is_valid_job_id, jobs_root, purge_job_dir
from genomescope.workspace import Workspace

STAGES: List[Dict[str, str]] = [
    {"key": "upload", "label": "Stage uploaded VCFs"},
    {"key": "index", "label": "Index gVCF"},
    {"key": "catalog", "label": "Build annotated catalog"},
    {"key": "clinvar", "label": "ClinVar release verification"},
    {"key": "package", "label": "Package results ZIP"},
]
LOG_TAIL_ON_ERROR = 30
SWEEP_INTERVAL_SECONDS = 15 * 60
ORPHAN_MIN_AGE_SECONDS = 6 * 3600


def initial_stages() -> List[Dict[str, Any]]:
    return [dict(s, status="done" if s["key"] == "upload" else "pending", message="") for s in STAGES]


class JobManager:
    def __init__(self) -> None:
        self._executor: Optional[ThreadPoolExecutor] = None
        self._procs: Dict[str, subprocess.Popen] = {}
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._sweeper: Optional[threading.Thread] = None

    # ------------------------------------------------------------------ lifecycle
    def start(self) -> None:
        self._stop.clear()
        self._executor = ThreadPoolExecutor(max_workers=get_settings().max_concurrent_jobs, thread_name_prefix="genomescope-job")
        self.recover_interrupted()
        self.sweep_expired()
        self._sweeper = threading.Thread(target=self._sweep_loop, name="genomescope-sweeper", daemon=True)
        self._sweeper.start()

    def shutdown(self) -> None:
        self._stop.set()
        with self._lock:
            running = dict(self._procs)
        for job_id, proc in running.items():
            job_store.update_if_status(job_id, ACTIVE_STATES, status=FAILED, completed_at=utcnow(),
                                       error_message="Interrupted: the server was stopped while this analysis was running.",
                                       message="Interrupted by server shutdown")
            self._kill(proc)
        if self._executor:
            self._executor.shutdown(wait=False, cancel_futures=True)
            self._executor = None

    def _sweep_loop(self) -> None:
        while not self._stop.wait(SWEEP_INTERVAL_SECONDS):
            try:
                self.sweep_expired()
            except Exception:  # pragma: no cover - never kill the sweeper
                pass

    # ------------------------------------------------------------------ queue
    def submit(self, job_id: str) -> None:
        if self._executor is None:
            raise RuntimeError("job manager is not running")
        self._executor.submit(self._run, job_id)

    def _run(self, job_id: str) -> None:
        job = job_store.get(job_id)
        if not job or job["status"] != QUEUED:
            return
        stages = job["stages"] or initial_stages()
        if not job_store.update_if_status(job_id, (QUEUED,), status=IN_PROGRESS, started_at=utcnow(), stage="index",
                                          progress=5, message="Starting analysis worker", stages=stages):
            return
        ws = Workspace(Path(job["temp_dir"]), job["label"])
        opts = job.get("options") or {}
        meta = {"id": job["id"], "label": job["label"], "created_at": job["created_at"], "options": opts}
        cmd = [
            sys.executable, "-m", "app.worker",
            "--job-dir", str(ws.root),
            "--label", job["label"],
            "--clinvar", str(opts.get("clinvar_mode", "cached")),
            "--job-meta", json.dumps(meta),
        ]
        env = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        result: Optional[Dict[str, Any]] = None
        error: Optional[str] = None
        tail: List[str] = []
        try:
            proc = subprocess.Popen(cmd, cwd=str(BACKEND_DIR), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, encoding="utf-8", errors="replace", bufsize=1, creationflags=flags)
        except OSError as e:
            job_store.update_if_status(job_id, (IN_PROGRESS,), status=FAILED, completed_at=utcnow(),
                                       error_message=f"Could not start the analysis worker: {e}", message="Failed")
            return
        with self._lock:
            self._procs[job_id] = proc
        try:
            ws.work_dir.mkdir(parents=True, exist_ok=True)
            with open(ws.log_path, "a", encoding="utf-8") as log:
                assert proc.stdout is not None
                for line in proc.stdout:
                    line = line.rstrip("\r\n")
                    if line.startswith("@@STAGE "):
                        self._on_stage(job_id, stages, line[8:])
                        continue
                    if line.startswith("@@RESULT "):
                        try:
                            result = json.loads(line[9:])
                        except ValueError:
                            result = {}
                        continue
                    if line.startswith("@@ERROR "):
                        error = line[8:]
                        log.write(f"ERROR: {error}\n")
                        continue
                    log.write(line + "\n")
                    log.flush()
                    tail = (tail + [line])[-LOG_TAIL_ON_ERROR:]
            rc = proc.wait()
        except Exception as e:  # pragma: no cover - log file vanished (discarded)
            rc = proc.poll() if proc.poll() is not None else -1
            error = error or f"{e.__class__.__name__}: {e}"
        finally:
            with self._lock:
                self._procs.pop(job_id, None)

        current = job_store.get(job_id)
        if not current or current["status"] == DISCARDED:
            # discard purged the folder; retry in case the log file was still open then
            purge_job_dir(ws.root)
            return
        zip_path = (result or {}).get("package", {}).get("zip_path")
        if rc == 0 and result is not None and zip_path and Path(zip_path).is_file():
            job_store.update_if_status(
                job_id, (IN_PROGRESS,), status=COMPLETED, progress=100, stage="done", completed_at=utcnow(),
                message="Analysis complete", zip_path=zip_path, result=result, stages=stages,
            )
        else:
            for s in stages:
                if s["status"] == "running":
                    s["status"] = "failed"
            msg = error or (tail[-1] if tail else f"worker exited with code {rc}")
            job_store.update_if_status(job_id, (IN_PROGRESS,), status=FAILED, completed_at=utcnow(),
                                       error_message=msg[:2000], message="Analysis failed", stages=stages)

    def _on_stage(self, job_id: str, stages: List[Dict[str, Any]], payload: str) -> None:
        try:
            ev = json.loads(payload)
        except ValueError:
            return
        for s in stages:
            if s["key"] == ev.get("key"):
                s["status"] = ev.get("status", s["status"])
                s["message"] = ev.get("message", "")
        job_store.update_if_status(job_id, (IN_PROGRESS,), stage=ev.get("key"), progress=int(ev.get("progress", 0)),
                                   message=ev.get("message") or "", stages=stages)

    # ------------------------------------------------------------------ actions
    @staticmethod
    def _kill(proc: subprocess.Popen) -> None:
        if proc.poll() is not None:
            return
        try:
            proc.kill()
            proc.wait(timeout=15)
        except Exception:  # pragma: no cover
            pass

    def discard(self, job_id: str) -> Dict[str, Any]:
        """Abort a queued / running analysis and purge its temp folder (the ledger row stays)."""
        job = job_store.get(job_id)
        if not job:
            raise KeyError(job_id)
        if job["status"] not in ACTIVE_STATES:
            raise ValueError(f"only queued or in-progress analyses can be discarded (status: {job['status']})")
        job_store.update_if_status(job_id, ACTIVE_STATES, status=DISCARDED, completed_at=utcnow(),
                                   message="Discarded; temporary files purged", progress=job.get("progress") or 0)
        with self._lock:
            proc = self._procs.get(job_id)
        if proc is not None:
            self._kill(proc)
        purged = purge_job_dir(job["temp_dir"])
        if not purged:
            job_store.update(job_id, message="Discarded; some temporary files could not be removed yet")
        return job_store.get(job_id)  # type: ignore[return-value]

    def delete(self, job_id: str) -> None:
        """Remove an analysis completely: temp folder and ledger row."""
        job = job_store.get(job_id)
        if not job:
            raise KeyError(job_id)
        if job["status"] in ACTIVE_STATES:
            self.discard(job_id)
        if not purge_job_dir(job["temp_dir"]):
            raise OSError("the analysis folder could not be removed (a file may still be open); try again")
        job_store.delete(job_id)

    def recover_interrupted(self) -> None:
        """After a restart: running jobs cannot be resumed (their worker is gone)."""
        for job in job_store.list_all():
            if job["status"] == IN_PROGRESS:
                job_store.update(job["id"], status=FAILED, completed_at=utcnow(), message="Interrupted",
                                 error_message="Interrupted: the server stopped while this analysis was running. Start it again from New Analysis.")
                purge_job_dir(job["temp_dir"])
            elif job["status"] == QUEUED:
                self.submit(job["id"])

    def sweep_expired(self) -> int:
        """Purge analyses past the TTL and orphaned job folders. Returns folders removed."""
        removed = 0
        ttl = get_settings().job_ttl_hours
        known = set()
        now = datetime.now(timezone.utc)
        for job in job_store.list_all():
            known.add(job["id"])
            if ttl <= 0 or job["status"] not in FINAL_STATES:
                continue
            stamp = job.get("completed_at") or job.get("updated_at") or job.get("created_at")
            try:
                when = datetime.fromisoformat(str(stamp).replace("Z", "+00:00"))
            except ValueError:
                continue
            if now - when > timedelta(hours=ttl):
                if purge_job_dir(job["temp_dir"]):
                    job_store.delete(job["id"])
                    removed += 1
        # Folders without a registry row (crash during an upload, registry deleted ...).
        # A folder is only an orphan once it is old enough not to be an upload in progress.
        for child in jobs_root().iterdir():
            if child.is_dir() and is_valid_job_id(child.name) and child.name not in known:
                try:
                    age = now.timestamp() - child.stat().st_mtime
                except OSError:
                    continue
                if age > ORPHAN_MIN_AGE_SECONDS and purge_job_dir(child):
                    removed += 1
        return removed

    def is_running(self, job_id: str) -> bool:
        with self._lock:
            return job_id in self._procs


manager = JobManager()
