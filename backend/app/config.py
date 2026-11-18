"""
Runtime settings, read from environment variables (and an optional repo-root ``.env``).

Nothing here is a secret: credentials such as NCBI_API_KEY are read where they are used
(genomescope.ncbi_client) and are never logged or returned by the API.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import List

BACKEND_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = BACKEND_DIR.parent
APP_VERSION = "2.1.0"


def load_dotenv(path: Path) -> None:
    """Minimal .env loader (KEY=VALUE lines, # comments). Existing variables win."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key.startswith("export "):
            key = key[len("export "):].strip()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
            value = value[1:-1]
        if key and key not in os.environ and value != "":
            os.environ[key] = value


load_dotenv(REPO_ROOT / ".env")


def _env_list(name: str, default: List[str]) -> List[str]:
    raw = os.environ.get(name)
    if not raw:
        return list(default)
    return [x.strip() for x in raw.split(",") if x.strip()]


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, "") or default)
    except ValueError:
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, "") or default)
    except ValueError:
        return default


@dataclass(frozen=True)
class Settings:
    # Root of all job folders and of the job registry (jobs.db). Under the system temp
    # directory by default so uploaded genomes never touch the repository.
    jobs_dir: Path
    max_upload_mb: int
    # Completed / failed analyses older than this are purged (0 = keep until deleted).
    job_ttl_hours: float
    max_concurrent_jobs: int
    cors_origins: List[str] = field(default_factory=list)
    allowed_hosts: List[str] = field(default_factory=list)

    @property
    def max_upload_bytes(self) -> int:
        return self.max_upload_mb * 1024 * 1024

    @property
    def jobs_db_path(self) -> Path:
        return self.jobs_dir / "jobs.db"


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    jobs_dir = (os.environ.get("GENOMESCOPE_JOBS_DIR") or "").strip()
    return Settings(
        jobs_dir=Path(jobs_dir).expanduser() if jobs_dir else Path(tempfile.gettempdir()) / "genomescope_jobs",
        max_upload_mb=max(1, _env_int("GENOMESCOPE_MAX_UPLOAD_MB", 4096)),
        job_ttl_hours=max(0.0, _env_float("GENOMESCOPE_JOB_TTL_HOURS", 168.0)),
        max_concurrent_jobs=max(1, _env_int("GENOMESCOPE_MAX_CONCURRENT_JOBS", 1)),
        # The API serves genomes: only the local frontend may call it cross-origin (the
        # Vite dev / preview servers proxy /api, so normally no CORS is needed at all).
        cors_origins=_env_list(
            "GENOMESCOPE_CORS_ORIGINS",
            ["http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:4173", "http://127.0.0.1:4173"],
        ),
        # Host-header allow-list blocks DNS-rebinding pages ("*" = any).
        allowed_hosts=_env_list("GENOMESCOPE_ALLOWED_HOSTS", ["127.0.0.1", "localhost", "::1", "[::1]", "testserver"]),
    )
