"""
Per-analysis workspace: where one uploaded genome and everything derived from it lives.

Every analysis (job) runs in an isolated folder under the system temp directory
(see ``app.utils.temp_manager``)::

    <job_dir>/
      inputs/    uploaded VCFs, stored under canonical names (SNP_VCF_NAME, ...)
                 plus manifest.json (original names, sizes, reference)
      work/      SQLite catalog, sparse gVCF index, catalog backups, pipeline log
      outputs/   reports and exports (ClinVar scan report, CSV exports, summary)
      archive/   the downloadable results ZIP

The active workspace is held in a ``contextvars.ContextVar`` so that concurrent
requests / worker threads each see their own genome. Modules never hard-code paths:
they call :func:`current` (or the helpers in ``indexer`` / ``vcf_engine``).
"""

from __future__ import annotations

import contextvars
import json
import os
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, Optional

# Canonical file names inside <job_dir>/inputs. Uploaded files are renamed to these so
# that no user-supplied file name (which often embeds a kit / sample id) ends up in the
# catalog's source_file column, in reports or in log lines.
SNP_VCF_NAME = "snp-indel.genome.vcf"
SV_VCF_NAME = "sv.vcf"
CNV_VCF_NAME = "cnv.vcf"
MANIFEST_NAME = "manifest.json"

DEFAULT_REFERENCE = "GRCh38"
DEFAULT_SAMPLE_LABEL = "Sample"


def offline_mode() -> bool:
    """GENOMESCOPE_OFFLINE=1: never contact MyVariant.info, mygene.info, NCBI or the ClinVar FTP.

    No genome coordinate or identifier leaves the machine; the catalog then holds only the
    structural / copy-number rows and direct VCF reads (plus a locally cached ClinVar scan).
    """
    return (os.environ.get("GENOMESCOPE_OFFLINE") or "").strip().lower() in ("1", "true", "yes", "on")


class NoActiveWorkspace(RuntimeError):
    """Raised when genome data is accessed without an active analysis workspace."""


@dataclass(frozen=True)
class Workspace:
    root: Path
    sample_label: str = DEFAULT_SAMPLE_LABEL

    # --- folders ---------------------------------------------------------
    @property
    def inputs_dir(self) -> Path:
        return self.root / "inputs"

    @property
    def work_dir(self) -> Path:
        return self.root / "work"

    @property
    def outputs_dir(self) -> Path:
        return self.root / "outputs"

    @property
    def archive_dir(self) -> Path:
        return self.root / "archive"

    # --- files -----------------------------------------------------------
    @property
    def snp_vcf_path(self) -> Path:
        return self.inputs_dir / SNP_VCF_NAME

    @property
    def sv_vcf_path(self) -> Path:
        return self.inputs_dir / SV_VCF_NAME

    @property
    def cnv_vcf_path(self) -> Path:
        return self.inputs_dir / CNV_VCF_NAME

    @property
    def manifest_path(self) -> Path:
        return self.inputs_dir / MANIFEST_NAME

    @property
    def db_path(self) -> Path:
        return self.work_dir / "genetic_variants.db"

    @property
    def sparse_index_path(self) -> Path:
        return self.work_dir / "vcf_sparse_index.json"

    @property
    def backup_dir(self) -> Path:
        return self.work_dir / "backups"

    @property
    def log_path(self) -> Path:
        return self.work_dir / "pipeline.log"

    @property
    def scan_report_path(self) -> Path:
        return self.outputs_dir / "clinvar_scan_report.json"

    @property
    def scan_carried_path(self) -> Path:
        return self.outputs_dir / "clinvar_scan_carried.tsv"

    # --- helpers ---------------------------------------------------------
    def ensure_dirs(self) -> "Workspace":
        for d in (self.inputs_dir, self.work_dir, self.outputs_dir, self.archive_dir):
            d.mkdir(parents=True, exist_ok=True)
        return self

    def manifest(self) -> Dict[str, Any]:
        """Upload manifest written at staging time ({} when absent / unreadable)."""
        try:
            return json.loads(self.manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}

    @property
    def reference(self) -> str:
        return str(self.manifest().get("reference") or DEFAULT_REFERENCE)


_current: contextvars.ContextVar[Optional[Workspace]] = contextvars.ContextVar("genomescope_workspace", default=None)


def get() -> Optional[Workspace]:
    """The active workspace, or None."""
    return _current.get()


def current() -> Workspace:
    """The active workspace; raises NoActiveWorkspace when no analysis is selected."""
    ws = _current.get()
    if ws is None:
        raise NoActiveWorkspace("No analysis workspace is active. Select or create an analysis first.")
    return ws


def set_current(ws: Optional[Workspace]) -> None:
    """Set the workspace for the rest of the current context (e.g. one HTTP request task)."""
    _current.set(ws)


@contextmanager
def activate(ws: Workspace) -> Iterator[Workspace]:
    """Activate `ws` for the duration of the with-block (worker threads, CLI, tests)."""
    token = _current.set(ws)
    try:
        yield ws
    finally:
        _current.reset(token)
