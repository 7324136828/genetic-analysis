"""
Analysis worker (one subprocess per job, so a job can be discarded by terminating it).

    python -m app.worker --job-dir <temp job folder> --label <label> --clinvar off|cached|update

Protocol on stdout (read by app.services.pipeline):
    @@STAGE {"key": ..., "status": "running|done|skipped|failed", "progress": 0-100, "message": ...}
    @@RESULT {...}            final result, printed once on success
    @@ERROR <message>         fatal error (exit code 1)
Every other line is a log line. Exit code 0 = success.
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path
from typing import Any, Dict

from app.services import analysis
from genomescope import workspace

# progress range of each stage (start, end)
STAGE_RANGES = {"index": (5, 20), "catalog": (20, 55), "clinvar": (55, 90), "package": (90, 100)}


def emit(tag: str, payload: Any) -> None:
    body = payload if isinstance(payload, str) else json.dumps(payload, default=str)
    print(f"@@{tag} {body}", flush=True)


def stage(key: str, status: str, message: str = "", progress: int | None = None) -> None:
    lo, hi = STAGE_RANGES[key]
    if progress is None:
        progress = lo if status == "running" else hi
    emit("STAGE", {"key": key, "status": status, "progress": progress, "message": message})


def log(msg: str) -> None:
    print(msg, flush=True)


def run(job_dir: Path, label: str, clinvar_mode: str, job_meta: Dict[str, Any]) -> Dict[str, Any]:
    ws = workspace.Workspace(job_dir, label)
    result: Dict[str, Any] = {}
    with workspace.activate(ws):
        stage("index", "running", "Indexing the gVCF (byte-offset index)")
        result["index"] = analysis.stage_index(ws, log)
        stage("index", "done", "gVCF indexed")

        stage("catalog", "running", "Building the annotated variant catalog")
        result["catalog"] = analysis.stage_catalog(ws, log)
        stage("catalog", "done", f"{result['catalog'].get('total_records', 0):,} catalog records")

        stage("clinvar", "running", "Verifying carried alleles against a local ClinVar release" if clinvar_mode != "off" else "Skipped")
        result["clinvar"] = analysis.stage_clinvar(ws, clinvar_mode, log)
        if result["clinvar"]["status"] == "skipped":
            stage("clinvar", "skipped", result["clinvar"].get("reason", "skipped"))
        else:
            stage("clinvar", "done", f"Verified against ClinVar release {result['clinvar'].get('clinvar_release')}")

        stage("package", "running", "Exporting CSV / JSON and packaging the results ZIP")
        result["package"] = analysis.stage_package(ws, log, job_meta)
        stage("package", "done", "Results ZIP ready")
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--job-dir", type=Path, required=True)
    ap.add_argument("--label", default=workspace.DEFAULT_SAMPLE_LABEL)
    ap.add_argument("--clinvar", choices=analysis.CLINVAR_MODES, default="cached")
    ap.add_argument("--job-meta", default="{}", help="JSON with id / label / created_at / options for the summary")
    args = ap.parse_args(argv)
    try:
        meta = json.loads(args.job_meta)
    except ValueError:
        meta = {}
    try:
        result = run(args.job_dir, args.label, args.clinvar, meta)
    except Exception as e:  # noqa: BLE001 - report every failure to the parent
        traceback.print_exc(file=sys.stdout)
        emit("ERROR", f"{e.__class__.__name__}: {e}")
        return 1
    emit("RESULT", result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
