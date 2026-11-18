#!/usr/bin/env python3
"""
Cross-platform setup for GenomeScope (invoked by setup.bat / setup.sh).

  1. checks Python >= 3.10 and Node.js / npm
  2. creates the virtual environment .venv and upgrades pip
  3. installs backend/requirements.txt into .venv
  4. installs the frontend dependencies (npm ci, or npm install without a lockfile)
  5. creates .env from .env.example when missing

Usage:  python setup.py [--skip-frontend] [--skip-backend]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import venv
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
BACKEND_DIR = ROOT_DIR / "backend"
FRONTEND_DIR = ROOT_DIR / "frontend"
VENV_DIR = ROOT_DIR / ".venv"
MIN_PYTHON = (3, 10)
MIN_NODE_MAJOR = 18


def log(msg: str) -> None:
    print(f"\n[SETUP] {msg}", flush=True)


def fail(msg: str) -> "None":
    print(f"\n[SETUP] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def venv_python() -> Path:
    return VENV_DIR / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def run(cmd, cwd: Path = ROOT_DIR) -> None:
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    rc = subprocess.call([str(c) for c in cmd], cwd=str(cwd))
    if rc != 0:
        fail(f"command failed with exit code {rc}: {' '.join(str(c) for c in cmd)}")


def check_prerequisites(need_node: bool) -> str:
    log("Checking prerequisites ...")
    if sys.version_info < MIN_PYTHON:
        fail(f"Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+ is required (found {sys.version.split()[0]}).")
    print(f"  Python {sys.version.split()[0]}  OK")
    npm = shutil.which("npm") or ""
    if need_node:
        node = shutil.which("node")
        if not node or not npm:
            fail("Node.js and npm are required. Install Node.js LTS from https://nodejs.org and re-run setup.")
        out = subprocess.run([node, "--version"], capture_output=True, text=True).stdout.strip()
        try:
            major = int(out.lstrip("v").split(".")[0])
        except ValueError:
            major = 0
        if major < MIN_NODE_MAJOR:
            fail(f"Node.js {MIN_NODE_MAJOR}+ is required (found {out or 'unknown'}).")
        print(f"  Node.js {out}  OK")
    return npm


def create_virtualenv() -> Path:
    py = venv_python()
    if py.exists():
        log(f"Using existing virtual environment {VENV_DIR.name}/")
    else:
        log(f"Creating virtual environment {VENV_DIR.name}/ ...")
        venv.EnvBuilder(with_pip=True, clear=False, upgrade_deps=False).create(VENV_DIR)
    if not py.exists():
        fail(f"virtual environment interpreter not found at {py}")
    return py


def install_backend(py: Path) -> None:
    log("Installing backend dependencies ...")
    run([py, "-m", "pip", "install", "--upgrade", "pip", "--disable-pip-version-check", "-q"])
    run([py, "-m", "pip", "install", "-r", BACKEND_DIR / "requirements.txt", "--disable-pip-version-check", "-q"])


def install_frontend(npm: str) -> None:
    log("Installing frontend dependencies ...")
    has_lock = (FRONTEND_DIR / "package-lock.json").exists()
    run([npm, "ci" if has_lock else "install", "--no-audit", "--no-fund"], cwd=FRONTEND_DIR)


def setup_env() -> None:
    log("Environment configuration ...")
    example, target = ROOT_DIR / ".env.example", ROOT_DIR / ".env"
    if example.exists() and not target.exists():
        shutil.copyfile(example, target)
        print("  Created .env from .env.example (edit it to set NCBI_API_KEY etc.; it is git-ignored).")
    else:
        print("  .env already present" if target.exists() else "  no .env.example found; skipping")


def main() -> None:
    ap = argparse.ArgumentParser(description="Set up the GenomeScope development / runtime environment.")
    ap.add_argument("--skip-frontend", action="store_true", help="do not install the frontend dependencies")
    ap.add_argument("--skip-backend", action="store_true", help="do not create .venv / install Python packages")
    args = ap.parse_args()

    npm = check_prerequisites(need_node=not args.skip_frontend)
    if not args.skip_backend:
        install_backend(create_virtualenv())
    if not args.skip_frontend:
        install_frontend(npm)
    setup_env()
    log("Setup completed. Start the app with run.bat (Windows) or ./run.sh (macOS / Linux).")


if __name__ == "__main__":
    main()
