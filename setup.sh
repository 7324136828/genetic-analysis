#!/usr/bin/env bash
# GenomeScope setup (macOS / Linux): creates .venv, installs backend + frontend dependencies.
set -euo pipefail
cd "$(dirname "$0")"

if command -v python3 >/dev/null 2>&1; then
    PY=python3
elif command -v python >/dev/null 2>&1; then
    PY=python
else
    echo "[SETUP.SH] Python 3.10+ was not found on PATH." >&2
    exit 1
fi

echo "[SETUP.SH] Dispatching setup.py with $PY ..."
"$PY" setup.py "$@"
echo "[SETUP.SH] Setup complete. Start the app with ./run.sh"
