#!/usr/bin/env bash
# Start GenomeScope: backend (FastAPI) + frontend (Vite) together.
#   ./run.sh              development servers
#   ./run.sh --build      production build served by vite preview
#   ./run.sh --no-browser do not open the browser
# Ports: first free port of each list in port.json. Ctrl+C stops both servers.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -x ".venv/bin/python" ]; then
    echo "[RUN.SH] Virtual environment not found. Running setup first ..."
    ./setup.sh
fi
if [ ! -d "frontend/node_modules" ]; then
    echo "[RUN.SH] Frontend dependencies missing. Running setup first ..."
    ./setup.sh --skip-backend
fi

# launch.py supervises both servers and stops them on Ctrl+C / SIGTERM
exec .venv/bin/python launch.py "$@"
