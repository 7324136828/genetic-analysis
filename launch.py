"""Launch the GenomeScope backend (FastAPI/uvicorn) and frontend (Vite) together.

Normally started through run.bat (Windows) or run.sh (macOS / Linux), which run setup
first when needed and then call this script with the project's .venv interpreter.

Ports are chosen dynamically: the first free port in each list of port.json
("backend" and "frontend") is used, and the choice is wired through to Vite
(listen port + /api proxy target) and the backend (CORS origin).

Usage:  python launch.py [--no-browser] [--build]
  --no-browser  do not open the browser when both servers are up
  --build       serve the production build (npm run build + vite preview)
                instead of the dev server

Ctrl+C stops both servers.
"""
import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
import webbrowser
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parent
BACKEND = ROOT / "backend"
FRONTEND = ROOT / "frontend"
PORT_FILE = ROOT / "port.json"
BACKEND_HOST = "127.0.0.1"
DEFAULT_PORTS = {
    "backend": list(range(8000, 8010)),
    "frontend": list(range(5173, 5183)),
}


def load_candidates() -> dict:
    """Read candidate port lists from port.json (falling back to defaults per key)."""
    ports = {k: list(v) for k, v in DEFAULT_PORTS.items()}
    if not PORT_FILE.exists():
        print(f"[launch] {PORT_FILE.name} not found - using default port ranges.")
        return ports
    try:
        data = json.loads(PORT_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        sys.exit(f"[launch] Could not read {PORT_FILE.name}: {e}")
    for key in ports:
        if key not in data:
            continue
        raw = data[key]
        raw = [raw] if isinstance(raw, int) else raw  # allow a single number too
        if (
            not isinstance(raw, list)
            or not raw
            or not all(isinstance(p, int) and not isinstance(p, bool) and 1024 <= p <= 65535 for p in raw)
        ):
            sys.exit(f'[launch] port.json: "{key}" must be a non-empty list of ports between 1024 and 65535.')
        ports[key] = list(dict.fromkeys(raw))  # de-duplicate, keep order
    return ports


def port_in_use(port: int) -> bool:
    """True if something listens on the port (IPv4 or IPv6 loopback) or it cannot be bound."""
    for family, addr in ((socket.AF_INET, "127.0.0.1"), (socket.AF_INET6, "::1")):
        try:
            with socket.socket(family, socket.SOCK_STREAM) as s:
                s.settimeout(0.5)
                if s.connect_ex((addr, port)) == 0:
                    return True
        except OSError:
            pass  # e.g. IPv6 unavailable
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("127.0.0.1", port))
    except OSError:
        return True
    return False


def pick_port(name: str, candidates: List[int], taken: Optional[int] = None) -> int:
    busy = []
    for p in candidates:
        if p == taken:
            continue
        if port_in_use(p):
            busy.append(p)
            continue
        if busy:
            print(f"[launch] {name}: port(s) {', '.join(map(str, busy))} busy -> using {p}")
        return p
    sys.exit(
        f"[launch] No free {name} port among {candidates} (all in use). "
        f"Free one up or add more ports to {PORT_FILE.name}."
    )


def kill_tree(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    else:
        proc.terminate()


def wait_for(url: str, timeout: float, procs) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        if any(p.poll() is not None for p in procs):
            return False
        try:
            with urllib.request.urlopen(url, timeout=2):
                return True
        except Exception:
            time.sleep(1)
    return False


def run_checked(cmd, cwd, label):
    print(f"[launch] {label} ...")
    if subprocess.run(cmd, cwd=cwd).returncode != 0:
        sys.exit(f"[launch] {label} failed.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--build", action="store_true")
    args = ap.parse_args()

    npm = shutil.which("npm")
    if not npm:
        return print("[launch] npm not found. Install Node.js from https://nodejs.org and retry.") or 1

    candidates = load_candidates()
    backend_port = pick_port("backend", candidates["backend"])
    frontend_port = pick_port("frontend", candidates["frontend"], taken=backend_port)

    # Dependencies are installed by setup.py (setup.bat / setup.sh) into .venv
    try:
        import fastapi  # noqa: F401
        import multipart  # noqa: F401  (python-multipart, needed for uploads)
        import uvicorn  # noqa: F401
        import httpx  # noqa: F401
    except ImportError as e:
        return print(
            f"[launch] Backend dependency missing ({e.name}). Run setup.bat (Windows) or ./setup.sh first, "
            "and start with run.bat / run.sh so the project's .venv is used."
        ) or 1
    if not (FRONTEND / "node_modules").exists():
        run_checked([npm, "install", "--no-audit", "--no-fund"], FRONTEND, "Installing frontend dependencies (npm install)")

    # Ports reach both servers through the environment: vite.config.js reads the first two,
    # the backend's CORS policy reads the third (only the chosen frontend origin is allowed).
    env = dict(os.environ)
    env["GENOMESCOPE_BACKEND_PORT"] = str(backend_port)
    env["GENOMESCOPE_FRONTEND_PORT"] = str(frontend_port)
    env["GENOMESCOPE_CORS_ORIGINS"] = f"http://localhost:{frontend_port},http://127.0.0.1:{frontend_port}"

    if args.build:
        print("[launch] Building frontend (npm run build) ...")
        if subprocess.run([npm, "run", "build"], cwd=FRONTEND, env=env).returncode != 0:
            sys.exit("[launch] Building frontend failed.")

    flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    procs = []
    try:
        print(f"[launch] Starting backend on http://{BACKEND_HOST}:{backend_port} ...")
        backend = subprocess.Popen(
            [sys.executable, "-m", "uvicorn", "app.main:app", "--host", BACKEND_HOST, "--port", str(backend_port)],
            cwd=BACKEND,
            env=env,
            creationflags=flags,
        )
        procs.append(backend)

        front_cmd = [npm, "run", "preview"] if args.build else [npm, "run", "dev"]
        print(f"[launch] Starting frontend on http://localhost:{frontend_port} ...")
        frontend = subprocess.Popen(front_cmd, cwd=FRONTEND, env=env, creationflags=flags)
        procs.append(frontend)

        if not wait_for(f"http://{BACKEND_HOST}:{backend_port}/api/health", 180, procs):
            print("[launch] Backend did not become healthy (it may have crashed - see output above).")
            return 1
        if not wait_for(f"http://localhost:{frontend_port}/", 60, procs):
            print("[launch] Frontend did not start (see output above).")
            return 1

        url = f"http://localhost:{frontend_port}/"
        print(
            f"\n[launch] Ready:  UI  {url}\n"
            f"[launch]         API  http://{BACKEND_HOST}:{backend_port}/docs\n"
            f"[launch] Press Ctrl+C to stop both servers.\n",
            flush=True,
        )
        if not args.no_browser:
            webbrowser.open(url)

        while all(p.poll() is None for p in procs):
            time.sleep(1)
        print("[launch] A server exited; shutting down the other one.")
        return 1
    except KeyboardInterrupt:
        print("\n[launch] Stopping servers ...")
        return 0
    finally:
        for p in procs:
            kill_tree(p)


if __name__ == "__main__":
    sys.exit(main())
