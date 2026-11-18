# GenomeScope - AGENTS.md

See `docs/ARCHITECTURE.md` for the architecture, job lifecycle, VCF indexing, ClinVar scan and API; `docs/NCBI_INTEGRATION.md` and `docs/FRONTEND.md` for details.

## Ground rules
- **Never commit genome data or anything derived from it** (VCF / BAM / CRAM, catalog `.db` files, scan reports, result ZIPs) and never put real sample IDs, kit IDs, genotypes or findings in code, tests, docs or commit messages. Tests use the synthetic genome in `backend/tests/synthetic.py` only.
- Never commit `.env` or credentials; configuration goes through environment variables (`.env.example` documents them). Record secret findings in `secrets.md` without values.
- Analyses run in the system temp folder (`genomescope_jobs`); code must not write user data inside the repository.

## Commands
- Setup: `setup.bat` / `./setup.sh` (creates `.venv`, installs backend + frontend deps, creates `.env`)
- Run both servers: `run.bat` / `./run.sh` (`--build`, `--no-browser`; ports from `port.json`)
- Backend tests (offline): from `backend/`, `../.venv/Scripts/python -m pytest` (POSIX `../.venv/bin/python -m pytest`)
- Live NCBI tests: `RUN_LIVE_NCBI=1 python -m pytest tests/test_ncbi.py -k live` (PowerShell: `$env:RUN_LIVE_NCBI='1'; python -m pytest tests/test_ncbi.py -k live`)
- Backend only (from `backend/`): `python -m uvicorn app.main:app --host 127.0.0.1 --port 8000`
- ClinVar scan for an existing analysis folder (from `backend/`): `python -m genomescope.clinvar_scan --workspace <job_dir> [--update] [--dry-run]`
- Frontend: `cd frontend && npm run dev` / `npm run build`
- No external network: set `GENOMESCOPE_OFFLINE=1`
