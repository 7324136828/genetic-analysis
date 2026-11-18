# genetic-analysis — GenomeScope

GenomeScope is a local, full-stack application for exploring your own whole-genome sequencing (WGS) data. Upload a single-sample **GRCh38** gVCF (and optionally its structural-variant and copy-number VCFs). GenomeScope then:

1. indexes the gVCF for millisecond lookups at any `chr:position`,
2. builds an annotated clinical variant catalog (ClinVar / dbSNP / gnomAD via MyVariant.info, genes via mygene.info),
3. optionally verifies every carried allele against a local, MD5-verified **ClinVar release** (allele-specific, not just position-matched),
4. packages the results (CSV catalog, carried variants, JSON summary, ClinVar reports) into a downloadable ZIP,

and lets you explore the result: an 11-field search across the catalog, raw VCF reads, and live **NCBI** ClinVar / dbSNP / Gene / PubMed lookups joined with your own genotype.

> Research-grade interpretation of sequencing data — not a medical diagnosis. Confirm clinically relevant findings with an accredited laboratory.

## Quick start

Requirements: **Python 3.10+** and **Node.js 18+** (with npm).

```bash
# Windows
setup.bat
run.bat

# macOS / Linux
./setup.sh
./run.sh
```

`setup` creates `.venv/`, installs the backend and frontend dependencies, and creates `.env` from `.env.example`. `run` starts the backend and frontend on the first free ports listed in `port.json`, opens the browser, and stops both on Ctrl+C. (`run` also runs setup automatically when needed.) Options: `--build` serves the production build, `--no-browser` skips opening the browser.

Then:
1. **New Analysis**: add your SNP & Indel gVCF (`.vcf` or `.vcf.gz`) by clicking, dragging, or pasting with Ctrl+V. Add the SV / CNV VCFs if you have them, choose the ClinVar option, and click **Start analysis**.
2. Follow the stages live (upload → index → catalog → ClinVar → package) with the log viewer; **Discard** stops the analysis and deletes its files.
3. **Open in Explorer** when it completes, or **Download ZIP**. Every analysis is listed under **History**.

## Privacy & data handling

- **Your files stay on your machine.** Each analysis runs in its own folder under the system temp directory (`<temp>/genomescope_jobs/<uuid>/`, owner-only on macOS / Linux), never inside the repository. Uploaded files are renamed to canonical names (`snp-indel.genome.vcf`, `sv.vcf`, `cnv.vcf`), so kit or sample IDs in file names do not end up in catalogs or reports.
- **Retention**: finished analyses are deleted automatically after `GENOMESCOPE_JOB_TTL_HOURS` (default 7 days). **Delete** in History removes one immediately, and **Discard** removes a running one.
- **External lookups**: for annotation the backend sends variant identifiers and genomic coordinates to MyVariant.info, mygene.info and NCBI, including the positions of some of your own calls (structural / copy-number intervals, NCBI gene scans). It never sends files or raw VCF lines. Set `GENOMESCOPE_OFFLINE=1` to disable all external lookups; the catalog then contains only structural / copy-number rows and direct VCF reads.
- **Local-only API**: the API answers only on `127.0.0.1` / `localhost` (Host allow-list against DNS rebinding), allows CORS only from the local frontend, requires an `X-Requested-With: GenomeScope` header on every state-changing request, and sends `Cache-Control: no-store`.
- The ClinVar release (public data, ~1.2 GB) is cached once per user in `%LOCALAPPDATA%/genomescope/clinvar` (Windows) or `~/.cache/genomescope/clinvar` (configurable via `GENOMESCOPE_CLINVAR_DIR`).

## Configuration

Copy `.env.example` to `.env` (setup does this) and adjust as needed. `.env` is git-ignored; never commit credentials.

| Variable | Default | Purpose |
|---|---|---|
| `NCBI_API_KEY` | – | NCBI rate limit 3 → 10 req/s |
| `NCBI_EMAIL` | – | contact address sent to NCBI E-utilities |
| `NCBI_CACHE_TTL_HOURS` | `168` | NCBI response cache lifetime |
| `GENOMESCOPE_OFFLINE` | `0` | `1` = never contact external services |
| `GENOMESCOPE_JOBS_DIR` | `<temp>/genomescope_jobs` | where analyses run |
| `GENOMESCOPE_JOB_TTL_HOURS` | `168` | auto-delete finished analyses (`0` = keep) |
| `GENOMESCOPE_MAX_UPLOAD_MB` | `4096` | upload size limit per analysis |
| `GENOMESCOPE_MAX_CONCURRENT_JOBS` | `1` | analyses processed in parallel |
| `GENOMESCOPE_CLINVAR_DIR` | user cache dir | shared ClinVar release cache |
| `GENOMESCOPE_CORS_ORIGINS` / `GENOMESCOPE_ALLOWED_HOSTS` | local only | network exposure |

## Project layout

```text
backend/
  app/            FastAPI production layer: jobs API, explorer API, pipeline, temp-folder manager, worker
  genomescope/    core logic: VCF engine, catalog builder, ClinVar release scan, NCBI client
  tests/          offline test suite on a synthetic genome
  requirements.txt
frontend/         React 18 + Vite + Tailwind UI (New Analysis, History, Explorer)
docs/             ARCHITECTURE.md, NCBI_INTEGRATION.md, FRONTEND.md
setup.py / setup.bat / setup.sh     environment setup
run.bat / run.sh / launch.py        start backend + frontend (ports from port.json)
secrets.md        secret / PII audit trail of the migration
```

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the pipeline, job lifecycle, API, ClinVar scan and NCBI integration, [docs/NCBI_INTEGRATION.md](docs/NCBI_INTEGRATION.md) for the explorer / NCBI API reference, and [docs/FRONTEND.md](docs/FRONTEND.md) for the UI.

## Development

```bash
# backend tests (offline, synthetic data only)
cd backend
../.venv/Scripts/python -m pytest        # Windows
../.venv/bin/python -m pytest            # macOS / Linux

# live NCBI tests (network)
RUN_LIVE_NCBI=1 ../.venv/bin/python -m pytest tests/test_ncbi.py -k live

# frontend
cd frontend && npm run dev               # or: npm run build
```

Run the backend alone (from `backend/`): `python -m uvicorn app.main:app --host 127.0.0.1 --port 8000` with the `.venv` interpreter; interactive API docs are at `/docs`.

## License

Apache License 2.0 — see [LICENSE](LICENSE).
