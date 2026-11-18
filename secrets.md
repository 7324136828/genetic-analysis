# Secrets & Personal-Data Audit

This file records potential secrets and personal data (PII / genetic data) found while
productionizing the original prototype into this repository, and what was done about them.

**No secret values, identifiers, genotypes or findings are stored in this report.**

- Audit date: 2026-11-17
- Scope: the original prototype (sources, docs, tests, build output, data folders) and every file of this repository (tracked and untracked, excluding `node_modules/`, `.venv/`, `frontend/dist/`)
- Tools: `detect-secrets` 1.5.0 (all default plugins), plus pattern checks (credential keywords, private-key headers, AWS / GitHub / OpenAI-style tokens, credentials in URLs, e-mail addresses, data file extensions, and the prototype's sample / kit / vendor / genotype / finding strings)

## 1. Credentials

| File | Line | Secret Type | Action | Status |
|---|---:|---|---|---|
| `backend/tests/test_ncbi.py` | 393 | "Secret Keyword" (test placeholder value for `NCBI_API_KEY`, not a credential) | None needed; placeholder kept | False positive |
| `docs/ARCHITECTURE.md` | 209 | "Secret Keyword" (prose: cache keys exclude `api_key`) | None needed | False positive |

- No hard-coded credentials, tokens, private keys, connection strings or `.env` files were found in the prototype or in this repository.
- Credentials are runtime configuration only: `NCBI_API_KEY` / `NCBI_EMAIL` are read from the environment (or the git-ignored `.env`); `.env.example` contains empty placeholders. The API key is never cached, logged or returned (`/api/ncbi/status` reports only whether one is configured).
- No credential rotation is required.
- Git history: this repository's only prior commit contains `README.md`, `LICENSE` and `.gitignore`; no secret or personal data was ever committed.

## 2. Personal and genetic data (PII) — not migrated

The prototype was built around one person's genome. None of the following was copied into this repository:

| Prototype location | Data category | Action | Status |
|---|---|---|---|
| `input/genome/` | Raw whole-genome VCFs (SNP/indel gVCF, SV, CNV) of one individual; file names embed a sample/kit identifier, vendor and date | Not migrated | Excluded |
| `backend/data/genetic_variants.db` | Variant catalog derived from that genome (genotypes, carried findings) | Not migrated | Excluded |
| `backend/data/backups/` | Backup copy of that catalog | Not migrated | Excluded |
| `backend/data/clinvar/scan_report.json`, `scan_carried.tsv` | ClinVar scan results listing the individual's carried variants | Not migrated | Excluded |
| `backend/data/vcf_sparse_index.json` | Byte-offset index of the individual's gVCF | Not migrated | Excluded |
| `backend/data/clinvar/variant_summary.txt.gz`, `clinvar_grch38.db` | Public ClinVar release (not personal; ~1.2 GB) | Not migrated; downloaded on demand into a per-user cache | Excluded |
| `frontend/dist/` | Build output embedding the sample identifier | Not migrated (rebuilt from sources) | Excluded |
| `.pytest_cache/`, `__pycache__/`, `node_modules/` | Caches / dependencies | Not migrated | Excluded |

## 3. Personal data embedded in source code, tests and docs — remediated

| Prototype location (migrated to) | Category | Action | Status |
|---|---|---|---|
| `backend/vcf_engine.py` (`backend/genomescope/vcf_engine.py`) | Sample/kit identifier, vendor and date in the hard-coded VCF file names; sample id in file metadata | Replaced with canonical names `snp-indel.genome.vcf` / `sv.vcf` / `cnv.vcf`; per-analysis paths | Remediated |
| `backend/main.py` (`backend/app/services/catalog.py`, `backend/app/main.py`) | Sample identifier, coverage and sequencing vendor in API responses; identifier in an SQL literal | Sample label now comes from the user's analysis name; parameterized query on the canonical name | Remediated |
| `backend/indexer.py` (`backend/genomescope/indexer.py`) | Sequencing vendor in a review-status string | Replaced with a neutral label | Remediated |
| `frontend/index.html`, `frontend/src/App.jsx` (`frontend/src/views/ExplorerView.jsx`) | Sample identifier in page title, header badge, KPI text, CSV file name; input file size | Replaced with the analysis label / neutral text | Remediated |
| `frontend/src/App.jsx` | Hard-coded list of the individual's carried genes and clinical findings (pharmacogenomic responses, risk factors, a CNV) | Removed; chips are now generated from the selected analysis' summary | Remediated |
| `frontend/src/components/ncbiShared.jsx` | Path of the personal input folder | Replaced | Remediated |
| `backend/clinvar_scan.py` docstring, root `agent.md` | Reference to a specific pathogenic variant / gene used as the motivating example (indicative of a personal finding) | Rewritten generically | Remediated |
| `backend/tests/test_api.py`, `test_ncbi.py` | Assertions on the individual's real genotypes, read depths and allele counts at named loci; sample identifier | Rewritten against a synthetic genome (`backend/tests/synthetic.py`) with deliberately different, fictitious genotypes | Remediated |
| `backend/tests/test_clinvar_scan.py` | Sample identifier as VCF sample column; a gVCF line, gene, variant id and condition mirroring a personal finding | Replaced with fictitious gene / variant / condition / metrics | Remediated |
| `agent.md`, `backend/agent.md`, `frontend/agent.md` (`docs/*.md`) | Sample identifier, vendor, input file names / sizes / record counts, the personal-finding example | Rewritten; personal details removed | Remediated |

## 4. Runtime protection (so new personal data is not leaked)

- Uploaded genomes are processed only in `<system temp>/genomescope_jobs/<uuid>/` (owner-only permissions on POSIX), never inside the repository; folders are purged on Discard / Delete and automatically after `GENOMESCOPE_JOB_TTL_HOURS`.
- `.gitignore` excludes `.env*` (except `.env.example`), VCF / BAM / CRAM / FASTQ files, SQLite databases, ClinVar files, scan reports, sparse indexes, result ZIPs and `data/` / `input/` folders.
- API responses carry `Cache-Control: no-store`; the API is local-only (Host allow-list, restricted CORS, CSRF header on state-changing routes); internal temp paths are never returned.
- TLS certificate verification, disabled in the prototype for MyVariant.info / mygene.info, is enabled again.
- `GENOMESCOPE_OFFLINE=1` disables every external lookup.

## 5. Verification

- Final pattern scan over all repository files: no sample/kit identifier, vendor/date file-name fragments, personal findings, real genotype metrics, user paths, user e-mail or account names; the only e-mail address is the reserved `example.org` test address.
- `detect-secrets`: only the two false positives listed in section 1.
- No data files (`.vcf*`, `.db`, `.tsv`, scan reports, indexes) are present in the repository.
- `secrets.md` (this file) contains no secret values or personal data.
