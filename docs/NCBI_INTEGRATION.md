# Backend Core, NCBI Integration & Explorer API Reference

## Responsibilities
The core package (`backend/genomescope/`), served by the FastAPI app in `backend/app/`, is responsible for:
1. **VCF File Indexing & Fast Binary-Search Reading** (`vcf_engine.py`):
   - Reads the VCFs of the active analysis (`<job>/inputs/`): `snp-indel.genome.vcf` (required), `sv.vcf`, `cnv.vcf` (optional). Uploaded files are stored under these canonical names.
   - Uses a 256 KB sparse byte-offset index (plus an exact entry at the first record of every contig, so small contigs such as chrM / chrY are never missed; `chr`-prefixed and bare contig names are equivalent) over `snp-indel.genome.vcf` to resolve any `(chr, position)` coordinate in `< 2 ms` via binary search, handling both explicit variant calls (`ALT != '.'`) and gVCF homozygous-reference blocks (`END=...`).
   - Parses all Structural Variants (`sv.vcf`) and Copy Number Variants (`cnv.vcf`) and maps their genomic intervals against GRCh38 gene coordinates.
2. **Clinical Annotation & SQLite Catalog** (`indexer.py`):
   - Stores pre-indexed ClinVar + gnomAD + dbSNP annotations matched against the actual VCF reads in the analysis' `work/genetic_variants.db` (table `variants`).
   - Performs live on-demand lookup against the raw VCF files and MyVariant.info (GRCh38) when a user searches for a new `rcv_id`, `variant_id` (`rsID`), `gene`, or `chr` + `position` coordinate not yet cached. If MyVariant adds nothing for an identifier search, `GET /api/variants` falls back to a live NCBI ClinVar import (`ncbi_service.fallback_import_for_search`).
   - Migrates the schema in place: `ensure_ncbi_columns()` (called from `init_db()`) adds `ncbi_variation_id`, `ncbi_last_evaluated`, `ncbi_synced_at`, `frequency_source` via guarded `ALTER TABLE` plus index `idx_variants_ncbi_vid`.
   - Upserts with `INSERT ... ON CONFLICT(unique_key) DO UPDATE` (row ids stay stable). Rows with `ncbi_synced_at` set keep their NCBI-synced fields (`all_rcv_ids`, `condition`, `classification`, `evidence`, `evidence_stars`, `review_status`, `risk_frequency`, `risk_frequency_text`) when a non-NCBI upsert arrives; NCBI provenance columns are only overwritten by non-NULL NCBI values.
3. **ClinVar Release Scan** (`alleles.py`, `clinvar_release.py`, `clinvar_scan.py`): complete, allele-specific scan of the gVCF against a locally indexed, MD5-verified ClinVar `variant_summary.txt.gz`; reconciles the catalog (adds missing carriers, re-adjudicates every existing snp-indel row, backs up first, checks post-conditions). Run as the `clinvar` pipeline stage, or by CLI: `python -m genomescope.clinvar_scan --workspace <job_dir> [--update] [--dry-run] [--include-vus]`. Full description, columns (`clinvar_variation_id`, `clinvar_release`, `clinvar_match`), limitations and outputs are in `docs/ARCHITECTURE.md`, section "ClinVar Release Scan".
4. **Live NCBI Integration** (`ncbi_client.py`, `ncbi_service.py`):
   - `ncbi_client.py` - three layers:
     1. Pure, I/O-free parsers / normalisers (unit-tested against fixtures): `parse_esearch`, `parse_esummary_result`, `parse_clinvar_summary`, `parse_dbsnp_summary`, `parse_gene_summary`, `parse_pubmed_summary`, `parse_elink_ids`, `parse_global_maf`, `parse_refsnp_frequencies`, `parse_spdi`, `spdi_to_vcf` (SPDI -> left-aligned VCF, indels need an anchor base), `alleles_equivalent`, `risk_allele_frequency_from_dbsnp` / `_from_refsnp`, `normalize_clinvar_classification`, `review_status_to_stars`, `translate_search_term`, `sanitize_phrase`, `is_gene_symbol`.
     2. Infrastructure: global thread-safe `RateLimiter` and persistent `NcbiCache` (table `ncbi_cache(key, payload, fetched_at)` in the active analysis' catalog DB; without an analysis the cache is simply bypassed).
     3. `NcbiClient` (httpx): `esearch`, `esummary` (chunks of 100 ids), `elink`, `fetch_sequence` (nuccore efetch FASTA), `refsnp` (Variation Services), `ping`. All traffic goes through `NcbiClient._transport`, the single monkeypatch target for tests. `get_client()` / `set_client()` manage the process-wide instance.
   - Base URLs: E-utilities `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`, Variation Services `https://api.ncbi.nlm.nih.gov/variation/v0/`. E-utilities requests carry `tool=genomescope_wgs_explorer`, plus `email` / `api_key` when configured.
   - `ncbi_service.py` - builds the `/api/ncbi/*` response objects: resolves ClinVar variation ids, fills VCF alleles for indels (anchor base via clustered nuccore efetch, max 6 fetches per request), fetches dbSNP / Gene / PubMed, reads the user's genotype at each NCBI locus (`build_your_genotype`), computes the ClinVar-ALT risk allele frequency, and imports / syncs records into the catalog. All NCBI failures are non-fatal and collected into each response's `errors` list.
5. **Explorer API Endpoints** (`backend/app/api/explorer.py`; all except the two status routes need the `X-GenomeScope-Job` header, see `docs/ARCHITECTURE.md`):
   - `GET /api/health` - Health check; with a selected analysis also catalog size and VCF file status; includes `ncbi: {api_key_configured, cache_entries}`
   - `GET /api/summary` - Aggregate statistics across VCF files, classifications, zygosity statuses, variant types, and top genes
   - `GET /api/variants` - Multi-parameter search & filter across `q`, `rcv_id`, `chr`, `position`, `variant_id`, `gene`, `risk_frequency`, `variant_type`, `your_status`, `condition`, `classification`, `evidence`, `source_file`, `carried_only`, `sort_by`, `sort_dir`, `page`, `page_size`. Zero results for an identifier absent from the catalog trigger MyVariant, then NCBI import (see below)
   - `GET /api/variants/{record_id}` - Detailed single variant read + surrounding VCF context reads
   - `GET /api/vcf/read` - Direct real-time VCF reader for any `chr` and `position` (or range) across the analysis' VCF files
   - `GET /api/filters` - Distinct filter option values and counts for UI dropdowns
   - `GET /api/clinvar/status` - State of the local ClinVar verification: `{state: verified|stale|never_run|index_missing|no_analysis, message, index: {release_date, ...}, last_scan: {coverage, clinvar_matches_by_tier, catalog, ...}, catalog: {rows_by_clinvar_match, carried_rows_without_verified_allele}}`
   - `GET /api/ncbi/status` - NCBI reachability and client configuration
   - `GET /api/ncbi/variant` - Live single-variant lookup across ClinVar / dbSNP / Gene / PubMed + the user's genotype
   - `GET /api/ncbi/gene/{symbol}` - Live ClinVar gene panel scan against the user's genome
   - `GET /api/ncbi/search` - Live ClinVar search with search-term translation
   - `POST /api/ncbi/import` - Import ClinVar variation ids into the catalog
   - `GET /api/variants/{record_id}/ncbi` - Read-only catalog-vs-NCBI comparison for a catalog row
   - `POST /api/variants/{record_id}/ncbi/sync` - Overwrite a catalog row's clinical fields with NCBI values

## NCBI Endpoint Reference
All NCBI bodies include `errors: [str]` (non-fatal upstream notices) and `fetched_at` (UTC ISO). Invalid parameters return `400`; missing catalog rows return `404`.

### `GET /api/ncbi/status`
Runs an uncached single-attempt ClinVar `esearch` probe (8 s timeout).
Response: `{reachable, api_key_configured, rate_limit_per_sec, cache_entries, cache_ttl_hours, services: {eutils, variation}, errors, fetched_at}`.

### `GET /api/ncbi/variant`
| Param | Notes |
|---|---|
| `rsid` | `rs671` or `671` (prefixed with `rs` automatically) |
| `rcv_id` | `RCV000019449` (optional `.version`) |
| `vcv_id` | `VCV000017865` or a bare variation id |
| `chr` + `position` | GRCh38, `chr` in `1-22, X, Y, MT` (optional `chr` prefix), integer position (commas / `chr:` prefix tolerated). Used only when no identifier is given |
| `include_publications` | default `true` |

At least one of `rsid`, `rcv_id`, `vcv_id` or `chr` + `position` is required. ClinVar records are filtered to those matching the supplied identifiers; for an rsID up to 100 (hard cap 500) records are checked, with a notice when ClinVar lists more. With no catalog-specific preference, the primary record is the simple (non-haplotype) record whose allele the user carries, else the most severe / best-reviewed one.

Response: `{query, resolved: {chr, position, rsid}, clinvar: [NcbiClinVarRecord, primary first], dbsnp, gene, publications (<= 10), your_genotype, risk_allele_frequency: {allele, frequency, study} | null, links: {clinvar, dbsnp, gene, pubmed}, errors, fetched_at}`.
- `NcbiClinVarRecord`: `variation_id, vcv_id, title, variant_type, rcv_accessions, classification, classification_normalized, somatic_classification, review_status, review_stars (0-4), last_evaluated, conditions [{name, medgen_id, omim_ids, url (MedGen)}], genes [{symbol, gene_id}], chr, start, stop, assembly ("GRCh38"), ref, alt, rsid, molecular_consequences, protein_change, url`.
- `dbsnp`: `rsid, chr, position, ref, alts, spdi, variant_class, function_class, clinical_significance, genes, frequencies [{study, allele, frequency, count}], preferred_frequency, url`.
- `gene`: `gene_id, symbol, name, summary, aliases, chr, start, stop, map_location, mim_ids, url`.
- `publications`: `[{pmid, title, journal, year, authors, url}]`.
- `your_genotype`: see "Genotype Matching" below (`null` when no GRCh38 locus resolves).

### `GET /api/ncbi/gene/{symbol}`
| Param | Notes |
|---|---|
| `symbol` | path; `^[A-Za-z0-9][A-Za-z0-9\-.]{0,29}$` |
| `significance` | `pathogenic` (default, pathogenic + likely pathogenic) / `risk` (risk factor, association, protective, affects) / `drug` (drug response) / `all` |
| `limit` | 1-100, default 40 |

Fetches the human NCBI Gene record, then a ClinVar pool (`max(limit, min(200, max(100, 2*limit)))` ids, supplemented with 2-4 star records when the gene has more hits). When ClinVar has more records than the pool, it also searches ClinVar at the user's own carried positions inside the gene (`user_carried_positions`: up to 180 positions, gene span <= 3 Mb, 2 kb padding, indels also contribute POS+1; queried in chunks of 60, up to 100 ids). Variants are genotyped and ranked carriers first.
Response: `{gene, significance, clinvar_total, variants: [{clinvar, your_genotype}], carried_count, errors, fetched_at}`.

### `GET /api/ncbi/search`
| Param | Notes |
|---|---|
| `term` | 1-200 chars: RCV, VCV, rsID, `chr:pos`, gene symbol, condition, or free text |
| `limit` | 1-100, default 25 |

`translate_search_term` maps the term to ordered ClinVar queries (`RCV...`, `VCV...[VACC]`, `rs...[VRID]`, `<chr>[chr] AND <pos>[chrpos38]`, `"<SYMBOL>"[gene]`, `"<text>"[dis]`, fallback `"word"[All Fields] AND ...`); the first query with hits wins.
Response: `{term, translated_query, total, items: [{clinvar, your_genotype}], errors, fetched_at}`.

### `POST /api/ncbi/import`
Requires header `X-Requested-With: GenomeScope`. Body `{"variation_ids": ["17865", "VCV000017865", ...]}` with 1-50 ids.
Response: `{imported, record_ids, skipped: [{variation_id, reason}], errors, fetched_at}`. Skip reasons: invalid ClinVar variation id, not found in NCBI ClinVar / NCBI ClinVar unavailable, no GRCh38 location in ClinVar, no VCF coverage at `chr:pos`.
Behaviour: bypasses the response cache (fresh data), 90 s budget, at most 10 Variation Services frequency fallbacks per call. Updates an existing row in place when it has the same `ncbi_variation_id` or is an unlinked row at the same chr / position / source file / variant id / ALT sharing an RCV with the record (haplotype / diplotype records never take over such rows); otherwise inserts a row with `unique_key = NCBI:{chr}:{pos}:VCV{9-digit id}`. `is_carried = 1` only when the PASS call carries exactly the ClinVar ALT (`carries_ncbi_alt is True`).

### `GET /api/variants/{record_id}/ncbi` and `POST /api/variants/{record_id}/ncbi/sync`
- `GET` is read-only. `GET ...?sync=true` returns `405` (sync moved to `POST`, so GETs stay side-effect free).
- `POST .../sync` requires `X-Requested-With: GenomeScope` and bypasses the response cache.
- The lookup uses the row's rsID when `variant_id` is one; otherwise its RCV and / or stored `ncbi_variation_id`; otherwise `chr` + `position`. It prefers the ClinVar record matching the row's RCVs (`rcv_id` + `all_rcv_ids`), stored `ncbi_variation_id`, or ALT (`risk_version`). A notice is added when the stored variation id is no longer linked to the rsID.
- Structural rows (non-SNP source file, symbolic ALT, `RCV-CNV-` / `RCV-SV-` / `RCV-LOH-` ids, CNV / SV / LOH types) are not looked up; the body explains why.
- Response: the `/api/ncbi/variant` body plus `{record_id, catalog_vs_ncbi: [{field, catalog, ncbi, match}], synced}` with fields `allele`, `classification`, `review_status`, `condition`, `rcv_id`, `risk_frequency`; `match` is `null` when NCBI has nothing to compare.
- Sync refusals (appended to `errors`, prefixed `Sync skipped:`): no ClinVar record resolved; the record is a haplotype / diplotype; no record matches the row's RCV / variation id / allele (the shown record describes another allele); the user's VCF call is a different ALT; structural / copy-number row.
- Fields written on success: `classification`, `review_status`, `evidence` (`... | NCBI ClinVar VCV...`), `evidence_stars`, `condition` (existing catalog text kept when it already names the NCBI condition, or NCBI only has a generic `-related genetic phenotype`), `all_rcv_ids`, `ncbi_variation_id`, `ncbi_last_evaluated`, `ncbi_synced_at`, and `risk_frequency` / `risk_frequency_text` (`<freq> [dbSNP <study>]`) / `frequency_source` (`NCBI dbSNP <study>`) when a frequency resolves.

### `GET /api/variants` NCBI fallback
Runs only when the search returns 0 rows, the queried identifier (RCV, rsID, VCV, gene symbol, `chr:position`; also parsed from `q`) is absent from the catalog regardless of other filters, and the MyVariant resolver added nothing. Single attempt per request, 20 s budget, up to 60 ClinVar ids searched, up to 20 located records imported (simple records first, max 3 Variation Services fallbacks). Free text and condition searches never use it.

## Genotype Matching (`build_your_genotype`)
Reads `snp-indel.genome.vcf` at the ClinVar / dbSNP locus, collecting every row at the position plus a flank (gVCFs split one site into several rows; a 1-bp deletion's VCF anchor sits one base earlier), and compares the NCBI SPDI to each overlapping row's ALT list with `alleles_equivalent` (SPDI vs left-aligned VCF indel representations).

`your_genotype`: `{found, chr, position, ref, alt, genotype, alleles_called, your_status, is_carried, dp, gq, ad, filter, is_reference_block, raw_vcf_line, source_file, carries_ncbi_alt, zygosity_for_ncbi_alt}`. For hom-ref blocks `ref` is the reference allele at the queried locus, not the block's first base.

| Situation | `zygosity_for_ncbi_alt` | `carries_ncbi_alt` |
|---|---|---|
| No call (`./.`) | `no_call` | `null` |
| ClinVar `=` (reference-allele) record, user has a REF copy / no REF copy | `reference_allele` / `absent` | `null` / `false` |
| Hom-ref block | `absent` | `false` |
| NCBI ALT matched: all copies / some copies / haploid call | `homozygous` / `heterozygous` / `hemizygous` (`homoplasmic` on MT) | `true` |
| Matched row lacks the NCBI ALT | `absent` (status `Different ALT Allele (Not the NCBI Variant)` if another ALT is carried) | `false` |
| dbSNP-only multi-allelic locus (no single target allele) | from GT | `null` |
| Representation undecidable | `unknown` | `null` |
| Matched row `FILTER != PASS` (not a no-call) | `low_quality` | `null` |
| Haplotype / diplotype / compound record (`obj_type` or >1 `variation_set`) that would otherwise be `true` | `unknown` | `null` |

Only `carries_ncbi_alt is True` counts as carried (ranking, `carried_count`, imported `is_carried`, sync / comparison allele match).

## NCBI Client Behaviour & Environment
| Variable | Default | Effect |
|---|---|---|
| `NCBI_API_KEY` | unset | 10 req/s instead of 3; sent as `api_key` on E-utilities |
| `NCBI_EMAIL` | unset | Sent as `email` on E-utilities |
| `NCBI_CACHE_TTL_HOURS` | `168` | Response cache TTL |
| `GENOMESCOPE_CORS_ORIGINS` | local Vite dev/preview origins | CORS allow-list (comma list) |
| `GENOMESCOPE_ALLOWED_HOSTS` | `127.0.0.1, localhost, ::1, [::1], testserver` | Host-header allow-list (`*` = any) |
| `GENOMESCOPE_OFFLINE` | unset | `1` = NCBI (and MyVariant / mygene / ClinVar FTP) are never contacted; NCBI calls fail fast with a notice |

- Rate limit: global minimum-interval limiter (safety factor 1.25 without key, ~0.42 s spacing; 1.1 with key). `X-RateLimit-Remaining: 0` or HTTP 429 -> 1 s cooldown for all threads.
- Retries: max 3 attempts, backoff 0.5 s x 2^n honouring `Retry-After` (<= 10 s); 4xx other than 429 are not retried. Timeouts 15 s read / 3 s connect.
- Circuit breaker: 2 consecutive transport failures -> fail fast for 60 s, serving stale cache when present; the status probe bypasses it.
- `NcbiClient.budget(seconds, attempts, fresh)` bounds every request in a block: lookup / search / catalog comparison 45 s, gene scan 60 s, import 90 s, fallback 20 s with 1 attempt; `fresh=True` (import, sync) skips cached responses but refreshes the cache.
- Cache key: URL + sorted params excluding `api_key`, `email`, `tool`. Empty esearch results expire after 1 h. Expired entries are served when NCBI is unreachable.
- Security: CORS limited to the local frontend, no credentials; `TrustedHostMiddleware` Host allow-list; `Cache-Control: no-store` on every API response; state-changing routes (`POST /api/jobs`, `POST /api/jobs/{id}/discard`, `DELETE /api/jobs/{id}`, `POST /api/ncbi/import`, `POST /api/variants/{id}/ncbi/sync`) require `X-Requested-With: GenomeScope`.

## Testing
- `python -m pytest` (from `backend/`) - offline. `test_ncbi.py` parser tests use trimmed public NCBI responses in `tests/fixtures/ncbi/`; endpoint tests replace `NcbiClient._transport` with a fixture router (`FakeTransport`) and run against a fresh copy of the synthetic analysis (`tests/synthetic.py`) with an empty `ncbi_cache` (fixture `temp_db`).
- `test_clinvar_scan.py` - offline; builds a synthetic ClinVar summary + gVCF and a throwaway analysis workspace (fixtures `ws`, `catalog`).
- Live tests (`test_live_rs671`, `test_live_rcv_and_gene`) query real NCBI and run only with `RUN_LIVE_NCBI=1`, e.g. `RUN_LIVE_NCBI=1 python -m pytest tests/test_ncbi.py -k live` (PowerShell: `$env:RUN_LIVE_NCBI='1'; python -m pytest tests/test_ncbi.py -k live`).
