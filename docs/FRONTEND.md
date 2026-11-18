# Frontend Guide

## Responsibilities
The ReactJS frontend (`frontend/`) provides a modernized clinical genomics dashboard and variant explorer built with React 18 (JavaScript), Vite, Tailwind CSS, and Lucide Icons. It talks only to the local backend (`/api`, proxied by Vite to `http://127.0.0.1:<backend port>`; the ports come from `GENOMESCOPE_BACKEND_PORT` / `GENOMESCOPE_FRONTEND_PORT`, set by `launch.py` from `port.json`, defaulting to 8000 / 5173; override the API base with `VITE_API_BASE`); all NCBI traffic goes through the backend.

## Source Layout (`frontend/src/`)
- `App.jsx` - app shell: top navigation (**New Analysis**, **History**, **Explorer**), the active analysis (remembered in `localStorage` as a per-viewer convenience; no genome data is stored in the browser), active-job counter
- `views/NewAnalysisView.jsx` - landing screen: three VCF slots, analysis name, ClinVar mode, privacy note, upload progress (XHR), then the live job monitor
- `views/HistoryView.jsx` - ledger of all analyses, auto-refresh every 4 s; Continue / Discard for queued or running jobs, Open / Download ZIP / Delete for completed ones
- `views/ExplorerView.jsx` - the explorer for one completed analysis (header incl. ClinVar / NCBI status pills, the three explorer tabs, catalog Explorer, Direct VCF Coordinate Streamer, inspector drawer); remounted per analysis
- `components/FileUpload.jsx` - per-slot dropzones (file picker, drag-and-drop) and `assignFiles` (multi-file drops and Ctrl+V paste are routed to the gVCF / SV / CNV slot by file name)
- `components/JobStatus.jsx` - job monitor: progress bar, five pipeline stages, error banner, actions; `StatusBadge`
- `components/LogViewer.jsx` - collapsible, live-polling tail of the pipeline log
- `api.js` - backend client. `setActiveJob(id)` makes every explorer / NCBI request carry `X-GenomeScope-Job: <id>`. Jobs: `fetchConfig`, `listJobs`, `getJob`, `getJobLogs`, `createJob` (multipart upload with progress, abortable), `discardJob`, `deleteJob`, `downloadZipUrl`. Catalog: `fetchSummary`, `fetchFilterOptions`, `fetchVariants`, `fetchVariantDetail`, `fetchDirectVcfReads`. NCBI: `fetchNcbiStatus`, `fetchNcbiVariant`, `fetchNcbiGene`, `searchNcbi`, `importNcbiVariants` (POST), `fetchVariantNcbi(id, { sync })` (GET comparison, or POST `/ncbi/sync` when `sync: true`). State-changing calls send `X-Requested-With: GenomeScope`; NCBI helpers surface FastAPI `detail` messages and "backend unreachable" errors.
- `components/NcbiStatusPill.jsx` - header pill
- `components/NcbiExplorer.jsx` - NCBI Live Explorer tab
- `components/NcbiVariantPanel.jsx` - drawer NCBI tab
- `components/ncbiShared.jsx` - shared NCBI UI: URL builders (`ncbiUrls`), `safeHref`, link components (`NcbiLink`, `RsidLink`, `RcvLink`, `RcvList`, `GeneLinks`), badges (`ClassificationBadge`, `ReviewStars`, `ZygosityBadge`), `GenotypeCell`, `ClinVarCard`, `FrequencyPanel`, `YourGenomeBox`, `GeneCard`, `PublicationsList`, `NcbiVariantDetail`, `ErrorsNotice`, `ErrorBanner`, `SkeletonCards`, `EmptyState`, and the `carriesNcbiAlt` / `isIndeterminateGenotype` helpers

## Screens
1. **New Analysis** - add the SNP & Indel gVCF (required) and optional SV / CNV VCFs by clicking, dropping, or pasting (Ctrl+V) files; `.vcf` and `.vcf.gz` are accepted and checked client-side (extension, empty file, size limit from `GET /api/config`). Choose the ClinVar mode (cached / download-refresh / skip) and start. The upload progress bar is followed by the job monitor (stages: Stage uploaded VCFs, Index gVCF, Build annotated catalog, ClinVar release verification, Package results ZIP) with a live log; **Discard** aborts and purges, **Open in Explorer** and **Download ZIP** appear when the analysis completes.
2. **History** - every analysis with name, original file names, size, created time, elapsed time and a status badge (animated while running). **Continue** reconnects the monitor to a running job.
3. **Explorer** - enabled once a completed analysis is opened; the features below.

## Key Explorer Features
1. **Multi-Field Search & Clinical Filter Matrix**:
   - Instant global search bar plus dedicated filter inputs/controls for all 11 required fields:
     - `RCV ID`
     - `Chromosome (chr)`
     - `Position` (exact or range)
     - `Variant ID` (rsID / Manta / Canvas)
     - `Gene`
     - `Risk Frequency`
     - `Variant Type`
     - `Your Status`
     - `Condition`
     - `Classification`
     - `Evidence`
   - Quick-filter presets (e.g., "Carried Clinical Findings", "Pathogenic / Likely Pathogenic", "Risk Factors", "Pharmacogenomics / Drug Response", "Copy Number & Structural Variants", "All VCF Reads").
   - Identifier searches absent from the catalog may be resolved live by the backend (MyVariant.info, then NCBI ClinVar import), so new rows can appear for an rsID / RCV / gene / locus that was not indexed.
2. **Results Table & Card Views**:
   - Displays Variant ID, RCV ID, Gene, Chromosome:Position, Your Data (genotype alleles + read depth `DP` / quality `GQ` / allelic depth `AD`), Risk Version (`ALT` vs `REF`), Risk Frequency, Variant Type, Your Status badge, Condition, Classification badge, and Evidence rating.
   - Outbound NCBI links: rsIDs link to dbSNP, RCV accessions to ClinVar, gene symbols to NCBI Gene (non-rsID / synthetic CNV-SV ids and non-gene tokens such as `Intergenic` render as plain text).
3. **Interactive Variant & Raw VCF Read Inspector Drawer**:
   - Clicking any row opens a comprehensive clinical & sequencing inspector with two tabs:
   - **VCF Read & Catalog**:
     - Allele read balance bar (`AD` ref vs alt read counts, `DP`, `GQ`, `QUAL`, `FILTER`)
     - Clinical significance, RCV accession links (ClinVar), condition details, and evidence stars; header rsID / RCV / gene link to dbSNP / ClinVar / NCBI Gene
     - Exact raw VCF line from the analysis' VCFs and neighboring genomic reads at that locus
   - **NCBI (ClinVar · dbSNP · Gene · PubMed)** (`NcbiVariantPanel`, `GET /api/variants/{id}/ncbi`):
     - Catalog vs NCBI comparison table (Allele (ClinVar ALT), Classification, Review status, Condition, RCV ID, Risk allele frequency) with match / mismatch / n/a per field
     - **Reload** and **Sync catalog from NCBI** (`POST /api/variants/{id}/ncbi/sync`; highlighted amber when any field mismatches, disabled when no ClinVar record resolved). The result notice shows success or the backend's `Sync skipped: ...` reason; after a successful sync the summary, the result list and the open drawer record are refreshed.
     - Full `NcbiVariantDetail` view (see below). Stale responses are discarded when the selected record changes.
4. **Direct VCF Coordinate Streamer** tab to inspect raw VCF reads at any `chr:position` directly from the analysis' VCFs.
5. **NCBI Live Explorer** tab (`NcbiExplorer`, kept mounted so results survive tab switches), three modes:
   - **ClinVar search** (`GET /api/ncbi/search`): free term (RCV, rsID, VCV, gene, condition, text), result limit, example chips; shows the translated E-utilities query sent to ClinVar.
   - **Gene panel scan** (`GET /api/ncbi/gene/{symbol}`): gene symbol, significance (Pathogenic + Likely pathogenic / Risk factor / Association / Drug response / All classifications), variant limit; NCBI Gene card and stat tiles (ClinVar total, Variants scanned, You carry); "Carried only" toggle; carried rows are highlighted and listed first.
   - **Variant lookup** (`GET /api/ncbi/variant`): rsID, RCV, VCV or `chr:position` parsed client-side (`parseLookupInput`), optional publications; renders `NcbiVariantDetail`.
   - Search and gene results share a table (variant / alleles, rsID, classification, review stars + last evaluated, conditions, your genotype) with row checkboxes and **Import selected into catalog** (max 50 per request, `POST /api/ncbi/import`). The import notice lists imported count, skipped ids with reasons and upstream notices, and offers **Open in Explorer** to jump to the catalog Explorer filtered by rsID / RCV / gene; imported rows get **View in Explorer**, and rows with an identifier or locus have **Details** (opens it in Variant lookup). Older responses arriving after a newer request are ignored.
6. **NCBI detail view** (`NcbiVariantDetail`, used by Variant lookup and the drawer NCBI tab):
   - Resolved GRCh38 locus + rsID, and outbound links to ClinVar / dbSNP / Gene / PubMed
   - Upstream `errors` shown as an "NCBI upstream notices" box
   - **Your genome vs NCBI** (`YourGenomeBox`): verdict sentence, zygosity badge, NCBI REF -> ALT, your genotype, status, carries NCBI ALT (Yes / No / Unknown), read metrics, FILTER, and the raw VCF line
   - ClinVar cards: classification badge (raw + normalized), review stars, last evaluated, locus, protein change, genes, conditions with MedGen and OMIM links
   - dbSNP frequency panel (placement, alleles, class, per-study frequencies with the preferred study and the ClinVar-ALT risk allele frequency), NCBI Gene card (cytoband, GRCh38 span, OMIM MIM links, aliases, summary), PubMed publications list
7. **NCBI status pill** (`NcbiStatusPill`, header): calls `GET /api/ncbi/status` once on load; shows `NCBI live` (green), `NCBI offline` (amber; cached data only) or `NCBI ?` (status call failed). Tooltip: rate limit, API key configured or not, cache entries and TTL, notes, check time.

## Genotype Display Semantics
- `carriesNcbiAlt(yourGenotype)` is true only when the backend reports `carries_ncbi_alt === true` (fallback: a carried zygosity with `FILTER = PASS` when the flag is missing), matching the backend's `carried_count`. Highlights, "Carried only" and "You carry" use it.
- `ZygosityBadge` labels: Homozygous ALT, Heterozygous, Hemizygous ALT, Homoplasmic ALT, ALT absent, Reference allele, FAIL – unconfirmed (`low_quality`), No call, Indeterminate (`unknown`, e.g. haplotype / diplotype records that a single-locus read cannot confirm).

## Link Safety
- Backend-provided URLs are used as `href` only when they are `http(s)` (`safeHref`); otherwise text is rendered without a link. External links open in a new tab.
- `ncbiUrls` builds ClinVar (RCV / variation), dbSNP, NCBI Gene (id or `symbol[sym] AND human[orgn]` search), MedGen, OMIM and PubMed URLs.
