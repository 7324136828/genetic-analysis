import React, { useState } from 'react';
import {
  AlertTriangle,
  BookOpen,
  ChevronDown,
  ChevronUp,
  Dna,
  ExternalLink,
  FlaskConical,
  Library,
  Star,
  User,
} from 'lucide-react';

/*
 * Shared presentational building blocks for NCBI (ClinVar / dbSNP / Gene / PubMed) data.
 * Every field coming from NCBI is treated as optional: NCBI records are sparse.
 * All NCBI text is rendered as React text nodes (never as HTML).
 */

// ---------------------------------------------------------------------------
// Identifier helpers & URL builders
// ---------------------------------------------------------------------------

const RCV_RE = /^RCV\d+(\.\d+)?$/i;
const VCV_RE = /^VCV\d+(\.\d+)?$/i;
const RSID_RE = /^rs\d+$/i;
// Same rule as the backend (_GENE_RE in backend/app/api/explorer.py). Rejects synthetic catalog
// labels such as CHR1_CNV / CHR10_LOCUS ('_'), so they render as plain text.
export const GENE_SYMBOL_RE = /^[A-Za-z0-9][A-Za-z0-9\-.]{0,29}$/;
const NON_GENE_TOKENS = new Set(['', '-', '.', 'NA', 'N/A', 'NONE', 'INTERGENIC', 'UNKNOWN', 'MULTIPLE']);

export const isRealRcv = (id) => typeof id === 'string' && RCV_RE.test(id.trim());
export const isRealVcv = (id) => typeof id === 'string' && VCV_RE.test(id.trim());
export const isRsid = (id) => typeof id === 'string' && RSID_RE.test(id.trim());

export const ncbiUrls = {
  clinvarRcv: (rcv) => `https://www.ncbi.nlm.nih.gov/clinvar/${encodeURIComponent(rcv.trim().toUpperCase())}/`,
  clinvarVariation: (vid) => `https://www.ncbi.nlm.nih.gov/clinvar/variation/${encodeURIComponent(String(vid))}/`,
  dbsnp: (rsid) => `https://www.ncbi.nlm.nih.gov/snp/${encodeURIComponent(rsid.trim().toLowerCase())}`,
  geneId: (geneId) => `https://www.ncbi.nlm.nih.gov/gene/${encodeURIComponent(String(geneId))}`,
  geneSymbol: (symbol) =>
    `https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(symbol.trim())}%5Bsym%5D+AND+human%5Borgn%5D`,
  medgen: (id) => `https://www.ncbi.nlm.nih.gov/medgen/${encodeURIComponent(String(id))}`,
  omim: (id) => `https://www.omim.org/entry/${encodeURIComponent(String(id))}`,
  pubmed: (pmid) => `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(String(pmid))}/`,
};

/** Only allow http(s) URLs coming from the backend / NCBI to be used as hrefs. */
export function safeHref(url) {
  if (typeof url !== 'string') return null;
  const u = url.trim();
  return /^https?:\/\//i.test(u) ? u : null;
}

export const asArray = (v) => (Array.isArray(v) ? v.filter((x) => x !== null && x !== undefined) : []);

export function splitGeneSymbols(geneField) {
  if (!geneField || typeof geneField !== 'string') return [];
  return geneField
    .split(/[,;/|\s]+/)
    .map((g) => g.trim())
    .filter((g) => g && GENE_SYMBOL_RE.test(g) && !NON_GENE_TOKENS.has(g.toUpperCase()));
}

export function splitRcvIds(field) {
  if (!field || typeof field !== 'string') return [];
  return Array.from(
    new Set(
      field
        .split(/[,;|\s]+/)
        .map((t) => t.trim())
        .filter(Boolean)
    )
  );
}

export function formatNumber(n) {
  const num = Number(n);
  return n === null || n === undefined || n === '' || Number.isNaN(num) ? null : num.toLocaleString();
}

export function formatFrequency(freq) {
  const f = Number(freq);
  if (freq === null || freq === undefined || Number.isNaN(f)) return '—';
  if (f === 0) return '0';
  if (f < 0.0001) return f.toExponential(2);
  if (f < 0.01) return `${(f * 100).toFixed(3)}%`;
  return `${(f * 100).toFixed(2)}%`;
}

export function formatLocus(chr, start, stop) {
  if (!chr && (start === null || start === undefined)) return null;
  const s = formatNumber(start);
  const e = formatNumber(stop);
  const chrLabel = chr ? `chr${String(chr).replace(/^chr/i, '')}` : 'chr?';
  if (!s) return chrLabel;
  return e && Number(stop) > Number(start) ? `${chrLabel}:${s}-${e}` : `${chrLabel}:${s}`;
}

export function formatFetchedAt(iso) {
  if (!iso) return null;
  const d = new Date(iso);
  return Number.isNaN(d.getTime()) ? String(iso) : d.toLocaleString();
}

const CARRIED_ZYGOSITIES = new Set(['homozygous', 'heterozygous', 'hemizygous', 'homoplasmic']);

/**
 * True only when the sample is confirmed to carry the NCBI/ClinVar ALT allele. Unknown
 * states (complex haplotype/diplotype records, FAIL calls, '=' reference-allele records)
 * are NOT carried, so highlights / "carried only" match the backend carried_count.
 */
export function carriesNcbiAlt(yourGenotype) {
  if (!yourGenotype) return false;
  if (yourGenotype.carries_ncbi_alt === true) return true;
  if (yourGenotype.carries_ncbi_alt === false) return false;
  return CARRIED_ZYGOSITIES.has(yourGenotype.zygosity_for_ncbi_alt) && yourGenotype.filter === 'PASS';
}

/** Neutral "could not be determined" state (complex record, low-quality call, ...). */
export function isIndeterminateGenotype(yourGenotype) {
  if (!yourGenotype || !yourGenotype.found) return false;
  const z = yourGenotype.zygosity_for_ncbi_alt;
  return z === 'unknown' || z === 'low_quality';
}

// ---------------------------------------------------------------------------
// Links
// ---------------------------------------------------------------------------

export function NcbiLink({ href, children, className = '', title, icon = true }) {
  const safe = safeHref(href);
  if (!safe) return <span className={className}>{children}</span>;
  return (
    <a
      href={safe}
      target="_blank"
      rel="noopener noreferrer"
      title={title}
      onClick={(e) => e.stopPropagation()}
      onKeyDown={(e) => e.stopPropagation()}
      className={`inline-flex items-center gap-1 hover:underline underline-offset-2 decoration-dotted ${className}`}
    >
      {children}
      {icon && <ExternalLink className="w-3 h-3 opacity-60 shrink-0" aria-hidden="true" />}
    </a>
  );
}

export function RsidLink({ rsid, className = 'text-cyan-300 hover:text-cyan-200' }) {
  if (!rsid) return null;
  if (!isRsid(String(rsid))) return <span className={className}>{rsid}</span>;
  return (
    <NcbiLink href={ncbiUrls.dbsnp(String(rsid))} className={className} title={`Open ${rsid} in NCBI dbSNP`}>
      {rsid}
    </NcbiLink>
  );
}

/** Links only genuine ClinVar RCV accessions; synthetic ids (RCV-CNV-..., VCF-DIRECT-READ) stay plain text. */
export function RcvLink({ rcv, className = 'text-slate-300 hover:text-cyan-200' }) {
  if (!rcv) return null;
  const id = String(rcv).trim();
  if (!isRealRcv(id)) return <span className={className}>{id}</span>;
  return (
    <NcbiLink href={ncbiUrls.clinvarRcv(id)} className={className} title={`Open ${id} in NCBI ClinVar`}>
      {id}
    </NcbiLink>
  );
}

export function RcvList({ ids, initial = 12, className = '' }) {
  const [expanded, setExpanded] = useState(false);
  const list = Array.isArray(ids) ? ids.filter(Boolean) : splitRcvIds(ids);
  if (list.length === 0) return <span className="text-slate-500">—</span>;
  const shown = expanded ? list : list.slice(0, initial);
  return (
    <div className={`flex flex-wrap items-center gap-x-2 gap-y-1 ${className}`}>
      {shown.map((id) => (
        <RcvLink key={id} rcv={id} className="text-cyan-300 hover:text-cyan-200 font-mono text-[11px]" />
      ))}
      {list.length > initial && (
        <button
          type="button"
          onClick={(e) => {
            e.stopPropagation();
            setExpanded((v) => !v);
          }}
          className="text-[11px] text-slate-400 hover:text-slate-200 underline decoration-dotted"
        >
          {expanded ? 'show fewer' : `+${list.length - initial} more`}
        </button>
      )}
    </div>
  );
}

/** Gene symbol(s) rendered as NCBI Gene search links (by symbol, human). */
export function GeneLinks({ gene, geneId, badge = true, className = '' }) {
  const symbols = splitGeneSymbols(gene);
  if (symbols.length === 0) {
    return gene ? (
      <span className={badge ? 'inline-block px-2 py-0.5 rounded-md bg-slate-800/90 border border-slate-700 font-mono font-bold text-white' : className}>
        {gene}
      </span>
    ) : null;
  }
  return (
    <span className={`inline-flex flex-wrap gap-1 ${className}`}>
      {symbols.map((sym) => (
        <NcbiLink
          key={sym}
          href={geneId && symbols.length === 1 ? ncbiUrls.geneId(geneId) : ncbiUrls.geneSymbol(sym)}
          title={`Open ${sym} in NCBI Gene`}
          icon={!badge}
          className={
            badge
              ? 'px-2 py-0.5 rounded-md bg-slate-800/90 border border-slate-700 font-mono font-bold text-white hover:border-cyan-500/60 hover:text-cyan-200'
              : 'font-mono text-cyan-300 hover:text-cyan-200'
          }
        >
          {sym}
        </NcbiLink>
      ))}
    </span>
  );
}

// ---------------------------------------------------------------------------
// Badges
// ---------------------------------------------------------------------------

export function classificationBadgeClass(classification) {
  const c = (classification || '').toLowerCase();
  if (c.includes('conflicting')) return 'bg-yellow-500/10 text-yellow-200 border-yellow-500/40';
  if (c.startsWith('pathogenic')) return 'bg-rose-500/15 text-rose-300 border-rose-500/40'; // incl. Pathogenic/Likely Pathogenic
  if (c.includes('likely pathogenic')) return 'bg-orange-500/15 text-orange-300 border-orange-500/40';
  if (c.includes('pathogenic')) return 'bg-rose-500/15 text-rose-300 border-rose-500/40';
  if (c.includes('risk factor') || c.includes('association')) return 'bg-amber-500/15 text-amber-300 border-amber-500/40';
  if (c.includes('drug response')) return 'bg-purple-500/15 text-purple-300 border-purple-500/40';
  if (c.includes('protective')) return 'bg-teal-500/15 text-teal-300 border-teal-500/40';
  if (c.includes('benign')) return 'bg-emerald-500/15 text-emerald-300 border-emerald-500/40';
  return 'bg-slate-700/50 text-slate-300 border-slate-600/60';
}

export function ClassificationBadge({ raw, normalized, size = 'sm' }) {
  const label = raw || normalized;
  if (!label) return <span className="text-slate-500 text-[11px]">not classified</span>;
  const pad = size === 'xs' ? 'px-1.5 py-0.5 text-[10px]' : 'px-2.5 py-1 text-[11px]';
  return (
    <span
      className={`inline-block rounded-lg border font-semibold ${pad} ${classificationBadgeClass(normalized || raw)}`}
      title={normalized && raw && normalized !== raw ? `Normalized: ${normalized}` : undefined}
    >
      {label}
    </span>
  );
}

export function ReviewStars({ stars, reviewStatus, showLabel = false }) {
  const n = Math.max(0, Math.min(4, Number(stars) || 0));
  return (
    <span className="inline-flex items-center gap-1.5" title={reviewStatus ? `${reviewStatus} (${n}/4)` : `${n} of 4 ClinVar review stars`}>
      <span className="inline-flex items-center gap-0.5">
        {[1, 2, 3, 4].map((idx) => (
          <Star
            key={idx}
            className={`w-3.5 h-3.5 ${idx <= n ? 'text-amber-400 fill-amber-400' : 'text-slate-700 fill-slate-800'}`}
          />
        ))}
      </span>
      {showLabel && reviewStatus && <span className="text-[11px] text-slate-400">{reviewStatus}</span>}
    </span>
  );
}

const ZYGOSITY_STYLES = {
  homozygous: { cls: 'bg-rose-500/20 text-rose-200 border-rose-500/40', label: 'Homozygous ALT' },
  heterozygous: { cls: 'bg-amber-500/20 text-amber-200 border-amber-500/40', label: 'Heterozygous' },
  hemizygous: { cls: 'bg-rose-500/20 text-rose-200 border-rose-500/40', label: 'Hemizygous ALT' },
  homoplasmic: { cls: 'bg-rose-500/20 text-rose-200 border-rose-500/40', label: 'Homoplasmic ALT' },
  absent: { cls: 'bg-slate-800 text-slate-400 border-slate-700', label: 'ALT absent' },
  reference_allele: { cls: 'bg-slate-800 text-slate-300 border-slate-700', label: 'Reference allele' },
  low_quality: { cls: 'bg-slate-700/40 text-amber-200/80 border-amber-500/30', label: 'FAIL – unconfirmed' },
  no_call: { cls: 'bg-fuchsia-500/15 text-fuchsia-200 border-fuchsia-500/40', label: 'No call' },
  unknown: { cls: 'bg-slate-700/40 text-slate-300 border-slate-600', label: 'Indeterminate' },
};

export function ZygosityBadge({ zygosity, className = '' }) {
  const style = ZYGOSITY_STYLES[zygosity] || ZYGOSITY_STYLES.unknown;
  return (
    <span className={`inline-block px-2 py-0.5 rounded-md border text-[10px] font-semibold uppercase tracking-wide ${style.cls} ${className}`}>
      {style.label}
    </span>
  );
}

/** Compact "your genotype" cell for result tables. */
export function GenotypeCell({ yourGenotype }) {
  if (!yourGenotype) {
    return <span className="text-[11px] text-slate-500">Not read</span>;
  }
  if (!yourGenotype.found) {
    return (
      <div className="space-y-1">
        <ZygosityBadge zygosity={yourGenotype.zygosity_for_ncbi_alt || 'unknown'} />
        <div className="text-[10px] text-slate-500">No VCF record at locus</div>
      </div>
    );
  }
  return (
    <div className="space-y-1">
      <div className="font-mono text-[11px] text-slate-100">
        {yourGenotype.alleles_called || '—'}
        {yourGenotype.genotype && <span className="text-cyan-400"> ({yourGenotype.genotype})</span>}
      </div>
      <ZygosityBadge zygosity={yourGenotype.zygosity_for_ncbi_alt || 'unknown'} />
      {yourGenotype.is_reference_block && <div className="text-[10px] text-slate-500">gVCF ref block</div>}
    </div>
  );
}

// ---------------------------------------------------------------------------
// States
// ---------------------------------------------------------------------------

export function ErrorsNotice({ errors, title = 'NCBI upstream notices' }) {
  const [open, setOpen] = useState(false);
  const list = asArray(errors).map(String);
  if (list.length === 0) return null;
  return (
    <div className="rounded-xl border border-amber-500/30 bg-amber-950/30 text-amber-200 text-xs">
      <button
        type="button"
        onClick={() => setOpen((v) => !v)}
        className="w-full flex items-center justify-between gap-2 px-3 py-2"
      >
        <span className="flex items-center gap-2">
          <AlertTriangle className="w-3.5 h-3.5" />
          {title}: {list.length} non-fatal issue{list.length === 1 ? '' : 's'} (results may be partial)
        </span>
        {open ? <ChevronUp className="w-3.5 h-3.5" /> : <ChevronDown className="w-3.5 h-3.5" />}
      </button>
      {open && (
        <ul className="px-4 pb-3 space-y-1 list-disc list-inside text-amber-200/90">
          {list.map((e, i) => (
            <li key={i} className="break-words">
              {e}
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}

export function ErrorBanner({ message, onRetry }) {
  if (!message) return null;
  return (
    <div className="p-3 rounded-xl bg-rose-950/50 border border-rose-500/40 text-rose-200 text-xs flex items-center justify-between gap-3">
      <span className="break-words">{message}</span>
      {onRetry && (
        <button
          type="button"
          onClick={onRetry}
          className="px-3 py-1 rounded-lg bg-rose-500/20 hover:bg-rose-500/30 text-xs font-semibold shrink-0"
        >
          Retry
        </button>
      )}
    </div>
  );
}

export function SkeletonCards({ count = 3 }) {
  return (
    <div className="space-y-3 animate-pulse" aria-busy="true" aria-label="Loading NCBI data">
      {Array.from({ length: count }).map((_, i) => (
        <div key={i} className="rounded-xl border border-slate-800 bg-slate-950 p-4 space-y-2">
          <div className="h-3 w-1/3 bg-slate-800 rounded" />
          <div className="h-3 w-2/3 bg-slate-800/80 rounded" />
          <div className="h-3 w-1/2 bg-slate-800/60 rounded" />
        </div>
      ))}
    </div>
  );
}

export function EmptyState({ title, children }) {
  return (
    <div className="p-6 rounded-xl border border-dashed border-slate-800 text-center text-slate-400 text-xs space-y-1">
      <p className="text-sm font-semibold text-slate-200">{title}</p>
      {children && <div>{children}</div>}
    </div>
  );
}

function SectionTitle({ icon: Icon, children, right }) {
  return (
    <div className="flex items-center justify-between gap-2 mb-2">
      <h4 className="text-xs font-semibold uppercase tracking-wider text-slate-400 flex items-center gap-1.5">
        {Icon && <Icon className="w-3.5 h-3.5 text-cyan-400" />}
        {children}
      </h4>
      {right}
    </div>
  );
}

function Field({ label, children }) {
  return (
    <div>
      <div className="text-[10px] uppercase tracking-wider text-slate-500">{label}</div>
      <div className="text-xs text-slate-200 mt-0.5">{children}</div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// ClinVar
// ---------------------------------------------------------------------------

export function ClinVarCard({ record }) {
  if (!record) return null;
  const rcvs = asArray(record.rcv_accessions);
  const conditions = asArray(record.conditions);
  const genes = asArray(record.genes);
  const consequences = asArray(record.molecular_consequences);
  const vcvUrl = safeHref(record.url) || (record.variation_id ? ncbiUrls.clinvarVariation(record.variation_id) : null);
  const locus = formatLocus(record.chr, record.start, record.stop);

  return (
    <div className="rounded-xl border border-slate-800 bg-slate-950 p-4 space-y-3">
      <div className="flex flex-wrap items-start justify-between gap-2">
        <div className="min-w-0">
          <div className="text-sm font-semibold text-white break-words">{record.title || `ClinVar variation ${record.variation_id || ''}`}</div>
          <div className="flex flex-wrap items-center gap-2 mt-1 text-[11px] font-mono">
            {(record.vcv_id || record.variation_id) && (
              <NcbiLink href={vcvUrl} className="text-cyan-300 hover:text-cyan-200" title="Open ClinVar variation page">
                {record.vcv_id || `Variation ${record.variation_id}`}
              </NcbiLink>
            )}
            {record.rsid && <RsidLink rsid={record.rsid} />}
            {record.variant_type && <span className="text-slate-400 font-sans">{record.variant_type}</span>}
          </div>
        </div>
        <ClassificationBadge raw={record.classification} normalized={record.classification_normalized} />
      </div>
      {record.somatic_classification && (
        <p className="text-[11px] text-slate-400">
          <span className="text-slate-500 uppercase tracking-wider text-[10px] mr-1.5">Somatic</span>
          {record.somatic_classification}
        </p>
      )}

      <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
        <Field label="Review status">
          <div className="flex flex-col gap-0.5">
            <ReviewStars stars={record.review_stars} reviewStatus={record.review_status} />
            <span className="text-[11px] text-slate-400">{record.review_status || '—'}</span>
          </div>
        </Field>
        <Field label="Last evaluated">{record.last_evaluated || '—'}</Field>
        <Field label={`Locus (${record.assembly || 'GRCh38'})`}>
          <span className="font-mono">{locus || '—'}</span>
          {(record.ref || record.alt) && (
            <span className="font-mono ml-1.5">
              {record.ref || '?'} &rarr; <span className="text-rose-300">{record.alt || '?'}</span>
            </span>
          )}
        </Field>
        {record.protein_change && (
          <Field label="Protein change">
            <span className="font-mono">{record.protein_change}</span>
          </Field>
        )}
        {genes.length > 0 && (
          <Field label="Gene(s)">
            <span className="flex flex-wrap gap-1.5">
              {genes.map((g, i) => (
                <NcbiLink
                  key={`${g.symbol || g.gene_id}-${i}`}
                  href={g.gene_id ? ncbiUrls.geneId(g.gene_id) : g.symbol ? ncbiUrls.geneSymbol(g.symbol) : null}
                  className="font-mono text-cyan-300 hover:text-cyan-200"
                >
                  {g.symbol || g.gene_id}
                </NcbiLink>
              ))}
            </span>
          </Field>
        )}
      </div>

      {consequences.length > 0 && (
        <div className="flex flex-wrap gap-1.5">
          {consequences.map((mc, i) => (
            <span key={`${mc}-${i}`} className="px-2 py-0.5 rounded-md bg-slate-800 border border-slate-700 text-[10px] text-slate-300">
              {mc}
            </span>
          ))}
        </div>
      )}

      <div>
        <div className="text-[10px] uppercase tracking-wider text-slate-500 mb-1">Conditions</div>
        {conditions.length === 0 ? (
          <span className="text-xs text-slate-500">No condition reported</span>
        ) : (
          <ul className="space-y-1">
            {conditions.map((c, i) => (
              <li key={`${c.name}-${i}`} className="text-xs text-slate-200 flex flex-wrap items-center gap-x-2 gap-y-0.5">
                <span className="font-medium">{c.name || 'Unnamed condition'}</span>
                {c.medgen_id && (
                  <NcbiLink href={safeHref(c.url) || ncbiUrls.medgen(c.medgen_id)} className="text-[10px] font-mono text-cyan-400">
                    MedGen {c.medgen_id}
                  </NcbiLink>
                )}
                {!c.medgen_id && safeHref(c.url) && (
                  <NcbiLink href={c.url} className="text-[10px] text-cyan-400">
                    details
                  </NcbiLink>
                )}
                {asArray(c.omim_ids).map((mim) => (
                  <NcbiLink key={mim} href={ncbiUrls.omim(mim)} className="text-[10px] font-mono text-indigo-300">
                    OMIM {mim}
                  </NcbiLink>
                ))}
              </li>
            ))}
          </ul>
        )}
      </div>

      <div>
        <div className="text-[10px] uppercase tracking-wider text-slate-500 mb-1">
          RCV accessions ({rcvs.length})
        </div>
        <RcvList ids={rcvs} />
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// dbSNP frequencies
// ---------------------------------------------------------------------------

export function FrequencyPanel({ dbsnp, riskAlleleFrequency, riskAllele }) {
  if (!dbsnp) return null;
  const freqs = asArray(dbsnp.frequencies).filter((f) => f && f.frequency !== null && f.frequency !== undefined && !Number.isNaN(Number(f.frequency)));
  const targetAllele = riskAlleleFrequency?.allele || riskAllele || asArray(dbsnp.alts)[0] || null;
  const matching = targetAllele ? freqs.filter((f) => f.allele === targetAllele) : [];
  const rows = matching.length > 0 ? matching : freqs;
  const showingOtherAlleles = matching.length === 0 && freqs.length > 0 && !!targetAllele;
  const preferredStudy = riskAlleleFrequency?.study || dbsnp.preferred_frequency?.study || null;
  const maxFreq = rows.reduce((m, f) => Math.max(m, Number(f.frequency)), 0);
  const scaleMax = maxFreq > 0 ? Math.min(1, maxFreq * 1.1) : 1;

  return (
    <div className="rounded-xl border border-slate-800 bg-slate-950 p-4 space-y-3">
      <SectionTitle
        icon={Library}
        right={dbsnp.rsid && <RsidLink rsid={dbsnp.rsid} className="text-[11px] font-mono text-cyan-300 hover:text-cyan-200" />}
      >
        dbSNP population frequencies
      </SectionTitle>

      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
        <Field label="Placement">
          <span className="font-mono">{formatLocus(dbsnp.chr, dbsnp.position) || '—'}</span>
        </Field>
        <Field label="Alleles">
          <span className="font-mono">
            {dbsnp.ref || '?'} &rarr; <span className="text-rose-300">{asArray(dbsnp.alts).join(', ') || '?'}</span>
          </span>
        </Field>
        <Field label="Class">{dbsnp.variant_class || '—'}</Field>
        <Field label={`Risk allele ${targetAllele ? `(${targetAllele})` : ''}`}>
          {riskAlleleFrequency ? (
            <span className="font-mono text-cyan-300">
              {formatFrequency(riskAlleleFrequency.frequency)}
              <span className="text-slate-500 font-sans"> · {riskAlleleFrequency.study}</span>
            </span>
          ) : (
            <span className="text-slate-500">not reported</span>
          )}
        </Field>
      </div>

      {(asArray(dbsnp.function_class).length > 0 || asArray(dbsnp.clinical_significance).length > 0) && (
        <div className="flex flex-wrap gap-1.5">
          {asArray(dbsnp.function_class).map((fc, i) => (
            <span key={`fc-${i}`} className="px-2 py-0.5 rounded-md bg-slate-800 border border-slate-700 text-[10px] text-slate-300">
              {fc}
            </span>
          ))}
          {asArray(dbsnp.clinical_significance).map((cs, i) => (
            <span key={`cs-${i}`} className={`px-2 py-0.5 rounded-md border text-[10px] ${classificationBadgeClass(cs)}`}>
              {cs}
            </span>
          ))}
        </div>
      )}

      {rows.length === 0 ? (
        <p className="text-xs text-slate-500">No per-study allele frequencies reported by dbSNP.</p>
      ) : (
        <div className="space-y-1.5">
          <div className="flex items-center justify-between text-[10px] text-slate-500">
            <span>
              {showingOtherAlleles
                ? freqs.some((f) => f.allele === '-' || String(f.allele || '').length !== String(targetAllele).length)
                  ? `dbSNP spells indel alleles differently from VCF ${targetAllele} (e.g. "-"); showing reported alleles`
                  : `No study reports allele ${targetAllele}; showing reported alleles`
                : `Frequency of allele ${targetAllele || 'ALT'} by study`}
            </span>
            <span>bar scale 0 – {formatFrequency(scaleMax)}</span>
          </div>
          {rows.map((f, i) => {
            const val = Number(f.frequency);
            const pct = scaleMax > 0 ? Math.max(0, Math.min(100, (val / scaleMax) * 100)) : 0;
            const isPreferred = preferredStudy && f.study === preferredStudy;
            return (
              <div
                key={`${f.study}-${f.allele}-${i}`}
                className={`grid grid-cols-[minmax(0,7.5rem)_1fr_auto] items-center gap-2 px-2 py-1 rounded-lg ${
                  isPreferred ? 'bg-cyan-500/10 ring-1 ring-cyan-500/40' : ''
                }`}
              >
                <span className={`text-[11px] truncate ${isPreferred ? 'text-cyan-200 font-semibold' : 'text-slate-300'}`} title={f.study}>
                  {f.study || 'unknown study'}
                  {showingOtherAlleles && f.allele ? <span className="font-mono text-slate-500"> ({f.allele})</span> : null}
                </span>
                <div className="h-2 bg-slate-800 rounded-full overflow-hidden" aria-hidden="true">
                  <div
                    className={`h-full rounded-full ${isPreferred ? 'bg-cyan-400' : 'bg-slate-500'}`}
                    style={{ width: `${pct}%`, minWidth: val > 0 ? '2px' : 0 }}
                  />
                </div>
                <span className="text-[11px] font-mono text-slate-200 text-right whitespace-nowrap">
                  {formatFrequency(val)}
                  {formatNumber(f.count) && <span className="text-slate-500"> · n={formatNumber(f.count)}</span>}
                  {isPreferred && <span className="ml-1 text-[9px] uppercase text-cyan-300">preferred</span>}
                </span>
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Your genome vs NCBI
// ---------------------------------------------------------------------------

export function YourGenomeBox({ yourGenotype, ncbiAlt, ncbiRef }) {
  const [showRaw, setShowRaw] = useState(false);
  if (!yourGenotype) {
    return (
      <div className="rounded-xl border border-slate-800 bg-slate-950 p-4">
        <SectionTitle icon={User}>Your genome vs NCBI</SectionTitle>
        <p className="text-xs text-slate-500">
          Your VCF was not read for this record (no resolvable GRCh38 locus).
        </p>
      </div>
    );
  }
  const yg = yourGenotype;
  const carried = carriesNcbiAlt(yg);
  const zyg = yg.zygosity_for_ncbi_alt || 'unknown';

  const altLabel = ncbiAlt ? ` (${ncbiAlt})` : '';
  let verdict;
  if (!yg.found) verdict = 'No record in your VCF files at this locus (no coverage or not called).';
  else if (zyg === 'homozygous') verdict = `You carry two copies of the NCBI ALT allele${altLabel}.`;
  else if (zyg === 'heterozygous') verdict = `You carry one copy of the NCBI ALT allele${altLabel}.`;
  else if (zyg === 'hemizygous')
    verdict = `You carry the NCBI ALT allele${altLabel} on your single copy of this chromosome (hemizygous).`;
  else if (zyg === 'homoplasmic') verdict = `Your mitochondrial DNA carries the NCBI ALT allele${altLabel} (homoplasmic call).`;
  else if (zyg === 'reference_allele')
    verdict = 'This ClinVar record describes the reference allele itself (“=”, no change); your call includes the reference allele. It is not an ALT allele you carry.';
  else if (zyg === 'low_quality')
    verdict = `Your VCF call here did not pass quality filters (FILTER ${yg.filter || '?'}), so it can neither confirm nor rule out the NCBI allele.`;
  else if (zyg === 'absent')
    verdict = yg.is_carried
      ? 'You carry a different ALT allele at this position, not the NCBI/ClinVar allele.'
      : 'You do not carry the NCBI ALT allele here (reference call).';
  else if (zyg === 'no_call') verdict = 'The genotype at this locus is a no-call (./.).';
  else
    verdict =
      'Indeterminate: this record (e.g. a haplotype / diplotype) needs several alleles and their phase, which a single-locus read cannot confirm.';

  return (
    <div
      className={`rounded-xl border p-4 space-y-3 ${
        carried ? 'border-amber-500/40 bg-amber-950/20' : 'border-slate-800 bg-slate-950'
      }`}
    >
      <SectionTitle icon={User} right={<ZygosityBadge zygosity={zyg} />}>
        Your genome vs NCBI
      </SectionTitle>
      <p className={`text-xs ${carried ? 'text-amber-100' : 'text-slate-300'}`}>{verdict}</p>

      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
        <Field label="Locus read">
          <span className="font-mono">{formatLocus(yg.chr, yg.position) || '—'}</span>
        </Field>
        <Field label="Your VCF REF > ALT">
          <span className="font-mono">
            {yg.ref || '—'} &rarr; {yg.alt || '—'}
          </span>
        </Field>
        <Field label="NCBI REF > ALT">
          <span className="font-mono">
            {ncbiRef || '?'} &rarr; <span className="text-rose-300">{ncbiAlt || '?'}</span>
          </span>
        </Field>
        <Field label="Called alleles (GT)">
          <span className="font-mono font-semibold text-white">
            {yg.alleles_called || '—'}
            {yg.genotype && <span className="text-cyan-400"> ({yg.genotype})</span>}
          </span>
        </Field>
        <Field label="Your status">{yg.your_status || '—'}</Field>
        <Field label="Carries NCBI ALT">
          {yg.carries_ncbi_alt === true ? (
            <span className="text-amber-300 font-semibold">Yes</span>
          ) : yg.carries_ncbi_alt === false ? (
            <span className="text-slate-300">No</span>
          ) : (
            <span className="text-slate-500">Unknown</span>
          )}
        </Field>
        <Field label="DP / GQ">
          <span className="font-mono">
            {yg.dp ?? '—'} / {yg.gq ?? '—'}
          </span>
        </Field>
        <Field label="AD / FILTER">
          <span className="font-mono">
            {yg.ad || '—'} /{' '}
            <span className={yg.filter === 'PASS' ? 'text-emerald-300' : 'text-amber-300'}>{yg.filter || '—'}</span>
          </span>
        </Field>
      </div>

      {yg.is_reference_block && (
        <p className="text-[11px] text-slate-400 bg-slate-900 border border-slate-800 rounded-lg px-3 py-2">
          This position falls inside a gVCF homozygous-reference block (END=…). The caller saw only reference
          reads across the block, so DP/GQ are block-level summaries rather than a site-specific call.
        </p>
      )}

      {(yg.raw_vcf_line || yg.source_file) && (
        <div>
          <button
            type="button"
            onClick={() => setShowRaw((v) => !v)}
            className="text-[11px] text-slate-400 hover:text-slate-200 flex items-center gap-1"
          >
            {showRaw ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
            Raw VCF line{yg.source_file ? ` (${yg.source_file})` : ''}
          </button>
          {showRaw && (
            <pre className="mt-1.5 p-2.5 rounded-lg bg-slate-900 border border-slate-800 text-[11px] font-mono text-cyan-200 overflow-x-auto">
              {yg.raw_vcf_line || '(no line)'}
            </pre>
          )}
        </div>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Gene & publications
// ---------------------------------------------------------------------------

export function GeneCard({ gene }) {
  const [expanded, setExpanded] = useState(false);
  if (!gene) return null;
  const aliases = asArray(gene.aliases);
  const mims = asArray(gene.mim_ids);
  const href = safeHref(gene.url) || (gene.gene_id ? ncbiUrls.geneId(gene.gene_id) : gene.symbol ? ncbiUrls.geneSymbol(gene.symbol) : null);
  const summary = typeof gene.summary === 'string' ? gene.summary.trim() : '';
  const longSummary = summary.length > 320;
  return (
    <div className="rounded-xl border border-slate-800 bg-slate-950 p-4 space-y-3">
      <SectionTitle icon={Dna} right={gene.gene_id && <span className="text-[10px] font-mono text-slate-500">Gene ID {gene.gene_id}</span>}>
        NCBI Gene
      </SectionTitle>
      <div className="flex flex-wrap items-baseline gap-2">
        <NcbiLink href={href} className="text-base font-bold font-mono text-white hover:text-cyan-200">
          {gene.symbol || gene.gene_id || 'Gene'}
        </NcbiLink>
        {gene.name && <span className="text-xs text-slate-300">{gene.name}</span>}
      </div>
      <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
        <Field label="Cytoband">{gene.map_location || '—'}</Field>
        <Field label="GRCh38 span">
          <span className="font-mono">{formatLocus(gene.chr, gene.start, gene.stop) || '—'}</span>
        </Field>
        <Field label="OMIM (MIM)">
          {mims.length === 0 ? (
            '—'
          ) : (
            <span className="flex flex-wrap gap-1.5">
              {mims.map((m) => (
                <NcbiLink key={m} href={ncbiUrls.omim(m)} className="font-mono text-indigo-300">
                  {m}
                </NcbiLink>
              ))}
            </span>
          )}
        </Field>
      </div>
      {aliases.length > 0 && (
        <div className="flex flex-wrap items-center gap-1.5">
          <span className="text-[10px] uppercase tracking-wider text-slate-500 mr-1">Aliases</span>
          {aliases.map((a, i) => (
            <span key={`${a}-${i}`} className="px-1.5 py-0.5 rounded bg-slate-800 border border-slate-700 text-[10px] font-mono text-slate-300">
              {a}
            </span>
          ))}
        </div>
      )}
      {summary && (
        <div>
          <p className={`text-xs text-slate-300 leading-relaxed whitespace-pre-line ${expanded || !longSummary ? '' : 'line-clamp-3'}`}>
            {summary}
          </p>
          {longSummary && (
            <button
              type="button"
              onClick={() => setExpanded((v) => !v)}
              className="mt-1 text-[11px] text-cyan-400 hover:text-cyan-300"
            >
              {expanded ? 'Collapse summary' : 'Read full RefSeq summary'}
            </button>
          )}
        </div>
      )}
    </div>
  );
}

export function PublicationsList({ publications, pubmedLink }) {
  const pubs = asArray(publications).slice(0, 10);
  return (
    <div className="rounded-xl border border-slate-800 bg-slate-950 p-4 space-y-2">
      <SectionTitle
        icon={BookOpen}
        right={
          safeHref(pubmedLink) && (
            <NcbiLink href={pubmedLink} className="text-[11px] text-cyan-400">
              All on PubMed
            </NcbiLink>
          )
        }
      >
        PubMed publications ({pubs.length})
      </SectionTitle>
      {pubs.length === 0 ? (
        <p className="text-xs text-slate-500">No linked publications returned.</p>
      ) : (
        <ol className="space-y-2">
          {pubs.map((p, i) => (
            <li key={p.pmid || i} className="text-xs">
              <NcbiLink
                href={safeHref(p.url) || (p.pmid ? ncbiUrls.pubmed(p.pmid) : null)}
                className="text-slate-100 hover:text-cyan-200 font-medium leading-snug"
              >
                {p.title || `PMID ${p.pmid}`}
              </NcbiLink>
              <div className="text-[11px] text-slate-400 mt-0.5">
                {[p.authors, p.journal, p.year].filter(Boolean).join(' · ')}
                {p.pmid && <span className="font-mono text-slate-500"> · PMID {p.pmid}</span>}
              </div>
            </li>
          ))}
        </ol>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Composite: full NCBI variant body (shared by inspector drawer & lookup form)
// ---------------------------------------------------------------------------

export function NcbiVariantDetail({ data, showPublications = true }) {
  if (!data) return null;
  const clinvar = asArray(data.clinvar);
  const primary = clinvar[0] || null;
  const links = data.links || {};
  const resolved = data.resolved || {};
  const nothing = clinvar.length === 0 && !data.dbsnp && !data.gene;

  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center justify-between gap-2 text-[11px]">
        <div className="flex flex-wrap items-center gap-2 text-slate-400">
          <span>Resolved:</span>
          <span className="font-mono text-slate-200">{formatLocus(resolved.chr, resolved.position) || 'no GRCh38 locus'}</span>
          {resolved.rsid && <RsidLink rsid={resolved.rsid} className="font-mono text-cyan-300 hover:text-cyan-200" />}
        </div>
        <div className="flex flex-wrap items-center gap-3">
          {safeHref(links.clinvar) && (
            <NcbiLink href={links.clinvar} className="text-cyan-400">
              ClinVar
            </NcbiLink>
          )}
          {safeHref(links.dbsnp) && (
            <NcbiLink href={links.dbsnp} className="text-cyan-400">
              dbSNP
            </NcbiLink>
          )}
          {safeHref(links.gene) && (
            <NcbiLink href={links.gene} className="text-cyan-400">
              Gene
            </NcbiLink>
          )}
          {safeHref(links.pubmed) && (
            <NcbiLink href={links.pubmed} className="text-cyan-400">
              PubMed
            </NcbiLink>
          )}
        </div>
      </div>

      <ErrorsNotice errors={data.errors} />

      {nothing ? (
        <EmptyState title="NCBI returned no ClinVar, dbSNP or Gene record for this query.">
          The variant may be absent from ClinVar, or the identifier may not map to GRCh38.
        </EmptyState>
      ) : (
        <>
          <YourGenomeBox
            yourGenotype={data.your_genotype}
            ncbiAlt={primary?.alt || data.risk_allele_frequency?.allele || null}
            ncbiRef={primary?.ref || data.dbsnp?.ref || null}
          />

          <div>
            <SectionTitle icon={FlaskConical}>ClinVar ({clinvar.length})</SectionTitle>
            {clinvar.length === 0 ? (
              <p className="text-xs text-slate-500">No ClinVar record at this locus.</p>
            ) : (
              <div className="space-y-3">
                {clinvar.map((rec, i) => (
                  <ClinVarCard key={rec.variation_id || i} record={rec} />
                ))}
              </div>
            )}
          </div>

          <FrequencyPanel dbsnp={data.dbsnp} riskAlleleFrequency={data.risk_allele_frequency} riskAllele={primary?.alt} />
          <GeneCard gene={data.gene} />
          {showPublications && <PublicationsList publications={data.publications} pubmedLink={links.pubmed} />}
        </>
      )}

      {data.fetched_at && (
        <p className="text-[10px] text-slate-500 text-right">Fetched from NCBI {formatFetchedAt(data.fetched_at)}</p>
      )}
    </div>
  );
}
