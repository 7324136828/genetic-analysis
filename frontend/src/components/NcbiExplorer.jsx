import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  ArrowRight,
  CheckCircle2,
  Crosshair,
  DatabaseZap,
  Dna,
  Download,
  Globe,
  Loader2,
  Search,
} from 'lucide-react';
import { fetchNcbiGene, fetchNcbiVariant, importNcbiVariants, searchNcbi } from '../api';
import {
  ClassificationBadge,
  EmptyState,
  ErrorBanner,
  ErrorsNotice,
  GENE_SYMBOL_RE,
  GeneCard,
  GenotypeCell,
  NcbiLink,
  NcbiVariantDetail,
  ReviewStars,
  RsidLink,
  SkeletonCards,
  asArray,
  carriesNcbiAlt,
  formatFetchedAt,
  formatLocus,
  isRealRcv,
  isRsid,
  ncbiUrls,
  safeHref,
} from './ncbiShared';

const SEARCH_EXAMPLES = ['RCV000019449', 'rs671', 'VCV000017865', 'HFE', 'hereditary hemochromatosis'];
const GENE_EXAMPLES = ['ALDH2', 'HFE', 'CYP2C19', 'BRCA1', 'MTHFR', 'APOE'];
const LOOKUP_EXAMPLES = ['rs671', 'RCV000019449', 'VCV000017865', '12:111803962'];

const SIGNIFICANCE_OPTIONS = [
  { value: 'pathogenic', label: 'Pathogenic + Likely pathogenic' },
  { value: 'risk', label: 'Risk factor / Association' },
  { value: 'drug', label: 'Drug response' },
  { value: 'all', label: 'All classifications' },
];

const inputCls =
  'w-full px-3 py-2 bg-slate-950 border border-slate-700/80 rounded-xl text-xs text-slate-100 placeholder-slate-500 focus:outline-none focus:border-cyan-500 focus:ring-1 focus:ring-cyan-500';
const selectCls =
  'px-3 py-2 bg-slate-950 border border-slate-700/80 rounded-xl text-xs text-slate-200 focus:outline-none focus:border-cyan-500';
const primaryBtn =
  'px-4 py-2 rounded-xl bg-cyan-500 hover:bg-cyan-400 text-slate-950 font-semibold text-xs transition flex items-center justify-center gap-1.5 disabled:opacity-50 disabled:cursor-not-allowed';

/** Parse the quick-lookup box into /api/ncbi/variant query params. */
export function parseLookupInput(raw) {
  const text = String(raw || '').trim();
  if (!text) return { error: 'Enter an rsID, RCV, VCV or chr:position.' };
  if (/^rs\d+$/i.test(text)) return { params: { rsid: text.toLowerCase() } };
  if (/^RCV\d+(\.\d+)?$/i.test(text)) return { params: { rcv_id: text.toUpperCase() } };
  if (/^VCV\d+(\.\d+)?$/i.test(text)) return { params: { vcv_id: text.toUpperCase() } };
  const m = text.match(/^(?:chr)?([0-9]{1,2}|X|Y|MT|M)\s*[:\s]\s*([\d,]+)$/i);
  if (m) {
    const chrRaw = m[1].toUpperCase();
    const chr = chrRaw === 'M' ? 'MT' : chrRaw;
    const position = parseInt(m[2].replace(/,/g, ''), 10);
    const chrNum = Number(chr);
    if (!Number.isNaN(chrNum) && (chrNum < 1 || chrNum > 22)) return { error: `Unknown chromosome "${m[1]}".` };
    if (!position || position < 1) return { error: 'Position must be a positive integer.' };
    return { params: { chr, position } };
  }
  return { error: `"${text}" is not an rsID (rs671), RCV (RCV000019449), VCV (VCV000017865) or chr:position (12:111803962).` };
}

/** Filter patch for jumping to the local catalog Explorer. */
function explorerPatchFor(clinvar) {
  if (!clinvar) return null;
  if (clinvar.rsid && isRsid(clinvar.rsid)) return { variant_id: clinvar.rsid };
  const rcv = asArray(clinvar.rcv_accessions).find(isRealRcv);
  if (rcv) return { rcv_id: rcv };
  const gene = asArray(clinvar.genes)[0]?.symbol;
  if (gene) return { gene };
  return null;
}

function lookupParamsFor(clinvar) {
  if (!clinvar) return null;
  if (clinvar.vcv_id) return { vcv_id: clinvar.vcv_id };
  const rcv = asArray(clinvar.rcv_accessions).find(isRealRcv);
  if (rcv) return { rcv_id: rcv };
  if (clinvar.rsid) return { rsid: clinvar.rsid };
  if (clinvar.chr && clinvar.start) return { chr: clinvar.chr, position: clinvar.start };
  return null;
}

function lookupLabelFor(params) {
  if (!params) return '';
  return params.vcv_id || params.rcv_id || params.rsid || (params.chr ? `${params.chr}:${params.position}` : '');
}

// ---------------------------------------------------------------------------
// Selection + import (shared by ClinVar search and Gene scan)
// ---------------------------------------------------------------------------

function useNcbiImport(onCatalogChanged) {
  const onChangedRef = useRef(onCatalogChanged);
  useEffect(() => {
    onChangedRef.current = onCatalogChanged;
  }, [onCatalogChanged]);
  const [selected, setSelected] = useState(() => new Set());
  const [importing, setImporting] = useState(false);
  const [notice, setNotice] = useState(null);
  const [importError, setImportError] = useState(null);
  const [importedIds, setImportedIds] = useState(() => new Set());

  const toggle = (id) =>
    setSelected((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });

  const setMany = (ids, on) =>
    setSelected((prev) => {
      const next = new Set(prev);
      ids.forEach((id) => (on ? next.add(id) : next.delete(id)));
      return next;
    });

  const reset = () => {
    setSelected(new Set());
    setNotice(null);
    setImportError(null);
  };

  const importSelected = async (items) => {
    const ids = Array.from(selected).slice(0, 50);
    if (ids.length === 0) return;
    setImporting(true);
    setImportError(null);
    setNotice(null);
    try {
      const res = await importNcbiVariants(ids);
      const skipped = asArray(res.skipped);
      const skippedIds = new Set(skipped.map((s) => String(s.variation_id)));
      const okIds = ids.filter((id) => !skippedIds.has(String(id)));
      const importedItems = asArray(items)
        .map((it) => it?.clinvar)
        .filter((cv) => cv && okIds.includes(String(cv.variation_id)));
      setImportedIds((prev) => {
        const next = new Set(prev);
        okIds.forEach((id) => next.add(String(id)));
        return next;
      });
      setSelected(new Set());
      setNotice({
        imported: Number(res.imported) || 0,
        recordIds: asArray(res.record_ids),
        skipped,
        errors: asArray(res.errors),
        importedItems,
        fetchedAt: res.fetched_at,
      });
      if ((Number(res.imported) || 0) > 0 && onChangedRef.current) onChangedRef.current();
    } catch (err) {
      setImportError(err.message || 'Import failed');
    } finally {
      setImporting(false);
    }
  };

  return { selected, toggle, setMany, reset, importing, notice, importError, importSelected, importedIds, setNotice };
}

function ImportBar({ importer, items, onOpenInExplorer }) {
  const { selected, importing, notice, importError, importSelected, setNotice } = importer;
  const count = selected.size;

  let jumpPatch = null;
  let jumpLabel = '';
  if (notice && notice.importedItems.length > 0) {
    const genes = new Set(notice.importedItems.map((cv) => asArray(cv.genes)[0]?.symbol).filter(Boolean));
    if (notice.importedItems.length > 1 && genes.size === 1) {
      jumpPatch = { gene: Array.from(genes)[0] };
      jumpLabel = `gene ${jumpPatch.gene}`;
    } else {
      jumpPatch = explorerPatchFor(notice.importedItems[0]);
      jumpLabel = jumpPatch ? Object.values(jumpPatch)[0] : '';
    }
  }

  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-center justify-between gap-2">
        <span className="text-xs text-slate-400">
          {count > 0 ? (
            <>
              <span className="font-semibold text-white">{count}</span> selected
              {count > 50 && <span className="text-amber-300"> (only the first 50 are imported per request)</span>}
            </>
          ) : (
            'Select ClinVar records to add them (with your VCF read) to the local catalog.'
          )}
        </span>
        <button
          type="button"
          disabled={count === 0 || importing}
          onClick={() => importSelected(items)}
          className={primaryBtn}
        >
          {importing ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : <Download className="w-3.5 h-3.5" />}
          {importing ? 'Importing…' : `Import selected into catalog${count ? ` (${Math.min(count, 50)})` : ''}`}
        </button>
      </div>

      <ErrorBanner message={importError} />

      {notice && (
        <div className="p-3 rounded-xl border border-emerald-500/30 bg-emerald-950/30 text-emerald-100 text-xs space-y-2">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <span className="flex items-center gap-2 font-semibold">
              <CheckCircle2 className="w-4 h-4 text-emerald-300" />
              Imported {notice.imported} record{notice.imported === 1 ? '' : 's'} into the local catalog
              {notice.skipped.length > 0 && <span className="text-amber-200 font-normal">· {notice.skipped.length} skipped</span>}
            </span>
            <div className="flex items-center gap-2">
              {jumpPatch && onOpenInExplorer && (
                <button
                  type="button"
                  onClick={() => onOpenInExplorer(jumpPatch)}
                  className="px-3 py-1 rounded-lg bg-emerald-500/20 hover:bg-emerald-500/30 border border-emerald-500/40 text-emerald-100 font-semibold flex items-center gap-1"
                >
                  Open in Explorer ({jumpLabel})
                  <ArrowRight className="w-3.5 h-3.5" />
                </button>
              )}
              <button type="button" onClick={() => setNotice(null)} className="text-emerald-300/70 hover:text-emerald-100">
                Dismiss
              </button>
            </div>
          </div>
          {notice.skipped.length > 0 && (
            <ul className="list-disc list-inside text-amber-200/90 space-y-0.5">
              {notice.skipped.map((s, i) => (
                <li key={`${s.variation_id}-${i}`}>
                  Variation <span className="font-mono">{s.variation_id}</span>: {s.reason || 'skipped'}
                </li>
              ))}
            </ul>
          )}
          <ErrorsNotice errors={notice.errors} title="Import notices" />
        </div>
      )}
    </div>
  );
}

function VariantResultsTable({ items, importer, onDetails, onOpenInExplorer, emptyTitle, emptyHint }) {
  const list = asArray(items).filter((it) => it && it.clinvar);
  const selectableIds = list.map((it) => it.clinvar.variation_id).filter(Boolean).map(String);
  const allSelected = selectableIds.length > 0 && selectableIds.every((id) => importer.selected.has(id));

  if (list.length === 0) {
    return <EmptyState title={emptyTitle || 'No ClinVar records'}>{emptyHint}</EmptyState>;
  }

  return (
    <div className="rounded-2xl border border-slate-800 overflow-hidden">
      <div className="overflow-x-auto">
        <table className="w-full text-left text-xs">
          <thead className="bg-slate-950 text-[10px] uppercase tracking-wider text-slate-400 border-b border-slate-800">
            <tr>
              <th className="py-3 px-3 w-8">
                <input
                  type="checkbox"
                  aria-label="Select all"
                  checked={allSelected}
                  onChange={(e) => importer.setMany(selectableIds, e.target.checked)}
                  className="rounded border-slate-600 bg-slate-900 text-cyan-500 focus:ring-cyan-500"
                />
              </th>
              <th className="py-3 px-3">ClinVar variant / VCV</th>
              <th className="py-3 px-3">Gene</th>
              <th className="py-3 px-3">GRCh38 locus</th>
              <th className="py-3 px-3">rsID</th>
              <th className="py-3 px-3">Classification</th>
              <th className="py-3 px-3">Review</th>
              <th className="py-3 px-3">Conditions</th>
              <th className="py-3 px-3">Your genotype</th>
              <th className="py-3 px-3 text-right">Actions</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-800/60">
            {list.map((it, idx) => {
              const cv = it.clinvar;
              const vid = cv.variation_id ? String(cv.variation_id) : null;
              const carried = carriesNcbiAlt(it.your_genotype);
              const conditions = asArray(cv.conditions);
              const genes = asArray(cv.genes);
              const imported = vid && importer.importedIds.has(vid);
              const patch = explorerPatchFor(cv);
              return (
                <tr
                  key={vid || idx}
                  className={`align-top transition-colors ${
                    carried ? 'bg-amber-500/[0.07] hover:bg-amber-500/10' : 'hover:bg-slate-800/40'
                  }`}
                >
                  <td className={`py-3 px-3 ${carried ? 'border-l-2 border-amber-400' : 'border-l-2 border-transparent'}`}>
                    <input
                      type="checkbox"
                      disabled={!vid}
                      aria-label={`Select ${cv.title || vid || 'variant'}`}
                      checked={!!vid && importer.selected.has(vid)}
                      onChange={() => vid && importer.toggle(vid)}
                      className="rounded border-slate-600 bg-slate-900 text-cyan-500 focus:ring-cyan-500 disabled:opacity-30"
                    />
                  </td>
                  <td className="py-3 px-3 max-w-[300px]">
                    <div className="text-slate-100 font-medium break-words leading-snug">{cv.title || '—'}</div>
                    <div className="mt-0.5 flex flex-wrap items-center gap-2 text-[11px] font-mono">
                      {(cv.vcv_id || vid) && (
                        <NcbiLink
                          href={safeHref(cv.url) || (vid ? ncbiUrls.clinvarVariation(vid) : null)}
                          className="text-cyan-300 hover:text-cyan-200"
                        >
                          {cv.vcv_id || `Variation ${vid}`}
                        </NcbiLink>
                      )}
                      {imported && (
                        <span className="px-1.5 py-0.5 rounded bg-emerald-500/15 border border-emerald-500/30 text-emerald-300 font-sans text-[10px]">
                          in catalog
                        </span>
                      )}
                    </div>
                  </td>
                  <td className="py-3 px-3">
                    {genes.length === 0 ? (
                      <span className="text-slate-500">—</span>
                    ) : (
                      <div className="flex flex-wrap gap-1">
                        {genes.slice(0, 3).map((g, i) => (
                          <NcbiLink
                            key={`${g.symbol}-${i}`}
                            href={g.gene_id ? ncbiUrls.geneId(g.gene_id) : g.symbol ? ncbiUrls.geneSymbol(g.symbol) : null}
                            icon={false}
                            className="px-1.5 py-0.5 rounded-md bg-slate-800 border border-slate-700 font-mono font-bold text-white hover:border-cyan-500/60"
                          >
                            {g.symbol || g.gene_id}
                          </NcbiLink>
                        ))}
                        {genes.length > 3 && <span className="text-[10px] text-slate-500">+{genes.length - 3}</span>}
                      </div>
                    )}
                  </td>
                  <td className="py-3 px-3 font-mono whitespace-nowrap">
                    <div className="text-slate-200">{formatLocus(cv.chr, cv.start, cv.stop) || '—'}</div>
                    {(cv.ref || cv.alt) && (
                      <div className="text-[11px] text-slate-400">
                        {cv.ref || '?'} &rarr; <span className="text-rose-300">{cv.alt || '?'}</span>
                      </div>
                    )}
                  </td>
                  <td className="py-3 px-3 font-mono">
                    {cv.rsid ? <RsidLink rsid={cv.rsid} /> : <span className="text-slate-500">—</span>}
                  </td>
                  <td className="py-3 px-3">
                    <ClassificationBadge raw={cv.classification} normalized={cv.classification_normalized} size="xs" />
                  </td>
                  <td className="py-3 px-3">
                    <ReviewStars stars={cv.review_stars} reviewStatus={cv.review_status} />
                    {cv.last_evaluated && <div className="text-[10px] text-slate-500 mt-0.5">eval. {cv.last_evaluated}</div>}
                  </td>
                  <td className="py-3 px-3 max-w-[220px]">
                    {conditions.length === 0 ? (
                      <span className="text-slate-500">—</span>
                    ) : (
                      <div className="text-slate-300 leading-snug">
                        {conditions
                          .slice(0, 2)
                          .map((c) => c?.name)
                          .filter(Boolean)
                          .join('; ')}
                        {conditions.length > 2 && <span className="text-slate-500"> +{conditions.length - 2} more</span>}
                      </div>
                    )}
                  </td>
                  <td className="py-3 px-3">
                    <GenotypeCell yourGenotype={it.your_genotype} />
                  </td>
                  <td className="py-3 px-3">
                    <div className="flex flex-col items-end gap-1">
                      {onDetails && lookupParamsFor(cv) && (
                        <button
                          type="button"
                          onClick={() => onDetails(lookupParamsFor(cv))}
                          className="px-2 py-1 rounded-md bg-slate-800 hover:bg-slate-700 text-slate-200 text-[11px] whitespace-nowrap"
                        >
                          Details
                        </button>
                      )}
                      {imported && patch && onOpenInExplorer && (
                        <button
                          type="button"
                          onClick={() => onOpenInExplorer(patch)}
                          className="px-2 py-1 rounded-md bg-emerald-500/15 hover:bg-emerald-500/25 text-emerald-200 text-[11px] whitespace-nowrap"
                        >
                          View in Explorer
                        </button>
                      )}
                    </div>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ExampleChips({ examples, onPick }) {
  return (
    <div className="flex flex-wrap items-center gap-1.5 text-[11px]">
      <span className="text-slate-500">Try:</span>
      {examples.map((ex) => (
        <button
          key={ex}
          type="button"
          onClick={() => onPick(ex)}
          className="px-2 py-0.5 rounded-full border border-slate-800 bg-slate-950 text-slate-300 hover:text-cyan-300 hover:border-slate-700 font-mono"
        >
          {ex}
        </button>
      ))}
    </div>
  );
}

function StatTile({ label, value, tone = 'text-white', hint }) {
  return (
    <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
      <div className="text-[10px] uppercase tracking-wider text-slate-500">{label}</div>
      <div className={`text-xl font-bold font-mono mt-0.5 ${tone}`}>{value}</div>
      {hint && <div className="text-[10px] text-slate-500 mt-0.5">{hint}</div>}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------

export default function NcbiExplorer({ onOpenInExplorer, onCatalogChanged }) {
  const [mode, setMode] = useState('search'); // 'search' | 'gene' | 'lookup'

  // ClinVar search
  const [term, setTerm] = useState('');
  const [searchLimit, setSearchLimit] = useState(25);
  const [searchData, setSearchData] = useState(null);
  const [searchLoading, setSearchLoading] = useState(false);
  const [searchError, setSearchError] = useState(null);
  const searchImporter = useNcbiImport(onCatalogChanged);

  // Gene panel scan
  const [geneSymbol, setGeneSymbol] = useState('');
  const [significance, setSignificance] = useState('pathogenic');
  const [geneLimit, setGeneLimit] = useState(40);
  const [carriedOnly, setCarriedOnly] = useState(false);
  const [geneData, setGeneData] = useState(null);
  const [geneLoading, setGeneLoading] = useState(false);
  const [geneError, setGeneError] = useState(null);
  const geneImporter = useNcbiImport(onCatalogChanged);

  // Variant lookup
  const [lookupText, setLookupText] = useState('');
  const [includePubs, setIncludePubs] = useState(true);
  const [lookupData, setLookupData] = useState(null);
  const [lookupLoading, setLookupLoading] = useState(false);
  const [lookupError, setLookupError] = useState(null);

  // Latest request per mode; older responses arriving later are ignored.
  const searchSeq = useRef(0);
  const geneSeq = useRef(0);
  const lookupSeq = useRef(0);

  const runSearch = async (overrideTerm) => {
    const t = String(overrideTerm ?? term).trim();
    if (!t) {
      setSearchError('Enter a ClinVar search term.');
      return;
    }
    const seq = ++searchSeq.current;
    setTerm(t);
    setSearchLoading(true);
    setSearchError(null);
    searchImporter.reset();
    try {
      const body = await searchNcbi(t, searchLimit);
      if (seq !== searchSeq.current) return;
      setSearchData(body);
    } catch (err) {
      if (seq !== searchSeq.current) return;
      setSearchError(err.message || 'ClinVar search failed');
    } finally {
      if (seq === searchSeq.current) setSearchLoading(false);
    }
  };

  const runGeneScan = async (overrideSymbol) => {
    const sym = String(overrideSymbol ?? geneSymbol).trim().toUpperCase();
    if (!sym) {
      setGeneError('Enter a gene symbol.');
      return;
    }
    if (!GENE_SYMBOL_RE.test(sym)) {
      setGeneError(`"${sym}" does not look like a gene symbol.`);
      return;
    }
    const seq = ++geneSeq.current;
    setGeneSymbol(sym);
    setGeneLoading(true);
    setGeneError(null);
    geneImporter.reset();
    try {
      const body = await fetchNcbiGene(sym, { significance, limit: geneLimit });
      if (seq !== geneSeq.current) return;
      setGeneData(body);
    } catch (err) {
      if (seq !== geneSeq.current) return;
      setGeneError(err.message || 'Gene scan failed');
    } finally {
      if (seq === geneSeq.current) setGeneLoading(false);
    }
  };

  const runLookup = async (params, label) => {
    const seq = ++lookupSeq.current;
    setMode('lookup');
    if (label !== undefined) setLookupText(label);
    setLookupLoading(true);
    setLookupError(null);
    try {
      const body = await fetchNcbiVariant({ ...params, include_publications: includePubs });
      if (seq !== lookupSeq.current) return;
      setLookupData(body);
    } catch (err) {
      if (seq !== lookupSeq.current) return;
      setLookupData(null);
      setLookupError(err.message || 'Variant lookup failed');
    } finally {
      if (seq === lookupSeq.current) setLookupLoading(false);
    }
  };

  const submitLookup = (e, override) => {
    if (e) e.preventDefault();
    const text = override ?? lookupText;
    const parsed = parseLookupInput(text);
    if (parsed.error) {
      setLookupError(parsed.error);
      return;
    }
    runLookup(parsed.params, String(text).trim());
  };

  const openDetails = (params) => runLookup(params, lookupLabelFor(params));

  const geneItems = useMemo(() => {
    const items = asArray(geneData?.variants).filter((v) => v && v.clinvar);
    const ranked = items
      .map((v, i) => ({ v, i, c: carriesNcbiAlt(v.your_genotype) ? 1 : 0 }))
      .sort((a, b) => b.c - a.c || a.i - b.i)
      .map((x) => x.v);
    return carriedOnly ? ranked.filter((v) => carriesNcbiAlt(v.your_genotype)) : ranked;
  }, [geneData, carriedOnly]);

  const searchItems = useMemo(() => asArray(searchData?.items).filter((v) => v && v.clinvar), [searchData]);

  const modeBtn = (key, Icon, label) => (
    <button
      type="button"
      onClick={() => setMode(key)}
      className={`px-3 py-1.5 rounded-md text-xs font-medium transition-all flex items-center gap-1.5 ${
        mode === key ? 'bg-slate-700 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'
      }`}
    >
      <Icon className="w-3.5 h-3.5" />
      {label}
    </button>
  );

  return (
    <section className="rounded-2xl bg-slate-900/80 border border-slate-800/90 p-4 sm:p-5 space-y-5 shadow-xl">
      <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-3 border-b border-slate-800 pb-4">
        <div>
          <h2 className="text-base font-bold text-white flex items-center gap-2">
            <Globe className="w-4 h-4 text-cyan-400" />
            NCBI Live Explorer
          </h2>
          <p className="text-xs text-slate-400 mt-0.5 max-w-2xl">
            Queries{' '}
            <NcbiLink href="https://www.ncbi.nlm.nih.gov/clinvar/" className="text-cyan-300">
              ClinVar
            </NcbiLink>
            ,{' '}
            <NcbiLink href="https://www.ncbi.nlm.nih.gov/snp/" className="text-cyan-300">
              dbSNP
            </NcbiLink>
            ,{' '}
            <NcbiLink href="https://www.ncbi.nlm.nih.gov/gene/" className="text-cyan-300">
              Gene
            </NcbiLink>{' '}
            and{' '}
            <NcbiLink href="https://pubmed.ncbi.nlm.nih.gov/" className="text-cyan-300">
              PubMed
            </NcbiLink>{' '}
            live, then reads your own VCF at each GRCh38 locus to show whether you carry the ClinVar allele.
          </p>
        </div>
        <div className="inline-flex rounded-lg bg-slate-950 p-1 border border-slate-800 self-start">
          {modeBtn('search', Search, 'ClinVar search')}
          {modeBtn('gene', Dna, 'Gene panel scan')}
          {modeBtn('lookup', Crosshair, 'Variant lookup')}
        </div>
      </div>

      {/* ---------------- ClinVar search ---------------- */}
      {mode === 'search' && (
        <div className="space-y-4">
          <form
            onSubmit={(e) => {
              e.preventDefault();
              runSearch();
            }}
            className="flex flex-col sm:flex-row gap-2"
          >
            <div className="relative flex-1">
              <Search className="w-4 h-4 text-slate-500 absolute left-3 top-1/2 -translate-y-1/2" />
              <input
                type="text"
                value={term}
                onChange={(e) => setTerm(e.target.value)}
                maxLength={300}
                placeholder="RCV000019449, rs671, VCV000017865, HFE, hereditary hemochromatosis…"
                className={`${inputCls} pl-9`}
              />
            </div>
            <select value={searchLimit} onChange={(e) => setSearchLimit(Number(e.target.value))} className={selectCls} aria-label="Result limit">
              {[10, 25, 50, 100].map((n) => (
                <option key={n} value={n}>
                  {n} results
                </option>
              ))}
            </select>
            <button type="submit" disabled={searchLoading} className={primaryBtn}>
              {searchLoading ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : <Search className="w-3.5 h-3.5" />}
              Search ClinVar
            </button>
          </form>
          <ExampleChips examples={SEARCH_EXAMPLES} onPick={(ex) => runSearch(ex)} />

          <ErrorBanner message={searchError} onRetry={term ? () => runSearch() : undefined} />

          {searchLoading ? (
            <SkeletonCards count={3} />
          ) : searchData ? (
            <div className="space-y-3">
              <div className="flex flex-wrap items-center justify-between gap-2 text-xs text-slate-400">
                <span>
                  <span className="font-semibold text-white font-mono">{Number(searchData.total || 0).toLocaleString()}</span> ClinVar
                  records for <span className="font-mono text-slate-200">{searchData.term || term}</span>
                  {searchItems.length < Number(searchData.total || 0) && <> · showing {searchItems.length}</>}
                </span>
                {searchData.translated_query && (
                  <span className="font-mono text-[10px] text-slate-500 break-all" title="E-utilities query sent to ClinVar">
                    query: {searchData.translated_query}
                  </span>
                )}
              </div>
              <ErrorsNotice errors={searchData.errors} />
              {searchItems.length > 0 && (
                <ImportBar importer={searchImporter} items={searchItems} onOpenInExplorer={onOpenInExplorer} />
              )}
              <VariantResultsTable
                items={searchItems}
                importer={searchImporter}
                onDetails={openDetails}
                onOpenInExplorer={onOpenInExplorer}
                emptyTitle="No ClinVar records matched."
                emptyHint="Try an RCV/VCV accession, an rsID, a gene symbol or a condition name."
              />
              {searchData.fetched_at && (
                <p className="text-[10px] text-slate-500 text-right">Fetched {formatFetchedAt(searchData.fetched_at)}</p>
              )}
            </div>
          ) : (
            <EmptyState title="Search NCBI ClinVar live">
              Results show the ClinVar germline classification alongside your own genotype at each locus.
            </EmptyState>
          )}
        </div>
      )}

      {/* ---------------- Gene panel scan ---------------- */}
      {mode === 'gene' && (
        <div className="space-y-4">
          <form
            onSubmit={(e) => {
              e.preventDefault();
              runGeneScan();
            }}
            className="flex flex-col sm:flex-row gap-2"
          >
            <div className="relative flex-1">
              <Dna className="w-4 h-4 text-slate-500 absolute left-3 top-1/2 -translate-y-1/2" />
              <input
                type="text"
                value={geneSymbol}
                onChange={(e) => setGeneSymbol(e.target.value)}
                maxLength={40}
                placeholder="Gene symbol, e.g. HFE, ALDH2, BRCA1"
                className={`${inputCls} pl-9 font-mono uppercase`}
              />
            </div>
            <select value={significance} onChange={(e) => setSignificance(e.target.value)} className={selectCls} aria-label="Significance">
              {SIGNIFICANCE_OPTIONS.map((o) => (
                <option key={o.value} value={o.value}>
                  {o.label}
                </option>
              ))}
            </select>
            <select value={geneLimit} onChange={(e) => setGeneLimit(Number(e.target.value))} className={selectCls} aria-label="Variant limit">
              {[20, 40, 100].map((n) => (
                <option key={n} value={n}>
                  scan {n}
                </option>
              ))}
            </select>
            <button type="submit" disabled={geneLoading} className={primaryBtn}>
              {geneLoading ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : <DatabaseZap className="w-3.5 h-3.5" />}
              Scan my genome
            </button>
          </form>
          <ExampleChips examples={GENE_EXAMPLES} onPick={(ex) => runGeneScan(ex)} />

          <ErrorBanner message={geneError} onRetry={geneSymbol ? () => runGeneScan() : undefined} />

          {geneLoading ? (
            <SkeletonCards count={3} />
          ) : geneData ? (
            <div className="space-y-4">
              <ErrorsNotice errors={geneData.errors} />
              <div className="grid grid-cols-1 lg:grid-cols-[minmax(0,2fr)_minmax(0,1fr)] gap-4">
                {geneData.gene ? (
                  <GeneCard gene={geneData.gene} />
                ) : (
                  <EmptyState title="No NCBI Gene record returned">The scan below may still include ClinVar variants.</EmptyState>
                )}
                <div className="grid grid-cols-3 lg:grid-cols-1 gap-3">
                  <StatTile
                    label="ClinVar total"
                    value={Number(geneData.clinvar_total || 0).toLocaleString()}
                    hint={SIGNIFICANCE_OPTIONS.find((o) => o.value === (geneData.significance || significance))?.label}
                  />
                  <StatTile label="Variants scanned" value={asArray(geneData.variants).length.toLocaleString()} tone="text-cyan-300" hint="read in your VCF" />
                  <StatTile
                    label="You carry"
                    value={Number(geneData.carried_count || 0).toLocaleString()}
                    tone={Number(geneData.carried_count) > 0 ? 'text-amber-300' : 'text-emerald-300'}
                    hint="ClinVar ALT allele present"
                  />
                </div>
              </div>

              <div className="flex flex-wrap items-center justify-between gap-2">
                <label className="inline-flex items-center gap-2 px-3 py-1.5 rounded-xl bg-slate-950 border border-slate-700/80 cursor-pointer select-none text-xs text-slate-200">
                  <input
                    type="checkbox"
                    checked={carriedOnly}
                    onChange={(e) => setCarriedOnly(e.target.checked)}
                    className="rounded border-slate-600 bg-slate-900 text-cyan-500 focus:ring-cyan-500"
                  />
                  Carried only
                </label>
                <span className="text-[11px] text-slate-500">Carried rows are highlighted and listed first.</span>
              </div>

              {geneItems.length > 0 && <ImportBar importer={geneImporter} items={geneItems} onOpenInExplorer={onOpenInExplorer} />}
              <VariantResultsTable
                items={geneItems}
                importer={geneImporter}
                onDetails={openDetails}
                onOpenInExplorer={onOpenInExplorer}
                emptyTitle={carriedOnly ? 'You carry none of the scanned ClinVar alleles.' : 'No ClinVar variants for this filter.'}
                emptyHint={carriedOnly ? 'Untick "Carried only" to see every scanned variant.' : 'Try another significance filter.'}
              />
              {geneData.fetched_at && (
                <p className="text-[10px] text-slate-500 text-right">Fetched {formatFetchedAt(geneData.fetched_at)}</p>
              )}
            </div>
          ) : (
            <EmptyState title="Scan a gene against your genome">
              Pulls ClinVar variants for the gene from NCBI and reads your VCF at every GRCh38 locus.
            </EmptyState>
          )}
        </div>
      )}

      {/* ---------------- Variant lookup ---------------- */}
      {mode === 'lookup' && (
        <div className="space-y-4">
          <form onSubmit={submitLookup} className="flex flex-col sm:flex-row gap-2">
            <div className="relative flex-1">
              <Crosshair className="w-4 h-4 text-slate-500 absolute left-3 top-1/2 -translate-y-1/2" />
              <input
                type="text"
                value={lookupText}
                onChange={(e) => setLookupText(e.target.value)}
                maxLength={60}
                placeholder="rsID (rs671), RCV000019449, VCV000017865 or chr:position (12:111803962)"
                className={`${inputCls} pl-9 font-mono`}
              />
            </div>
            <label className="inline-flex items-center gap-2 px-3 py-2 rounded-xl bg-slate-950 border border-slate-700/80 cursor-pointer select-none text-xs text-slate-200">
              <input
                type="checkbox"
                checked={includePubs}
                onChange={(e) => setIncludePubs(e.target.checked)}
                className="rounded border-slate-600 bg-slate-900 text-cyan-500 focus:ring-cyan-500"
              />
              PubMed
            </label>
            <button type="submit" disabled={lookupLoading} className={primaryBtn}>
              {lookupLoading ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : <Crosshair className="w-3.5 h-3.5" />}
              Look up
            </button>
          </form>
          <ExampleChips
            examples={LOOKUP_EXAMPLES}
            onPick={(ex) => {
              setLookupText(ex);
              submitLookup(null, ex);
            }}
          />

          <ErrorBanner message={lookupError} />

          {lookupLoading ? (
            <SkeletonCards count={4} />
          ) : lookupData ? (
            <div className="max-w-4xl">
              {(() => {
                const patch = explorerPatchFor(asArray(lookupData.clinvar)[0]) ||
                  (lookupData.resolved?.rsid ? { variant_id: lookupData.resolved.rsid } : null);
                return patch && onOpenInExplorer ? (
                  <div className="mb-3 flex justify-end">
                    <button
                      type="button"
                      onClick={() => onOpenInExplorer(patch)}
                      className="px-3 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-xs flex items-center gap-1"
                    >
                      Search local catalog for {Object.values(patch)[0]}
                      <ArrowRight className="w-3.5 h-3.5" />
                    </button>
                  </div>
                ) : null;
              })()}
              <NcbiVariantDetail data={lookupData} showPublications={includePubs} />
            </div>
          ) : (
            <EmptyState title="Look up one variant across NCBI">
              ClinVar, dbSNP frequencies, NCBI Gene, PubMed and your own genotype in one view.
            </EmptyState>
          )}
        </div>
      )}
    </section>
  );
}
