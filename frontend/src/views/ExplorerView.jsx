import React, { useState, useEffect, useCallback, useRef } from 'react';
import {
  Dna,
  Search,
  Filter,
  RotateCcw,
  Download,
  ChevronLeft,
  ChevronRight,
  Activity,
  ShieldAlert,
  Pill,
  Layers,
  CheckCircle2,
  AlertTriangle,
  Info,
  Star,
  FileSpreadsheet,
  X,
  Terminal,
  ExternalLink,
  SlidersHorizontal,
  Sparkles,
  Database,
  Globe,
} from 'lucide-react';
import {
  fetchSummary,
  fetchFilterOptions,
  fetchVariants,
  fetchVariantDetail,
  fetchDirectVcfReads,
} from '../api';
import NcbiExplorer from '../components/NcbiExplorer';
import NcbiVariantPanel from '../components/NcbiVariantPanel';
import NcbiStatusPill from '../components/NcbiStatusPill';
import ClinvarStatusPill from '../components/ClinvarStatusPill';
import { GeneLinks, RcvLink, RcvList, RsidLink } from '../components/ncbiShared';

const INITIAL_FILTERS = {
  q: '',
  rcv_id: '',
  chr: '',
  position: '',
  variant_id: '',
  gene: '',
  risk_frequency: '',
  variant_type: '',
  your_status: '',
  condition: '',
  classification: '',
  evidence: '',
  source_file: '',
  carried_only: true,
  sort_by: 'priority',
  sort_dir: 'desc',
};

function getClassificationBadge(classification) {
  const c = (classification || '').toLowerCase();
  if (c.includes('conflicting')) {
    return 'bg-yellow-500/10 text-yellow-200 border-yellow-500/40';
  }
  // 'Pathogenic', 'Pathogenic / Risk Factor', 'Pathogenic/Likely Pathogenic', ...
  if (c.startsWith('pathogenic')) {
    return 'bg-rose-500/15 text-rose-300 border-rose-500/40';
  }
  if (c.includes('likely pathogenic')) {
    return 'bg-orange-500/15 text-orange-300 border-orange-500/40';
  }
  if (c.includes('risk factor') || c.includes('association')) {
    return 'bg-amber-500/15 text-amber-300 border-amber-500/40';
  }
  if (c.includes('drug response')) {
    return 'bg-purple-500/15 text-purple-300 border-purple-500/40';
  }
  if (c.includes('protective')) {
    return 'bg-teal-500/15 text-teal-300 border-teal-500/40';
  }
  if (c.includes('benign')) {
    return 'bg-emerald-500/15 text-emerald-300 border-emerald-500/40';
  }
  return 'bg-slate-700/50 text-slate-300 border-slate-600/60';
}

function getStatusBadge(status) {
  const s = (status || '').toLowerCase();
  if (s.includes('homozygous') || s.includes('gain')) {
    return 'bg-rose-500/20 text-rose-200 border-rose-500/40';
  }
  if (s.includes('carrier') || s.includes('heterozygous')) {
    return 'bg-amber-500/20 text-amber-200 border-amber-500/40';
  }
  if (s.includes('loss of heterozygosity') || s.includes('loh')) {
    return 'bg-fuchsia-500/20 text-fuchsia-200 border-fuchsia-500/40';
  }
  if (s.includes('reference') || s.includes('not detected')) {
    return 'bg-slate-800 text-slate-400 border-slate-700';
  }
  return 'bg-slate-700/40 text-slate-300 border-slate-600';
}

function getVariantTypeBadge(vtype) {
  const v = (vtype || '').toLowerCase();
  if (v.includes('cnv')) return 'bg-indigo-500/15 text-indigo-300 border-indigo-500/30';
  if (v.includes('loh')) return 'bg-fuchsia-500/15 text-fuchsia-300 border-fuchsia-500/30';
  if (v.includes('sv')) return 'bg-cyan-500/15 text-cyan-300 border-cyan-500/30';
  if (v.includes('del') || v.includes('ins') || v.includes('indel')) {
    return 'bg-sky-500/15 text-sky-300 border-sky-500/30';
  }
  return 'bg-slate-800/90 text-slate-300 border-slate-700';
}

function renderStars(count = 0) {
  const stars = Math.max(0, Math.min(4, Number(count) || 0));
  return (
    <div className="inline-flex items-center gap-0.5" title={`${stars} of 4 ClinVar review stars`}>
      {[1, 2, 3, 4].map((idx) => (
        <Star
          key={idx}
          className={`w-3.5 h-3.5 ${
            idx <= stars ? 'text-amber-400 fill-amber-400' : 'text-slate-700 fill-slate-800'
          }`}
        />
      ))}
    </div>
  );
}

function parseAllelicDepth(adStr) {
  if (!adStr || adStr === '.' || adStr.startsWith('PR:')) return null;
  const parts = String(adStr)
    .split(',')
    .map((n) => parseInt(n, 10))
    .filter((n) => !Number.isNaN(n));
  if (parts.length < 2) return null;
  const refCount = parts[0];
  const altCount = parts.slice(1).reduce((a, b) => a + b, 0);
  const total = refCount + altCount;
  if (total <= 0) return null;
  return {
    refCount,
    altCount,
    total,
    altPct: Math.round((altCount / total) * 100),
  };
}

/**
 * Clinical & variant explorer for ONE completed analysis. The analysis is selected in
 * the app shell; every API call carries its id (see api.js setActiveJob).
 */
export default function ExplorerView({ job }) {
  const sampleLabel = job?.label || 'this sample';
  const [summary, setSummary] = useState(null);
  const [filterOptions, setFilterOptions] = useState({
    chromosomes: [],
    variant_types: [],
    statuses: [],
    classifications: [],
    evidence_options: [],
    source_files: [],
    popular_genes: [],
    risk_frequency_presets: [],
  });

  const [filters, setFilters] = useState(INITIAL_FILTERS);
  const [activePreset, setActivePreset] = useState('carried');
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(25);

  const [results, setResults] = useState({
    total: 0,
    page: 1,
    page_size: 25,
    total_pages: 1,
    items: [],
  });
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  // Selected variant drawer state
  const [selectedVariant, setSelectedVariant] = useState(null);
  const [detailLoading, setDetailLoading] = useState(false);
  const [drawerTab, setDrawerTab] = useState('vcf'); // 'vcf' | 'ncbi'

  // Direct VCF Coordinate Reader state
  const [activeTab, setActiveTab] = useState('explorer'); // 'explorer' | 'raw_vcf' | 'ncbi'
  const [rawVcfQuery, setRawVcfQuery] = useState({
    chr: '12',
    position: '111803900-111804050',
    variants_only: false,
    limit: 50,
  });
  const [rawVcfData, setRawVcfData] = useState(null);
  const [rawVcfLoading, setRawVcfLoading] = useState(false);
  const [rawVcfError, setRawVcfError] = useState(null);

  // Load initial summary and filter metadata
  useEffect(() => {
    Promise.all([fetchSummary(), fetchFilterOptions()])
      .then(([sumData, optData]) => {
        setSummary(sumData);
        setFilterOptions(optData);
      })
      .catch((err) => {
        setError(err.message);
      });
  }, []);

  // Only the most recent catalog query may update the table (responses can arrive out of order).
  const variantsSeq = useRef(0);
  const loadVariants = useCallback(async () => {
    const seq = ++variantsSeq.current;
    setLoading(true);
    setError(null);
    try {
      const data = await fetchVariants({
        ...filters,
        page,
        page_size: pageSize,
      });
      if (seq !== variantsSeq.current) return;
      setResults(data);
    } catch (err) {
      if (seq !== variantsSeq.current) return;
      setError(err.message || 'Failed to query genetic variants');
    } finally {
      if (seq === variantsSeq.current) setLoading(false);
    }
  }, [filters, page, pageSize]);
  // Callbacks from long-running NCBI actions always reload with the *current* filters.
  const loadVariantsRef = useRef(loadVariants);
  useEffect(() => {
    loadVariantsRef.current = loadVariants;
  }, [loadVariants]);

  useEffect(() => {
    const timer = setTimeout(() => {
      loadVariants();
    }, 180);
    return () => clearTimeout(timer);
  }, [loadVariants]);

  const updateFilter = (field, value) => {
    setActivePreset('custom');
    setPage(1);
    setFilters((prev) => ({
      ...prev,
      [field]: value,
    }));
  };

  const applyPreset = (presetKey) => {
    setActivePreset(presetKey);
    setPage(1);
    if (presetKey === 'carried') {
      setFilters({ ...INITIAL_FILTERS, carried_only: true });
    } else if (presetKey === 'pathogenic') {
      setFilters({
        ...INITIAL_FILTERS,
        carried_only: true,
        classification: 'Pathogenic',
      });
    } else if (presetKey === 'risk_factor') {
      setFilters({
        ...INITIAL_FILTERS,
        carried_only: true,
        classification: 'Risk Factor',
      });
    } else if (presetKey === 'pharmacogenomics') {
      setFilters({
        ...INITIAL_FILTERS,
        carried_only: true,
        classification: 'Drug Response',
      });
    } else if (presetKey === 'structural_cnv') {
      setFilters({
        ...INITIAL_FILTERS,
        carried_only: true,
        variant_type: 'CNV',
      });
    } else if (presetKey === 'sv_all') {
      setFilters({
        ...INITIAL_FILTERS,
        carried_only: true,
        variant_type: 'SV',
      });
    } else if (presetKey === 'all_loci') {
      setFilters({
        ...INITIAL_FILTERS,
        carried_only: false,
      });
    }
  };

  const resetAllFilters = () => {
    setActivePreset('carried');
    setPage(1);
    setFilters(INITIAL_FILTERS);
  };

  const handleInspectVariant = async (item) => {
    setSelectedVariant(item);
    setDetailLoading(true);
    try {
      const fullDetail = await fetchVariantDetail(item.id);
      setSelectedVariant(fullDetail);
    } catch (_e) {
      // Keep row data if detail call fails
    } finally {
      setDetailLoading(false);
    }
  };

  const refreshSummary = useCallback(() => {
    fetchSummary()
      .then(setSummary)
      .catch(() => {});
  }, []);

  // Called by the NCBI explorer after importing records into the catalog
  const handleCatalogChanged = useCallback(() => {
    refreshSummary();
    loadVariantsRef.current();
  }, [refreshSummary]);

  // Jump from the NCBI explorer to the local catalog Explorer with a filter applied
  const openInExplorer = useCallback((patch = {}) => {
    setActivePreset('custom');
    setPage(1);
    setFilters({ ...INITIAL_FILTERS, carried_only: false, ...patch });
    setActiveTab('explorer');
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }, []);

  // Called by the drawer NCBI panel after "Sync catalog from NCBI"
  const handleNcbiSynced = useCallback(
    (body) => {
      const recordId = body?.record_id;
      refreshSummary();
      loadVariantsRef.current();
      if (recordId !== undefined && recordId !== null) {
        fetchVariantDetail(recordId)
          .then((fresh) => {
            setSelectedVariant((cur) => (cur && cur.id === recordId ? fresh : cur));
          })
          .catch(() => {});
      }
    },
    [refreshSummary]
  );

  const handleRunRawVcfRead = async (e) => {
    if (e) e.preventDefault();
    setRawVcfLoading(true);
    setRawVcfError(null);
    try {
      const data = await fetchDirectVcfReads(rawVcfQuery);
      setRawVcfData(data);
    } catch (err) {
      setRawVcfError(err.message);
    } finally {
      setRawVcfLoading(false);
    }
  };

  useEffect(() => {
    if (activeTab === 'raw_vcf' && !rawVcfData && !rawVcfLoading) {
      handleRunRawVcfRead();
    }
  }, [activeTab]);

  const exportCsv = () => {
    if (!results.items || results.items.length === 0) return;
    const headers = [
      'RCV ID',
      'Chromosome',
      'Position',
      'Variant ID',
      'Gene',
      'Your Data',
      'Risk Version',
      'Risk Frequency',
      'Variant Type',
      'Your Status',
      'Condition',
      'Classification',
      'Evidence',
      'Source VCF',
    ];
    const rows = results.items.map((r) => [
      r.rcv_id,
      r.chr,
      r.position,
      r.variant_id,
      r.gene,
      r.your_data,
      r.risk_version,
      r.risk_frequency_text,
      r.variant_type,
      r.your_status,
      r.condition,
      r.classification,
      r.evidence,
      r.source_file,
    ]);
    const csvContent =
      'data:text/csv;charset=utf-8,' +
      [headers, ...rows]
        .map((row) =>
          row.map((cell) => `"${String(cell ?? '').replace(/"/g, '""')}"`).join(',')
        )
        .join('\n');
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement('a');
    link.setAttribute('href', encodedUri);
    const safeLabel = String(sampleLabel).replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'analysis';
    link.setAttribute('download', `genetic_analysis_${safeLabel}_page${page}.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const activeFilterTags = [
    { key: 'q', label: 'Search' },
    { key: 'rcv_id', label: 'RCV ID' },
    { key: 'chr', label: 'Chr' },
    { key: 'position', label: 'Position' },
    { key: 'variant_id', label: 'Variant ID' },
    { key: 'gene', label: 'Gene' },
    { key: 'risk_frequency', label: 'Risk Freq' },
    { key: 'variant_type', label: 'Variant Type' },
    { key: 'your_status', label: 'Your Status' },
    { key: 'condition', label: 'Condition' },
    { key: 'classification', label: 'Classification' },
    { key: 'evidence', label: 'Evidence' },
    { key: 'source_file', label: 'VCF File' },
  ].filter((f) => filters[f.key] && String(filters[f.key]).trim() !== '');

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 flex flex-col">
      {/* Top Navigation & Sample Header */}
      <header className="sticky top-14 z-30 border-b border-slate-800/80 bg-slate-950/90 backdrop-blur-md">
        <div className="max-w-[1600px] mx-auto px-4 sm:px-6 py-3.5 flex flex-wrap items-center justify-between gap-4">
          <div className="flex items-center gap-3.5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center shadow-lg shadow-cyan-500/20">
              <Dna className="w-6 h-6 text-white" />
            </div>
            <div>
              <div className="flex items-center gap-2.5">
                <h1 className="text-lg font-bold tracking-tight text-white">
                  GenomeScope WGS Explorer
                </h1>
                <span className="px-2 py-0.5 text-xs font-mono font-semibold rounded-md bg-cyan-500/15 text-cyan-300 border border-cyan-500/30">
                  Sample: {summary?.sample_id || sampleLabel}
                </span>
                <span className="hidden sm:inline-block px-2 py-0.5 text-xs font-mono rounded-md bg-slate-800 text-slate-300 border border-slate-700">
                  {summary?.reference_genome || 'GRCh38'}
                </span>
              </div>
              <p className="text-xs text-slate-400">
                Direct VCF Variant &amp; Clinical Interpretation Engine • Analysis{' '}
                <code className="text-cyan-400 font-mono">{job?.id ? job.id.slice(0, 8) : '—'}</code>{' '}
                (isolated temp workspace)
              </p>
            </div>
          </div>

          {/* Mode Switcher & VCF File Indicators */}
          <div className="flex flex-wrap items-center gap-2.5">
            <ClinvarStatusPill />
            <NcbiStatusPill />
            <div className="inline-flex flex-wrap rounded-lg bg-slate-900 p-1 border border-slate-800">
              <button
                type="button"
                onClick={() => setActiveTab('explorer')}
                className={`px-3 py-1.5 rounded-md text-xs font-medium transition-all flex items-center gap-1.5 ${
                  activeTab === 'explorer'
                    ? 'bg-cyan-500 text-slate-950 font-semibold shadow-sm'
                    : 'text-slate-400 hover:text-slate-200'
                }`}
              >
                <Sparkles className="w-3.5 h-3.5" />
                Clinical &amp; Variant Explorer
              </button>
              <button
                type="button"
                onClick={() => setActiveTab('raw_vcf')}
                className={`px-3 py-1.5 rounded-md text-xs font-medium transition-all flex items-center gap-1.5 ${
                  activeTab === 'raw_vcf'
                    ? 'bg-cyan-500 text-slate-950 font-semibold shadow-sm'
                    : 'text-slate-400 hover:text-slate-200'
                }`}
              >
                <Terminal className="w-3.5 h-3.5" />
                Direct VCF Coordinate Streamer
              </button>
              <button
                type="button"
                onClick={() => setActiveTab('ncbi')}
                className={`px-3 py-1.5 rounded-md text-xs font-medium transition-all flex items-center gap-1.5 ${
                  activeTab === 'ncbi'
                    ? 'bg-cyan-500 text-slate-950 font-semibold shadow-sm'
                    : 'text-slate-400 hover:text-slate-200'
                }`}
              >
                <Globe className="w-3.5 h-3.5" />
                NCBI Live Explorer
              </button>
            </div>

            <button
              type="button"
              onClick={exportCsv}
              className="px-3 py-2 rounded-lg bg-slate-900 hover:bg-slate-800 text-slate-200 border border-slate-700 text-xs font-medium flex items-center gap-1.5 transition"
              title="Export current VCF results to CSV"
            >
              <Download className="w-3.5 h-3.5 text-cyan-400" />
              Export CSV
            </button>
          </div>
        </div>
      </header>

      {/* Main Content Container */}
      <main className="flex-1 max-w-[1600px] w-full mx-auto px-4 sm:px-6 py-5 space-y-5">
        {/* KPI Summary Cards */}
        <section className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-3">
          <button
            type="button"
            onClick={() => applyPreset('carried')}
            className={`text-left p-3.5 rounded-xl border transition-all ${
              activePreset === 'carried'
                ? 'bg-cyan-950/40 border-cyan-500/60 shadow-lg shadow-cyan-950/50'
                : 'bg-slate-900/70 border-slate-800/80 hover:border-slate-700'
            }`}
          >
            <div className="flex items-center justify-between text-xs text-slate-400 mb-1">
              <span>Carried Variants</span>
              <Dna className="w-4 h-4 text-cyan-400" />
            </div>
            <div className="text-2xl font-bold text-white font-mono">
              {summary ? summary.carried_variants.toLocaleString() : '—'}
            </div>
            <div className="text-[11px] text-cyan-400/90 mt-0.5">
              Detected in this sample's VCFs
            </div>
          </button>

          <button
            type="button"
            onClick={() => applyPreset('pathogenic')}
            className={`text-left p-3.5 rounded-xl border transition-all ${
              activePreset === 'pathogenic'
                ? 'bg-rose-950/40 border-rose-500/60 shadow-lg shadow-rose-950/50'
                : 'bg-slate-900/70 border-slate-800/80 hover:border-slate-700'
            }`}
          >
            <div className="flex items-center justify-between text-xs text-slate-400 mb-1">
              <span>Pathogenic / Likely</span>
              <ShieldAlert className="w-4 h-4 text-rose-400" />
            </div>
            <div className="text-2xl font-bold text-rose-300 font-mono">
              {summary ? summary.pathogenic_or_likely_carried.toLocaleString() : '—'}
            </div>
            <div className="text-[11px] text-rose-400/80 mt-0.5">
              High-priority clinical reads
            </div>
          </button>

          <button
            type="button"
            onClick={() => applyPreset('risk_factor')}
            className={`text-left p-3.5 rounded-xl border transition-all ${
              activePreset === 'risk_factor'
                ? 'bg-amber-950/40 border-amber-500/60 shadow-lg shadow-amber-950/50'
                : 'bg-slate-900/70 border-slate-800/80 hover:border-slate-700'
            }`}
          >
            <div className="flex items-center justify-between text-xs text-slate-400 mb-1">
              <span>Risk Factors</span>
              <AlertTriangle className="w-4 h-4 text-amber-400" />
            </div>
            <div className="text-2xl font-bold text-amber-300 font-mono">
              {summary ? summary.risk_factor_carried.toLocaleString() : '—'}
            </div>
            <div className="text-[11px] text-amber-400/80 mt-0.5">
              Disease susceptibility alleles
            </div>
          </button>

          <button
            type="button"
            onClick={() => applyPreset('pharmacogenomics')}
            className={`text-left p-3.5 rounded-xl border transition-all ${
              activePreset === 'pharmacogenomics'
                ? 'bg-purple-950/40 border-purple-500/60 shadow-lg shadow-purple-950/50'
                : 'bg-slate-900/70 border-slate-800/80 hover:border-slate-700'
            }`}
          >
            <div className="flex items-center justify-between text-xs text-slate-400 mb-1">
              <span>Drug Response</span>
              <Pill className="w-4 h-4 text-purple-400" />
            </div>
            <div className="text-2xl font-bold text-purple-300 font-mono">
              {summary ? summary.drug_response_carried.toLocaleString() : '—'}
            </div>
            <div className="text-[11px] text-purple-400/80 mt-0.5">
              Pharmacogenomic alleles
            </div>
          </button>

          <button
            type="button"
            onClick={() => applyPreset('structural_cnv')}
            className={`text-left p-3.5 rounded-xl border transition-all ${
              activePreset === 'structural_cnv'
                ? 'bg-indigo-950/40 border-indigo-500/60 shadow-lg shadow-indigo-950/50'
                : 'bg-slate-900/70 border-slate-800/80 hover:border-slate-700'
            }`}
          >
            <div className="flex items-center justify-between text-xs text-slate-400 mb-1">
              <span>CNV &amp; Structural</span>
              <Layers className="w-4 h-4 text-indigo-400" />
            </div>
            <div className="text-2xl font-bold text-indigo-300 font-mono">
              {summary ? summary.structural_cnv_carried.toLocaleString() : '—'}
            </div>
            <div className="text-[11px] text-indigo-400/80 mt-0.5">
              PASS Canvas CNV &amp; Manta SV
            </div>
          </button>

          <button
            type="button"
            onClick={() => applyPreset('all_loci')}
            className={`text-left p-3.5 rounded-xl border transition-all ${
              activePreset === 'all_loci'
                ? 'bg-emerald-950/40 border-emerald-500/60 shadow-lg shadow-emerald-950/50'
                : 'bg-slate-900/70 border-slate-800/80 hover:border-slate-700'
            }`}
          >
            <div className="flex items-center justify-between text-xs text-slate-400 mb-1">
              <span>Total Indexed Loci</span>
              <Database className="w-4 h-4 text-emerald-400" />
            </div>
            <div className="text-2xl font-bold text-emerald-300 font-mono">
              {summary ? summary.total_indexed_loci.toLocaleString() : '—'}
            </div>
            <div className="text-[11px] text-emerald-400/80 mt-0.5">
              ClinVar + Reference blocks
            </div>
          </button>
        </section>

        {activeTab === 'explorer' ? (
          <>
            {/* 11-Parameter Search & Filter Matrix */}
            <section className="rounded-2xl bg-slate-900/80 border border-slate-800/90 p-4 sm:p-5 shadow-xl space-y-4">
              {/* Top Row: Global Search, VCF Source Filter, Carried Toggle, Reset */}
              <div className="flex flex-col lg:flex-row items-stretch lg:items-center justify-between gap-3 pb-3 border-b border-slate-800/80">
                <div className="relative flex-1">
                  <Search className="w-4 h-4 text-slate-400 absolute left-3.5 top-1/2 -translate-y-1/2" />
                  <input
                    type="text"
                    value={filters.q}
                    onChange={(e) => updateFilter('q', e.target.value)}
                    placeholder="Quick search across RCV ID (e.g. RCV000782431), Variant ID (rs671), Gene (CYP2C19, ALDH2, BRCA1), Chr:Pos (12:111803962), Condition..."
                    className="w-full pl-10 pr-9 py-2.5 bg-slate-950 border border-slate-700/80 rounded-xl text-sm text-slate-100 placeholder-slate-500 focus:outline-none focus:border-cyan-500 focus:ring-1 focus:ring-cyan-500"
                  />
                  {filters.q && (
                    <button
                      type="button"
                      onClick={() => updateFilter('q', '')}
                      className="absolute right-3 top-1/2 -translate-y-1/2 text-slate-400 hover:text-slate-200"
                    >
                      <X className="w-4 h-4" />
                    </button>
                  )}
                </div>

                <div className="flex flex-wrap items-center gap-2.5">
                  {/* VCF File Selector */}
                  <select
                    value={filters.source_file}
                    onChange={(e) => updateFilter('source_file', e.target.value)}
                    className="px-3 py-2 bg-slate-950 border border-slate-700/80 rounded-xl text-xs text-slate-200 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="">All uploaded VCF files</option>
                    <option value="snp-indel.genome.vcf">
                      snp-indel.genome.vcf (SNP / indel gVCF)
                    </option>
                    <option value="cnv.vcf">cnv.vcf (Copy Number Variants)</option>
                    <option value="sv.vcf">sv.vcf (Structural Variants)</option>
                  </select>

                  {/* Carried Only Toggle */}
                  <label className="inline-flex items-center gap-2 px-3 py-2 rounded-xl bg-slate-950 border border-slate-700/80 cursor-pointer select-none text-xs font-medium text-slate-200">
                    <input
                      type="checkbox"
                      checked={filters.carried_only}
                      onChange={(e) => updateFilter('carried_only', e.target.checked)}
                      className="rounded border-slate-600 bg-slate-900 text-cyan-500 focus:ring-cyan-500"
                    />
                    <span>Carried / Detected in My Genome Only</span>
                  </label>

                  <button
                    type="button"
                    onClick={resetAllFilters}
                    className="px-3 py-2 rounded-xl bg-slate-800 hover:bg-slate-700 text-slate-200 text-xs font-medium flex items-center gap-1.5 transition"
                  >
                    <RotateCcw className="w-3.5 h-3.5" />
                    Reset Filters
                  </button>
                </div>
              </div>

              {/* All 11 Dedicated Parameter Inputs */}
              <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-3">
                {/* 1. RCV ID */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    1. RCV ID
                  </label>
                  <input
                    type="text"
                    value={filters.rcv_id}
                    onChange={(e) => updateFilter('rcv_id', e.target.value)}
                    placeholder="e.g. RCV000782431"
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 placeholder-slate-600 font-mono focus:outline-none focus:border-cyan-500"
                  />
                </div>

                {/* 2. Chromosome (chr) */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    2. Chromosome (chr)
                  </label>
                  <select
                    value={filters.chr}
                    onChange={(e) => updateFilter('chr', e.target.value)}
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="">All Chromosomes (1-22, X, Y, MT)</option>
                    {filterOptions.chromosomes.map((c) => (
                      <option key={c} value={c}>
                        Chromosome {c} (chr{c})
                      </option>
                    ))}
                  </select>
                </div>

                {/* 3. Position */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    3. Position (GRCh38)
                  </label>
                  <input
                    type="text"
                    value={filters.position}
                    onChange={(e) => updateFilter('position', e.target.value)}
                    placeholder="e.g. 111803962 or start-end"
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 placeholder-slate-600 font-mono focus:outline-none focus:border-cyan-500"
                  />
                </div>

                {/* 4. Variant ID */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    4. Variant ID
                  </label>
                  <input
                    type="text"
                    value={filters.variant_id}
                    onChange={(e) => updateFilter('variant_id', e.target.value)}
                    placeholder="e.g. rs671, rs4244285, Canvas:..."
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 placeholder-slate-600 font-mono focus:outline-none focus:border-cyan-500"
                  />
                </div>

                {/* 5. Gene */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    5. Gene Symbol
                  </label>
                  <input
                    type="text"
                    list="gene-suggestions"
                    value={filters.gene}
                    onChange={(e) => updateFilter('gene', e.target.value)}
                    placeholder="e.g. ALDH2, CYP2C19, BRCA1"
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 placeholder-slate-600 font-mono uppercase focus:outline-none focus:border-cyan-500"
                  />
                  <datalist id="gene-suggestions">
                    {filterOptions.popular_genes.map((g) => (
                      <option key={g} value={g} />
                    ))}
                  </datalist>
                </div>

                {/* 6. Risk Frequency */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    6. Risk Frequency
                  </label>
                  <input
                    type="text"
                    list="freq-presets"
                    value={filters.risk_frequency}
                    onChange={(e) => updateFilter('risk_frequency', e.target.value)}
                    placeholder="e.g. <0.01, <0.05, 0.01-0.10"
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 placeholder-slate-600 font-mono focus:outline-none focus:border-cyan-500"
                  />
                  <datalist id="freq-presets">
                    <option value="<0.001">Very Rare (&lt; 0.1%)</option>
                    <option value="<0.01">Rare (&lt; 1%)</option>
                    <option value="<0.05">Low Frequency (&lt; 5%)</option>
                    <option value="0.01-0.10">Moderate (1% - 10%)</option>
                    <option value=">=0.05">Common (&gt;= 5%)</option>
                  </datalist>
                </div>

                {/* 7. Variant Type */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    7. Variant Type
                  </label>
                  <select
                    value={filters.variant_type}
                    onChange={(e) => updateFilter('variant_type', e.target.value)}
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="">All Variant Types</option>
                    {filterOptions.variant_types.map((vt) => (
                      <option key={vt} value={vt}>
                        {vt}
                      </option>
                    ))}
                  </select>
                </div>

                {/* 8. Your Status */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    8. Your Status (VCF Call)
                  </label>
                  <select
                    value={filters.your_status}
                    onChange={(e) => {
                      const val = e.target.value;
                      setActivePreset('custom');
                      setPage(1);
                      setFilters((prev) => ({
                        ...prev,
                        your_status: val,
                        carried_only:
                          val.toLowerCase().includes('reference') ||
                          val.toLowerCase().includes('low quality')
                            ? false
                            : prev.carried_only,
                      }));
                    }}
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="">All Zygosity / Carrier Statuses</option>
                    {filterOptions.statuses.map((st) => (
                      <option key={st} value={st}>
                        {st}
                      </option>
                    ))}
                  </select>
                </div>

                {/* 9. Condition */}
                <div className="sm:col-span-2 lg:col-span-1">
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    9. Condition / Disease
                  </label>
                  <input
                    type="text"
                    value={filters.condition}
                    onChange={(e) => updateFilter('condition', e.target.value)}
                    placeholder="e.g. Clopidogrel, Cancer, Diabetes"
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 placeholder-slate-600 focus:outline-none focus:border-cyan-500"
                  />
                </div>

                {/* 10. Classification */}
                <div>
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    10. Classification
                  </label>
                  <select
                    value={filters.classification}
                    onChange={(e) => updateFilter('classification', e.target.value)}
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="">All Clinical Classifications</option>
                    {filterOptions.classifications.map((cl) => (
                      <option key={cl} value={cl}>
                        {cl}
                      </option>
                    ))}
                  </select>
                </div>

                {/* 11. Evidence */}
                <div className="sm:col-span-2 lg:col-span-2">
                  <label className="block text-[11px] font-semibold uppercase tracking-wider text-slate-400 mb-1">
                    11. Evidence (ClinVar Review Stars / VCF Quality)
                  </label>
                  <select
                    value={filters.evidence}
                    onChange={(e) => updateFilter('evidence', e.target.value)}
                    className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-lg text-xs text-slate-100 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="">All Evidence Levels &amp; Read Qualities</option>
                    <option value="Practice Guideline (4 Stars)">
                      ★★★★ Practice Guideline (4 Stars)
                    </option>
                    <option value="Expert Panel Reviewed (3 Stars)">
                      ★★★☆ Expert Panel Reviewed (3 Stars)
                    </option>
                    <option value="Strong - Multiple Submitters (2 Stars)">
                      ★★☆☆ Strong - Multiple Submitters (2 Stars)
                    </option>
                    <option value="Moderate - Single Submitter (1 Star)">
                      ★☆☆☆ Moderate - Single Submitter (1 Star)
                    </option>
                    <option value="Conflicting Interpretations">
                      ★☆☆☆ Conflicting Interpretations
                    </option>
                    <option value="WGS Canvas">WGS Canvas CNV Read Evidence</option>
                    <option value="WGS Manta">WGS Manta SV Read Evidence</option>
                    <option value="VCF PASS">VCF Filter = PASS</option>
                  </select>
                </div>
              </div>

              {/* Quick Gene & Example Search Chips */}
              <div className="flex flex-wrap items-center justify-between gap-2 pt-2 border-t border-slate-800/60 text-xs">
                <div className="flex flex-wrap items-center gap-1.5">
                  <span className="text-slate-400 font-medium mr-1">
                    Top carried genes in this sample:
                  </span>
                  {(summary?.top_carried_genes || []).length === 0 && (
                    <span className="text-slate-500 italic">none yet</span>
                  )}
                  {(summary?.top_carried_genes || []).slice(0, 10).map((g) => ({
                    gene: g.gene,
                    note: `${g.carried_count} carried`,
                  })).map((item) => (
                    <button
                      key={item.gene}
                      type="button"
                      onClick={() => {
                        setActivePreset('custom');
                        setPage(1);
                        setFilters((prev) => ({
                          ...prev,
                          gene: item.gene,
                          carried_only: false,
                        }));
                      }}
                      className={`px-2.5 py-1 rounded-full border text-[11px] font-mono transition ${
                        filters.gene.toUpperCase() === item.gene
                          ? 'bg-cyan-500/20 border-cyan-500 text-cyan-200 font-semibold'
                          : 'bg-slate-950/90 border-slate-800 text-slate-300 hover:border-slate-700 hover:text-cyan-300'
                      }`}
                    >
                      {item.gene}{' '}
                      <span className="text-slate-500 font-sans">({item.note})</span>
                    </button>
                  ))}
                </div>

                {/* Active filter chips */}
                {activeFilterTags.length > 0 && (
                  <div className="flex flex-wrap items-center gap-1.5">
                    {activeFilterTags.map((tag) => (
                      <span
                        key={tag.key}
                        className="inline-flex items-center gap-1 px-2 py-0.5 rounded-md bg-cyan-500/15 border border-cyan-500/30 text-cyan-200 text-[11px]"
                      >
                        <span className="text-cyan-400/80">{tag.label}:</span>
                        <span className="font-mono font-semibold">{filters[tag.key]}</span>
                        <button
                          type="button"
                          onClick={() => updateFilter(tag.key, '')}
                          className="hover:text-white"
                        >
                          <X className="w-3 h-3" />
                        </button>
                      </span>
                    ))}
                  </div>
                )}
              </div>
            </section>

            {/* Results Header Bar */}
            <div className="flex flex-wrap items-center justify-between gap-3 px-1">
              <div className="flex items-center gap-2 text-sm text-slate-300">
                <span className="font-semibold text-white font-mono">
                  {results.total.toLocaleString()}
                </span>
                <span>matching VCF variant reads</span>
                {filters.carried_only && (
                  <span className="px-2 py-0.5 text-xs rounded-full bg-cyan-500/15 text-cyan-300 border border-cyan-500/30">
                    Carried in {summary?.sample_id || sampleLabel}
                  </span>
                )}
              </div>

              <div className="flex items-center gap-3 text-xs">
                <div className="flex items-center gap-1.5">
                  <span className="text-slate-400">Sort by:</span>
                  <select
                    value={filters.sort_by}
                    onChange={(e) => updateFilter('sort_by', e.target.value)}
                    className="px-2.5 py-1.5 bg-slate-900 border border-slate-800 rounded-lg text-slate-200 focus:outline-none focus:border-cyan-500"
                  >
                    <option value="priority">Clinical Priority &amp; Evidence</option>
                    <option value="evidence">Evidence Rating (Stars)</option>
                    <option value="risk_frequency">Risk Allele Frequency</option>
                    <option value="position">Genomic Coordinate (Chr:Pos)</option>
                    <option value="gene">Gene Symbol (A-Z)</option>
                    <option value="classification">Classification</option>
                    <option value="variant_id">Variant ID</option>
                    <option value="rcv_id">RCV ID</option>
                  </select>
                </div>

                <div className="flex items-center gap-1.5">
                  <span className="text-slate-400">Rows:</span>
                  <select
                    value={pageSize}
                    onChange={(e) => {
                      setPageSize(Number(e.target.value));
                      setPage(1);
                    }}
                    className="px-2 py-1.5 bg-slate-900 border border-slate-800 rounded-lg text-slate-200 focus:outline-none focus:border-cyan-500"
                  >
                    <option value={15}>15</option>
                    <option value={25}>25</option>
                    <option value={50}>50</option>
                    <option value={100}>100</option>
                  </select>
                </div>
              </div>
            </div>

            {/* Error Banner */}
            {error && (
              <div className="p-4 rounded-xl bg-rose-950/50 border border-rose-500/40 text-rose-200 text-sm flex items-center justify-between">
                <span>{error}</span>
                <button
                  type="button"
                  onClick={loadVariants}
                  className="px-3 py-1 rounded-lg bg-rose-500/20 hover:bg-rose-500/30 text-xs font-semibold"
                >
                  Retry
                </button>
              </div>
            )}

            {/* Main Variant Reads Table */}
            <section className="rounded-2xl bg-slate-900/80 border border-slate-800/90 overflow-hidden shadow-2xl">
              <div className="overflow-x-auto">
                <table className="w-full text-left border-collapse">
                  <thead>
                    <tr className="border-b border-slate-800 bg-slate-950/80 text-[11px] font-semibold uppercase tracking-wider text-slate-400">
                      <th className="py-3.5 px-4">Variant ID &amp; RCV</th>
                      <th className="py-3.5 px-3">Chr : Position</th>
                      <th className="py-3.5 px-3">Gene</th>
                      <th className="py-3.5 px-4">Your Data (VCF Read)</th>
                      <th className="py-3.5 px-3">Risk Version &amp; Freq</th>
                      <th className="py-3.5 px-3">Variant Type</th>
                      <th className="py-3.5 px-3">Your Status</th>
                      <th className="py-3.5 px-4">Condition (Genetic Issue)</th>
                      <th className="py-3.5 px-3">Classification</th>
                      <th className="py-3.5 px-4">Evidence</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-800/60 text-xs">
                    {loading ? (
                      Array.from({ length: 8 }).map((_, i) => (
                        <tr key={i} className="animate-pulse">
                          <td colSpan={10} className="py-4 px-4">
                            <div className="h-5 bg-slate-800/70 rounded w-full" />
                          </td>
                        </tr>
                      ))
                    ) : results.items.length === 0 ? (
                      <tr>
                        <td colSpan={10} className="py-12 text-center text-slate-400">
                          <div className="max-w-md mx-auto space-y-2">
                            <p className="text-sm font-semibold text-slate-200">
                              No matching VCF variant reads found with current filters.
                            </p>
                            <p className="text-xs text-slate-400">
                              Try unchecking &ldquo;Carried / Detected in My Genome Only&rdquo; to view homozygous reference blocks (0/0), or reset filters.
                            </p>
                            <div className="pt-2 flex justify-center gap-2">
                              {filters.carried_only && (
                                <button
                                  type="button"
                                  onClick={() => updateFilter('carried_only', false)}
                                  className="px-3 py-1.5 rounded-lg bg-cyan-500/20 text-cyan-300 border border-cyan-500/30 text-xs font-medium"
                                >
                                  Include Reference Reads (0/0)
                                </button>
                              )}
                              <button
                                type="button"
                                onClick={resetAllFilters}
                                className="px-3 py-1.5 rounded-lg bg-slate-800 text-slate-200 text-xs font-medium"
                              >
                                Reset All Filters
                              </button>
                            </div>
                          </div>
                        </td>
                      </tr>
                    ) : (
                      results.items.map((row) => {
                        const adParsed = parseAllelicDepth(row.ad);
                        return (
                          <tr
                            key={row.id}
                            onClick={() => handleInspectVariant(row)}
                            className="hover:bg-slate-800/50 cursor-pointer transition-colors group"
                          >
                            {/* Variant ID & RCV ID */}
                            <td className="py-3.5 px-4 align-top">
                              <div className="font-mono font-semibold text-cyan-300 group-hover:text-cyan-200 flex items-center gap-1.5">
                                <RsidLink
                                  rsid={row.variant_id}
                                  className="text-cyan-300 group-hover:text-cyan-200 hover:text-white"
                                />
                              </div>
                              <div className="text-[11px] font-mono text-slate-400 mt-0.5">
                                <RcvLink rcv={row.rcv_id} className="text-slate-400 hover:text-cyan-200" />
                              </div>
                            </td>

                            {/* Chr : Position */}
                            <td className="py-3.5 px-3 align-top font-mono whitespace-nowrap">
                              <span className="px-1.5 py-0.5 rounded bg-slate-800 text-slate-200 font-semibold">
                                chr{row.chr}
                              </span>
                              <span className="text-slate-300 ml-1">
                                {Number(row.position).toLocaleString()}
                              </span>
                              {row.end_position && row.end_position > row.position && (
                                <div className="text-[10px] text-slate-500 mt-0.5">
                                  to {Number(row.end_position).toLocaleString()}
                                </div>
                              )}
                            </td>

                            {/* Gene */}
                            <td className="py-3.5 px-3 align-top">
                              <GeneLinks gene={row.gene} />
                              {row.hgvs && (
                                <div
                                  className="text-[10px] font-mono text-slate-400 mt-1 max-w-[150px] truncate"
                                  title={row.hgvs}
                                >
                                  {row.hgvs}
                                </div>
                              )}
                            </td>

                            {/* Your Data (VCF Read) */}
                            <td className="py-3.5 px-4 align-top">
                              <div className="font-mono font-semibold text-slate-100">
                                {row.alleles_called}{' '}
                                <span className="text-cyan-400">({row.genotype})</span>
                              </div>
                              <div className="text-[11px] font-mono text-slate-400 mt-0.5">
                                {row.dp ? `DP: ${row.dp}` : ''}{' '}
                                {row.gq ? `• GQ: ${row.gq}` : ''}{' '}
                                {row.ad && row.ad !== '.' ? `• AD: ${row.ad}` : ''}
                              </div>
                              {adParsed && (
                                <div
                                  className="w-28 h-1.5 bg-slate-800 rounded-full overflow-hidden flex mt-1.5"
                                  title={`Ref reads: ${adParsed.refCount}, Alt reads: ${adParsed.altCount} (${adParsed.altPct}% Alt)`}
                                >
                                  <div
                                    className="bg-slate-500 h-full"
                                    style={{ width: `${100 - adParsed.altPct}%` }}
                                  />
                                  <div
                                    className="bg-cyan-400 h-full"
                                    style={{ width: `${adParsed.altPct}%` }}
                                  />
                                </div>
                              )}
                            </td>

                            {/* Risk Version & Risk Frequency */}
                            <td className="py-3.5 px-3 align-top">
                              <div className="flex items-center gap-1.5 font-mono">
                                <span
                                  className="px-1.5 py-0.5 rounded bg-rose-500/15 border border-rose-500/30 text-rose-300 font-bold"
                                  title="Risk / Alternate Allele"
                                >
                                  {row.risk_version}
                                </span>
                                <span className="text-slate-500 text-[11px]">
                                  (Ref: {row.ref_allele})
                                </span>
                              </div>
                              <div className="text-[11px] font-mono text-slate-400 mt-1">
                                {row.risk_frequency_text}
                              </div>
                            </td>

                            {/* Variant Type */}
                            <td className="py-3.5 px-3 align-top">
                              <span
                                className={`inline-block px-2 py-0.5 rounded-md border text-[11px] font-medium ${getVariantTypeBadge(
                                  row.variant_type
                                )}`}
                              >
                                {row.variant_type}
                              </span>
                            </td>

                            {/* Your Status */}
                            <td className="py-3.5 px-3 align-top">
                              <span
                                className={`inline-block px-2.5 py-1 rounded-lg border text-[11px] font-semibold ${getStatusBadge(
                                  row.your_status
                                )}`}
                              >
                                {row.your_status}
                              </span>
                            </td>

                            {/* Condition */}
                            <td className="py-3.5 px-4 align-top max-w-[280px]">
                              <div className="text-slate-200 font-medium line-clamp-2 leading-snug">
                                {row.condition}
                              </div>
                            </td>

                            {/* Classification */}
                            <td className="py-3.5 px-3 align-top">
                              <span
                                className={`inline-block px-2.5 py-1 rounded-lg border text-[11px] font-semibold ${getClassificationBadge(
                                  row.classification
                                )}`}
                              >
                                {row.classification}
                              </span>
                            </td>

                            {/* Evidence */}
                            <td className="py-3.5 px-4 align-top">
                              <div className="flex items-center gap-1.5">
                                {renderStars(row.evidence_stars)}
                                <span
                                  className={`px-1.5 py-0.5 rounded text-[10px] font-mono ${
                                    row.filter_status === 'PASS'
                                      ? 'bg-emerald-500/15 text-emerald-300'
                                      : 'bg-amber-500/15 text-amber-300'
                                  }`}
                                >
                                  {row.filter_status}
                                </span>
                              </div>
                              <div className="text-[11px] text-slate-400 mt-1 line-clamp-1">
                                {row.evidence.split('|')[0]}
                              </div>
                            </td>
                          </tr>
                        );
                      })
                    )}
                  </tbody>
                </table>
              </div>

              {/* Pagination Footer */}
              <div className="px-4 py-3 bg-slate-950/90 border-t border-slate-800 flex flex-wrap items-center justify-between gap-3 text-xs text-slate-400">
                <div>
                  Showing page{' '}
                  <span className="font-semibold text-white">{results.page}</span> of{' '}
                  <span className="font-semibold text-white">{results.total_pages}</span>{' '}
                  ({results.total.toLocaleString()} total records)
                </div>
                <div className="flex items-center gap-2">
                  <button
                    type="button"
                    disabled={page <= 1}
                    onClick={() => setPage((p) => Math.max(1, p - 1))}
                    className="px-3 py-1.5 rounded-lg bg-slate-900 border border-slate-800 text-slate-200 disabled:opacity-40 hover:bg-slate-800 flex items-center gap-1"
                  >
                    <ChevronLeft className="w-3.5 h-3.5" />
                    Previous
                  </button>
                  <button
                    type="button"
                    disabled={page >= results.total_pages}
                    onClick={() => setPage((p) => Math.min(results.total_pages, p + 1))}
                    className="px-3 py-1.5 rounded-lg bg-slate-900 border border-slate-800 text-slate-200 disabled:opacity-40 hover:bg-slate-800 flex items-center gap-1"
                  >
                    Next
                    <ChevronRight className="w-3.5 h-3.5" />
                  </button>
                </div>
              </div>
            </section>
          </>
        ) : activeTab === 'raw_vcf' ? (
          /* Direct Raw VCF Coordinate Streamer Mode */
          <section className="rounded-2xl bg-slate-900/80 border border-slate-800/90 p-5 space-y-5 shadow-xl">
            <div className="flex flex-col md:flex-row md:items-center justify-between gap-3 border-b border-slate-800 pb-4">
              <div>
                <h2 className="text-base font-bold text-white flex items-center gap-2">
                  <Terminal className="w-4 h-4 text-cyan-400" />
                  Direct VCF Coordinate &amp; Region Reader
                </h2>
                <p className="text-xs text-slate-400 mt-0.5">
                  Binary-searches the VCF files uploaded for this analysis (
                  <code>snp-indel.genome.vcf</code>, <code>sv.vcf</code> and{' '}
                  <code>cnv.vcf</code>) using a sparse byte-offset index.
                </p>
              </div>
            </div>

            <form
              onSubmit={handleRunRawVcfRead}
              className="grid grid-cols-1 sm:grid-cols-4 gap-3 items-end"
            >
              <div>
                <label className="block text-xs text-slate-400 mb-1">
                  Chromosome (chr)
                </label>
                <select
                  value={rawVcfQuery.chr}
                  onChange={(e) =>
                    setRawVcfQuery((p) => ({ ...p, chr: e.target.value }))
                  }
                  className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-xl text-xs text-white"
                >
                  {filterOptions.chromosomes.map((c) => (
                    <option key={c} value={c}>
                      Chromosome {c}
                    </option>
                  ))}
                </select>
              </div>
              <div>
                <label className="block text-xs text-slate-400 mb-1">
                  Genomic Position or Range (GRCh38)
                </label>
                <input
                  type="text"
                  value={rawVcfQuery.position}
                  onChange={(e) =>
                    setRawVcfQuery((p) => ({ ...p, position: e.target.value }))
                  }
                  placeholder="e.g. 111803962 or 111803000-111805000"
                  className="w-full px-3 py-2 bg-slate-950 border border-slate-800 rounded-xl text-xs text-white font-mono"
                />
              </div>
              <div className="flex items-center gap-3">
                <label className="inline-flex items-center gap-2 text-xs text-slate-300 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={rawVcfQuery.variants_only}
                    onChange={(e) =>
                      setRawVcfQuery((p) => ({
                        ...p,
                        variants_only: e.target.checked,
                      }))
                    }
                    className="rounded border-slate-700 bg-slate-950 text-cyan-500"
                  />
                  <span>Skip Hom-Ref Blocks (ALT = .)</span>
                </label>
              </div>
              <div>
                <button
                  type="submit"
                  className="w-full py-2 px-4 rounded-xl bg-cyan-500 hover:bg-cyan-400 text-slate-950 font-semibold text-xs transition"
                >
                  Read VCF Locus
                </button>
              </div>
            </form>

            {rawVcfError && (
              <div className="p-3 rounded-xl bg-rose-950/60 border border-rose-500/40 text-rose-200 text-xs">
                {rawVcfError}
              </div>
            )}

            {rawVcfData && (
              <div className="overflow-x-auto rounded-xl border border-slate-800">
                <table className="w-full text-left text-xs">
                  <thead className="bg-slate-950 text-slate-400 uppercase text-[11px]">
                    <tr>
                      <th className="p-3">Source File</th>
                      <th className="p-3">Chr:Pos</th>
                      <th className="p-3">Variant ID</th>
                      <th className="p-3">Ref &gt; Alt</th>
                      <th className="p-3">Type</th>
                      <th className="p-3">Your Data (VCF Sample Call)</th>
                      <th className="p-3">Status</th>
                      <th className="p-3">Raw VCF Line</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-800/70 font-mono">
                    {rawVcfData.reads.map((r, i) => (
                      <tr key={i} className="hover:bg-slate-800/40">
                        <td className="p-3 text-cyan-300">{r.source_file}</td>
                        <td className="p-3">
                          chr{r.chr}:{r.position.toLocaleString()}
                          {r.end_position > r.position
                            ? `-${r.end_position.toLocaleString()}`
                            : ''}
                        </td>
                        <td className="p-3 text-white font-semibold">{r.variant_id}</td>
                        <td className="p-3">
                          {r.ref} &rarr; <span className="text-rose-300">{r.alt}</span>
                        </td>
                        <td className="p-3 font-sans">{r.variant_type}</td>
                        <td className="p-3 text-slate-200">{r.your_data}</td>
                        <td className="p-3 font-sans">
                          <span
                            className={`px-2 py-0.5 rounded border text-[11px] ${getStatusBadge(
                              r.your_status
                            )}`}
                          >
                            {r.your_status}
                          </span>
                        </td>
                        <td className="p-3 text-[11px] text-slate-400 max-w-md truncate">
                          {r.raw_vcf_line}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
          </section>
        ) : null}

        {/* NCBI Live Explorer (kept mounted so results survive tab switches) */}
        <div className={activeTab === 'ncbi' ? 'block' : 'hidden'}>
          <NcbiExplorer
            onOpenInExplorer={openInExplorer}
            onCatalogChanged={handleCatalogChanged}
          />
        </div>
      </main>

      {/* Slide-Over Variant & VCF Read Inspector Drawer */}
      {selectedVariant && (
        <div className="fixed inset-0 z-50 flex justify-end bg-slate-950/75 backdrop-blur-sm">
          <div className="w-full max-w-3xl bg-slate-900 border-l border-slate-800 h-full overflow-y-auto p-6 space-y-6 shadow-2xl">
            <div className="flex items-start justify-between gap-4 border-b border-slate-800 pb-4">
              <div>
                <div className="flex flex-wrap items-center gap-2">
                  <span className="px-2.5 py-1 rounded-lg bg-cyan-500/15 border border-cyan-500/30 font-mono font-bold text-cyan-300 text-sm">
                    <RsidLink rsid={selectedVariant.variant_id} className="text-cyan-300 hover:text-cyan-100" />
                  </span>
                  <span className="px-2.5 py-1 rounded-lg bg-slate-800 font-mono text-slate-200 text-xs">
                    <RcvLink rcv={selectedVariant.rcv_id} className="text-slate-200 hover:text-cyan-200" />
                  </span>
                  <span
                    className={`px-2.5 py-1 rounded-lg border text-xs font-semibold ${getClassificationBadge(
                      selectedVariant.classification
                    )}`}
                  >
                    {selectedVariant.classification}
                  </span>
                  <span
                    className={`px-2.5 py-1 rounded-lg border text-xs font-semibold ${getStatusBadge(
                      selectedVariant.your_status
                    )}`}
                  >
                    {selectedVariant.your_status}
                  </span>
                </div>
                <h3 className="text-lg font-bold text-white mt-2 flex flex-wrap items-center gap-x-2 gap-y-1">
                  <GeneLinks gene={selectedVariant.gene} badge={false} />
                  <span>— {selectedVariant.condition}</span>
                </h3>
                <p className="text-xs font-mono text-slate-400 mt-0.5">
                  GRCh38 Coordinate: chr{selectedVariant.chr}:
                  {Number(selectedVariant.position).toLocaleString()}
                  {selectedVariant.hgvs ? ` • ${selectedVariant.hgvs}` : ''}
                </p>
              </div>
              <button
                type="button"
                onClick={() => setSelectedVariant(null)}
                className="p-2 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-300"
              >
                <X className="w-5 h-5" />
              </button>
            </div>

            {/* Drawer tabs: local VCF/catalog view vs live NCBI view */}
            <div className="inline-flex flex-wrap rounded-lg bg-slate-950 p-1 border border-slate-800">
              <button
                type="button"
                onClick={() => setDrawerTab('vcf')}
                className={`px-3 py-1.5 rounded-md text-xs font-medium transition-all flex items-center gap-1.5 ${
                  drawerTab === 'vcf'
                    ? 'bg-cyan-500 text-slate-950 font-semibold'
                    : 'text-slate-400 hover:text-slate-200'
                }`}
              >
                <Terminal className="w-3.5 h-3.5" />
                VCF Read &amp; Catalog
              </button>
              <button
                type="button"
                onClick={() => setDrawerTab('ncbi')}
                className={`px-3 py-1.5 rounded-md text-xs font-medium transition-all flex items-center gap-1.5 ${
                  drawerTab === 'ncbi'
                    ? 'bg-cyan-500 text-slate-950 font-semibold'
                    : 'text-slate-400 hover:text-slate-200'
                }`}
              >
                <Globe className="w-3.5 h-3.5" />
                NCBI (ClinVar · dbSNP · Gene · PubMed)
              </button>
            </div>

            {drawerTab === 'ncbi' && selectedVariant.id !== undefined && selectedVariant.id !== null && (
              <NcbiVariantPanel
                key={selectedVariant.id}
                recordId={selectedVariant.id}
                onSynced={handleNcbiSynced}
              />
            )}

            {drawerTab === 'vcf' && (
            <>
            {/* Key 11-Field Summary Grid */}
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
              <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
                <div className="text-[11px] text-slate-400 uppercase">Your VCF Data</div>
                <div className="text-sm font-mono font-bold text-white mt-1">
                  {selectedVariant.your_data}
                </div>
              </div>
              <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
                <div className="text-[11px] text-slate-400 uppercase">
                  Risk Version (Allele)
                </div>
                <div className="text-sm font-mono font-bold text-rose-300 mt-1">
                  {selectedVariant.risk_version}{' '}
                  <span className="text-xs text-slate-400 font-normal">
                    (Ref: {selectedVariant.ref_allele})
                  </span>
                </div>
              </div>
              <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
                <div className="text-[11px] text-slate-400 uppercase">
                  Population Risk Freq
                </div>
                <div className="text-sm font-mono font-bold text-cyan-300 mt-1">
                  {selectedVariant.risk_frequency_text}
                </div>
              </div>
              <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
                <div className="text-[11px] text-slate-400 uppercase">Variant Type</div>
                <div className="text-sm font-semibold text-slate-200 mt-1">
                  {selectedVariant.variant_type}
                </div>
              </div>
              <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
                <div className="text-[11px] text-slate-400 uppercase">Evidence</div>
                <div className="mt-1 flex items-center gap-2">
                  {renderStars(selectedVariant.evidence_stars)}
                  <span className="text-xs text-slate-300">
                    ({selectedVariant.evidence_stars}/4★)
                  </span>
                </div>
                <div className="text-[11px] text-slate-400 mt-0.5">
                  {selectedVariant.evidence}
                </div>
              </div>
              <div className="p-3 rounded-xl bg-slate-950 border border-slate-800">
                <div className="text-[11px] text-slate-400 uppercase">
                  All Associated RCV IDs
                </div>
                <div className="text-xs font-mono text-slate-300 mt-1 break-words">
                  <RcvList ids={selectedVariant.all_rcv_ids || selectedVariant.rcv_id} />
                </div>
              </div>
            </div>

            {/* Exact raw VCF record */}
            <div className="space-y-2">
              <div className="flex items-center justify-between">
                <h4 className="text-xs font-semibold uppercase tracking-wider text-slate-400">
                  Exact VCF Record from{' '}
                  <span className="text-cyan-400 font-mono">
                    {selectedVariant.source_file}
                  </span>
                </h4>
              </div>
              <pre className="p-3.5 rounded-xl bg-slate-950 border border-slate-800 text-xs font-mono text-cyan-200 overflow-x-auto">
                #CHROM{'\t'}POS{'\t'}ID{'\t'}REF{'\t'}ALT{'\t'}QUAL{'\t'}FILTER{'\t'}INFO{'\t'}FORMAT{'\t'}SAMPLE{'\n'}
                {selectedVariant.raw_vcf_line}
              </pre>
            </div>

            {/* Series of Neighboring VCF Reads Around This Locus */}
            <div className="space-y-2">
              <h4 className="text-xs font-semibold uppercase tracking-wider text-slate-400">
                Series of Neighboring VCF Reads (±2,500 bp around chr
                {selectedVariant.chr}:{Number(selectedVariant.position).toLocaleString()})
              </h4>
              {detailLoading ? (
                <div className="p-4 text-xs text-slate-400">
                  Reading neighboring locus records from VCF...
                </div>
              ) : selectedVariant.context_vcf_reads &&
                selectedVariant.context_vcf_reads.length > 0 ? (
                <div className="rounded-xl border border-slate-800 overflow-x-auto">
                  <table className="w-full text-left text-xs font-mono">
                    <thead className="bg-slate-950 text-slate-400 text-[11px]">
                      <tr>
                        <th className="p-2.5">Coordinate</th>
                        <th className="p-2.5">ID</th>
                        <th className="p-2.5">Ref &gt; Alt</th>
                        <th className="p-2.5">Your Data</th>
                        <th className="p-2.5">Status</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-800/70">
                      {selectedVariant.context_vcf_reads.map((r, idx) => {
                        const isTarget =
                          Number(r.position) === Number(selectedVariant.position) &&
                          r.alt !== '.';
                        return (
                          <tr
                            key={idx}
                            className={
                              isTarget ? 'bg-cyan-500/15 text-cyan-200' : 'text-slate-300'
                            }
                          >
                            <td className="p-2.5">
                              chr{r.chr}:{Number(r.position).toLocaleString()}
                              {r.end_position > r.position
                                ? `-${Number(r.end_position).toLocaleString()}`
                                : ''}
                            </td>
                            <td className="p-2.5 font-semibold">{r.variant_id}</td>
                            <td className="p-2.5">
                              {r.ref} &rarr; {r.alt}
                            </td>
                            <td className="p-2.5">{r.your_data}</td>
                            <td className="p-2.5 font-sans">{r.your_status}</td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              ) : (
                <div className="p-3 text-xs text-slate-500">
                  No additional context reads found.
                </div>
              )}
            </div>
            </>
            )}
          </div>
        </div>
      )}
    </div>
  );
}
