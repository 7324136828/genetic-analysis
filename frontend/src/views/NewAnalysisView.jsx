import React, { useEffect, useRef, useState } from 'react';
import { Info, Loader2, Play, PlusCircle, ShieldCheck, Upload, X } from 'lucide-react';
import { createJob } from '../api';
import FileUpload, { SLOTS, assignFiles, formatBytes, isVcfName } from '../components/FileUpload';
import JobStatus from '../components/JobStatus';

const CLINVAR_OPTIONS = [
  {
    value: 'cached',
    label: 'Use the cached ClinVar release',
    hint: 'Verify carried alleles against the local release if one has been downloaded; otherwise skip.',
  },
  {
    value: 'update',
    label: 'Download / refresh the ClinVar release',
    hint: 'Fetches NCBI variant_summary.txt.gz (~460 MB, MD5-verified) once, then reuses it.',
  },
  { value: 'off', label: 'Skip ClinVar verification', hint: 'Fastest; the catalog is not cross-checked allele-by-allele.' },
];

/**
 * Landing screen: pick / drop / paste VCFs, choose options, upload, then follow the
 * analysis live (stages, progress, log) until it can be opened in the explorer.
 */
export default function NewAnalysisView({ config, monitoredJobId, onMonitor, onOpenExplorer, onJobsChanged }) {
  const [files, setFiles] = useState({});
  const [label, setLabel] = useState('');
  const [clinvarMode, setClinvarMode] = useState('cached');
  const [phase, setPhase] = useState('idle'); // idle | uploading
  const [uploadPct, setUploadPct] = useState(0);
  const [error, setError] = useState(null);
  const [fileErrors, setFileErrors] = useState({});
  const uploadRef = useRef(null);

  const maxBytes = (config?.max_upload_mb || 4096) * 1024 * 1024;

  // Ctrl+V: paste VCF files copied in the OS file manager
  useEffect(() => {
    if (monitoredJobId) return undefined;
    const onPaste = (e) => {
      if (phase !== 'idle') return;
      const pasted = Array.from(e.clipboardData?.files || []);
      if (!pasted.length) return;
      e.preventDefault();
      setFiles((prev) => assignFiles(pasted, prev));
    };
    window.addEventListener('paste', onPaste);
    return () => window.removeEventListener('paste', onPaste);
  }, [phase, monitoredJobId]);

  useEffect(() => {
    const errs = {};
    let total = 0;
    SLOTS.forEach((s) => {
      const f = files[s.key];
      if (!f) return;
      total += f.size;
      if (!isVcfName(f.name)) errs[s.key] = 'Expected a .vcf or .vcf.gz file';
      else if (f.size === 0) errs[s.key] = 'The file is empty';
    });
    if (total > maxBytes) errs.snp_indel = `Total upload ${formatBytes(total)} exceeds the ${config?.max_upload_mb} MB limit`;
    setFileErrors(errs);
  }, [files, maxBytes, config]);

  const totalSize = SLOTS.reduce((n, s) => n + (files[s.key]?.size || 0), 0);
  const canStart = phase === 'idle' && files.snp_indel && Object.keys(fileErrors).length === 0;

  const start = async () => {
    if (!canStart) return;
    setError(null);
    setPhase('uploading');
    setUploadPct(0);
    const { promise, abort } = createJob(files, { label: label.trim(), clinvar_mode: clinvarMode }, setUploadPct);
    uploadRef.current = abort;
    try {
      const job = await promise;
      setFiles({});
      setLabel('');
      onMonitor(job.id);
      onJobsChanged?.();
    } catch (e) {
      setError(e.message);
    } finally {
      uploadRef.current = null;
      setPhase('idle');
    }
  };

  if (monitoredJobId) {
    return (
      <div className="space-y-4">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <h2 className="text-lg font-bold text-white">Analysis progress</h2>
          <button
            type="button"
            onClick={() => onMonitor(null)}
            className="px-3 py-2 rounded-lg border border-slate-700 bg-slate-900 hover:bg-slate-800 text-xs font-semibold text-slate-200 flex items-center gap-1.5"
          >
            <PlusCircle className="w-3.5 h-3.5 text-cyan-400" /> Start another analysis
          </button>
        </div>
        <JobStatus key={monitoredJobId} jobId={monitoredJobId} onOpenExplorer={onOpenExplorer} onJobChange={onJobsChanged} />
      </div>
    );
  }

  const cache = config?.clinvar_cache;

  return (
    <div className="space-y-5">
      <section className="rounded-2xl border border-slate-800 bg-slate-900/70 p-5 space-y-5 shadow-xl">
        <div>
          <h2 className="text-lg font-bold text-white flex items-center gap-2">
            <Upload className="w-5 h-5 text-cyan-400" /> New genome analysis
          </h2>
          <p className="text-xs text-slate-400 mt-1">
            Upload a single-sample GRCh38 gVCF (and optionally its structural-variant and copy-number VCFs). GenomeScope indexes it,
            builds an annotated clinical variant catalog, optionally verifies it against a local ClinVar release, and packages the
            results — all inside an isolated folder in your system temp directory.
          </p>
        </div>

        <FileUpload files={files} errors={fileErrors} disabled={phase !== 'idle'} onChange={setFiles} />

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
          <label className="block">
            <span className="text-xs font-semibold text-slate-300">Analysis name (optional)</span>
            <input
              type="text"
              value={label}
              maxLength={80}
              disabled={phase !== 'idle'}
              onChange={(e) => setLabel(e.target.value)}
              placeholder="e.g. My genome – 2026"
              className="mt-1 w-full px-3 py-2 bg-slate-950 border border-slate-700 rounded-lg text-sm text-slate-100 focus:outline-none focus:border-cyan-500"
            />
            <span className="text-[11px] text-slate-500">Shown in History and on exports. Avoid putting identifying details here.</span>
          </label>
          <fieldset className="lg:col-span-2">
            <legend className="text-xs font-semibold text-slate-300">ClinVar release verification</legend>
            <div className="mt-1 grid grid-cols-1 md:grid-cols-3 gap-2">
              {CLINVAR_OPTIONS.map((o) => (
                <label
                  key={o.value}
                  className={`rounded-lg border p-2.5 cursor-pointer text-xs ${
                    clinvarMode === o.value ? 'border-cyan-500 bg-cyan-500/10' : 'border-slate-700 bg-slate-950 hover:border-slate-500'
                  }`}
                >
                  <input
                    type="radio"
                    name="clinvar_mode"
                    value={o.value}
                    checked={clinvarMode === o.value}
                    disabled={phase !== 'idle'}
                    onChange={() => setClinvarMode(o.value)}
                    className="sr-only"
                  />
                  <div className="font-semibold text-slate-100">{o.label}</div>
                  <div className="text-slate-400 mt-0.5 leading-snug">{o.hint}</div>
                </label>
              ))}
            </div>
            <div className="mt-1 text-[11px] text-slate-500">
              {cache
                ? `Cached ClinVar release: ${cache.release_date} (${(cache.rows_indexed || 0).toLocaleString()} GRCh38 alleles).`
                : 'No ClinVar release is cached yet.'}
            </div>
          </fieldset>
        </div>

        <div className="rounded-xl border border-slate-800 bg-slate-950/60 p-3 text-[11px] text-slate-400 flex gap-2">
          <ShieldCheck className="w-4 h-4 text-emerald-400 shrink-0" />
          <div>
            Your VCF files are processed locally and are never uploaded anywhere else. To annotate them, the backend sends variant
            identifiers and genomic coordinates to MyVariant.info, mygene.info and NCBI — including the positions of some of your
            own variant calls (structural / copy-number intervals, gene scans). Start the backend with GENOMESCOPE_OFFLINE=1 to
            disable every external lookup. Analyses are kept in the temp folder
            {config?.job_ttl_hours ? ` for ${config.job_ttl_hours} h after they finish` : ' until you delete them'}; delete them any
            time from History.
          </div>
        </div>

        {error && (
          <div className="rounded-xl border border-rose-500/40 bg-rose-500/10 p-3 text-xs text-rose-100 flex gap-2">
            <Info className="w-4 h-4 shrink-0" /> {error}
          </div>
        )}

        <div className="flex flex-wrap items-center justify-between gap-3">
          <div className="text-xs text-slate-400">
            {files.snp_indel ? `${Object.keys(files).length} file(s) • ${formatBytes(totalSize)}` : 'Add at least the SNP & Indel gVCF.'}
            <span className="text-slate-600"> • limit {config?.max_upload_mb?.toLocaleString() || '…'} MB</span>
          </div>
          {phase === 'uploading' ? (
            <div className="flex items-center gap-3 min-w-[260px]">
              <div className="flex-1">
                <div className="flex justify-between text-[11px] text-slate-400 mb-1">
                  <span className="flex items-center gap-1">
                    <Loader2 className="w-3 h-3 animate-spin" /> {uploadPct < 100 ? 'Uploading' : 'Validating & staging'}
                  </span>
                  <span className="font-mono">{uploadPct}%</span>
                </div>
                <div className="h-2 rounded-full bg-slate-800 overflow-hidden">
                  <div className="h-full bg-cyan-500 transition-all" style={{ width: `${uploadPct}%` }} />
                </div>
              </div>
              <button
                type="button"
                onClick={() => uploadRef.current?.()}
                className="p-2 rounded-lg border border-slate-700 text-slate-300 hover:text-white"
                aria-label="Cancel upload"
              >
                <X className="w-4 h-4" />
              </button>
            </div>
          ) : (
            <button
              type="button"
              disabled={!canStart}
              onClick={start}
              className="px-4 py-2.5 rounded-lg bg-cyan-500 hover:bg-cyan-400 disabled:bg-slate-700 disabled:text-slate-400 text-slate-950 text-sm font-bold flex items-center gap-2"
            >
              <Play className="w-4 h-4" /> Start analysis
            </button>
          )}
        </div>
      </section>
    </div>
  );
}
