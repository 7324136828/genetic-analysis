import React, { useEffect, useState } from 'react';
import {
  AlertTriangle,
  CheckCircle2,
  Circle,
  Download,
  Loader2,
  MinusCircle,
  Sparkles,
  Trash2,
  XCircle,
} from 'lucide-react';
import { discardJob, downloadZipUrl, getJob } from '../api';
import LogViewer from './LogViewer';
import { formatBytes } from './FileUpload';

export const ACTIVE = ['queued', 'in_progress'];

const STATUS_STYLES = {
  queued: 'bg-slate-700/40 text-slate-200 border-slate-600',
  in_progress: 'bg-cyan-500/15 text-cyan-200 border-cyan-500/40',
  completed: 'bg-emerald-500/15 text-emerald-200 border-emerald-500/40',
  failed: 'bg-rose-500/15 text-rose-200 border-rose-500/40',
  discarded: 'bg-slate-800 text-slate-400 border-slate-700',
};
const STATUS_LABELS = {
  queued: 'Queued',
  in_progress: 'In progress',
  completed: 'Completed',
  failed: 'Failed',
  discarded: 'Discarded',
};

export function StatusBadge({ status, progress }) {
  return (
    <span
      className={`inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full border text-[11px] font-semibold ${
        STATUS_STYLES[status] || STATUS_STYLES.queued
      }`}
    >
      {ACTIVE.includes(status) && <span className="w-1.5 h-1.5 rounded-full bg-cyan-300 animate-pulse" />}
      {STATUS_LABELS[status] || status}
      {status === 'in_progress' && typeof progress === 'number' ? ` ${progress}%` : ''}
    </span>
  );
}

export function formatElapsed(seconds) {
  if (seconds == null) return '—';
  const s = Math.max(0, Math.round(seconds));
  if (s < 60) return `${s}s`;
  const m = Math.floor(s / 60);
  if (m < 60) return `${m}m ${s % 60}s`;
  return `${Math.floor(m / 60)}h ${m % 60}m`;
}

function StageIcon({ status }) {
  if (status === 'done') return <CheckCircle2 className="w-4 h-4 text-emerald-400" />;
  if (status === 'running') return <Loader2 className="w-4 h-4 text-cyan-300 animate-spin" />;
  if (status === 'failed') return <XCircle className="w-4 h-4 text-rose-400" />;
  if (status === 'skipped') return <MinusCircle className="w-4 h-4 text-slate-500" />;
  return <Circle className="w-4 h-4 text-slate-600" />;
}

/**
 * Live monitor of one analysis: multi-stage progress, messages, log, and the actions that
 * fit its state (discard while running; open / download when completed).
 */
export default function JobStatus({ jobId, initialJob = null, onOpenExplorer, onJobChange }) {
  const [job, setJob] = useState(initialJob);
  const [error, setError] = useState(null);
  const [busy, setBusy] = useState(false);

  useEffect(() => {
    let alive = true;
    let timer = null;
    const load = async () => {
      try {
        const j = await getJob(jobId);
        if (!alive) return;
        setJob(j);
        setError(null);
        onJobChange?.(j);
        if (ACTIVE.includes(j.status)) timer = setTimeout(load, 1500);
      } catch (e) {
        if (!alive) return;
        setError(e.message);
        timer = setTimeout(load, 4000);
      }
    };
    load();
    return () => {
      alive = false;
      if (timer) clearTimeout(timer);
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId]);

  const handleDiscard = async () => {
    if (!window.confirm('Discard this analysis? The running step is stopped and all of its temporary files are deleted.')) return;
    setBusy(true);
    try {
      const j = await discardJob(jobId);
      setJob(j);
      onJobChange?.(j);
    } catch (e) {
      setError(e.message);
    } finally {
      setBusy(false);
    }
  };

  if (!job) {
    return (
      <div className="rounded-2xl border border-slate-800 bg-slate-900/70 p-6 text-sm text-slate-400 flex items-center gap-2">
        <Loader2 className="w-4 h-4 animate-spin" /> {error ? `Loading analysis… (${error})` : 'Loading analysis…'}
      </div>
    );
  }

  const active = ACTIVE.includes(job.status);
  const pkg = job.result?.package || {};

  return (
    <div className="rounded-2xl border border-slate-800 bg-slate-900/70 p-5 space-y-4 shadow-xl">
      <div className="flex flex-wrap items-start justify-between gap-3">
        <div className="min-w-0">
          <div className="flex items-center gap-2 flex-wrap">
            <h2 className="text-base font-bold text-white truncate">{job.label}</h2>
            <StatusBadge status={job.status} progress={job.progress} />
          </div>
          <p className="text-xs text-slate-400 mt-0.5">
            {job.files.map((f) => f.original_filename).join(' • ') || 'No files'} • {formatBytes(job.file_size)} •
            elapsed {formatElapsed(job.elapsed_seconds)}
          </p>
        </div>
        <div className="flex flex-wrap items-center gap-2">
          {active && (
            <button
              type="button"
              disabled={busy}
              onClick={handleDiscard}
              className="px-3 py-2 rounded-lg border border-rose-500/50 text-rose-200 hover:bg-rose-500/10 text-xs font-semibold flex items-center gap-1.5 disabled:opacity-50"
            >
              <Trash2 className="w-3.5 h-3.5" /> Discard
            </button>
          )}
          {job.can_open && (
            <button
              type="button"
              onClick={() => onOpenExplorer?.(job)}
              className="px-3 py-2 rounded-lg bg-cyan-500 text-slate-950 hover:bg-cyan-400 text-xs font-semibold flex items-center gap-1.5"
            >
              <Sparkles className="w-3.5 h-3.5" /> Open in Explorer
            </button>
          )}
          {job.can_download && (
            <a
              href={downloadZipUrl(job.id)}
              download
              className="px-3 py-2 rounded-lg border border-slate-700 bg-slate-950 hover:bg-slate-800 text-slate-100 text-xs font-semibold flex items-center gap-1.5"
            >
              <Download className="w-3.5 h-3.5 text-cyan-400" /> Download ZIP
            </a>
          )}
        </div>
      </div>

      <div>
        <div className="flex items-center justify-between text-[11px] text-slate-400 mb-1">
          <span>{job.message || STATUS_LABELS[job.status]}</span>
          <span className="font-mono">{job.progress}%</span>
        </div>
        <div className="h-2 rounded-full bg-slate-800 overflow-hidden" role="progressbar" aria-valuenow={job.progress} aria-valuemin={0} aria-valuemax={100}>
          <div
            className={`h-full transition-all duration-500 ${
              job.status === 'failed' ? 'bg-rose-500' : job.status === 'discarded' ? 'bg-slate-600' : 'bg-gradient-to-r from-cyan-500 to-blue-500'
            } ${active ? 'animate-pulse' : ''}`}
            style={{ width: `${Math.max(2, job.progress)}%` }}
          />
        </div>
      </div>

      <ol className="grid grid-cols-1 sm:grid-cols-5 gap-2">
        {job.stages.map((s, i) => (
          <li
            key={s.key}
            className={`rounded-lg border px-3 py-2 ${
              s.status === 'running' ? 'border-cyan-500/50 bg-cyan-500/5' : 'border-slate-800 bg-slate-950/50'
            }`}
          >
            <div className="flex items-center gap-1.5 text-xs font-semibold text-slate-200">
              <StageIcon status={s.status} />
              <span>
                {i + 1}. {s.label}
              </span>
            </div>
            {s.message && <div className="mt-1 text-[11px] text-slate-400 leading-snug">{s.message}</div>}
          </li>
        ))}
      </ol>

      {job.status === 'failed' && (
        <div className="rounded-xl border border-rose-500/40 bg-rose-500/10 p-3 text-xs text-rose-100 flex gap-2">
          <AlertTriangle className="w-4 h-4 shrink-0 text-rose-300" />
          <div>
            <div className="font-semibold">The analysis failed.</div>
            <div className="mt-0.5 font-mono break-words">{job.error_message || 'Unknown error'}</div>
          </div>
        </div>
      )}
      {job.status === 'completed' && (
        <div className="rounded-xl border border-emerald-500/30 bg-emerald-500/5 p-3 text-xs text-emerald-100">
          Catalog: <b>{(pkg.catalog_rows ?? 0).toLocaleString()}</b> records, <b>{(pkg.carried_rows ?? 0).toLocaleString()}</b>{' '}
          carried ({pkg.pathogenic_or_likely_carried ?? 0} pathogenic / likely pathogenic, {pkg.drug_response_carried ?? 0} drug
          response). ClinVar verification: {job.result?.clinvar?.status === 'done' ? `release ${job.result.clinvar.clinvar_release}` : 'skipped'}.
        </div>
      )}
      {error && <div className="text-xs text-amber-300">{error}</div>}

      <LogViewer jobId={job.id} active={active} defaultOpen={job.status === 'failed'} />
    </div>
  );
}
