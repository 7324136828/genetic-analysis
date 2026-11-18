import React, { useCallback, useEffect, useState } from 'react';
import { Download, History, Loader2, PlayCircle, RefreshCw, Sparkles, Trash2, XCircle } from 'lucide-react';
import { deleteJob, discardJob, downloadZipUrl, listJobs } from '../api';
import { ACTIVE, StatusBadge, formatElapsed } from '../components/JobStatus';
import { formatBytes } from '../components/FileUpload';

const POLL_MS = 4000;

function formatDate(ts) {
  if (!ts) return '—';
  const d = new Date(ts);
  return Number.isNaN(d.getTime()) ? ts : d.toLocaleString();
}

/**
 * Ledger of all analyses (auto-refreshing). In-progress: Continue / Discard.
 * Completed: Open / Download ZIP / Delete. Failed / discarded: Delete.
 */
export default function HistoryView({ activeJobId, onContinue, onOpenExplorer, onDeleted, onJobsChanged }) {
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [busyId, setBusyId] = useState(null);

  const load = useCallback(async () => {
    try {
      const data = await listJobs();
      setJobs(data);
      setError(null);
      onJobsChanged?.(data);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }, [onJobsChanged]);

  useEffect(() => {
    load();
    const timer = setInterval(load, POLL_MS);
    return () => clearInterval(timer);
  }, [load]);

  const act = async (job, fn, confirmText) => {
    if (confirmText && !window.confirm(confirmText)) return;
    setBusyId(job.id);
    try {
      await fn(job.id);
      if (fn === deleteJob) onDeleted?.(job.id);
      await load();
    } catch (e) {
      setError(e.message);
    } finally {
      setBusyId(null);
    }
  };

  return (
    <section className="rounded-2xl border border-slate-800 bg-slate-900/70 p-5 space-y-4 shadow-xl">
      <div className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <h2 className="text-lg font-bold text-white flex items-center gap-2">
            <History className="w-5 h-5 text-cyan-400" /> Analysis history
          </h2>
          <p className="text-xs text-slate-400">Refreshes automatically every few seconds.</p>
        </div>
        <button
          type="button"
          onClick={load}
          className="px-3 py-2 rounded-lg border border-slate-700 bg-slate-950 hover:bg-slate-800 text-xs text-slate-200 flex items-center gap-1.5"
        >
          <RefreshCw className="w-3.5 h-3.5" /> Refresh
        </button>
      </div>

      {error && <div className="rounded-lg border border-rose-500/40 bg-rose-500/10 p-2.5 text-xs text-rose-100">{error}</div>}

      {loading ? (
        <div className="text-sm text-slate-400 flex items-center gap-2">
          <Loader2 className="w-4 h-4 animate-spin" /> Loading history…
        </div>
      ) : jobs.length === 0 ? (
        <div className="text-sm text-slate-400 py-8 text-center">No analyses yet. Start one from “New Analysis”.</div>
      ) : (
        <div className="overflow-x-auto">
          <table className="w-full text-xs">
            <thead>
              <tr className="text-left text-slate-400 border-b border-slate-800">
                <th className="py-2 pr-3 font-semibold">Analysis</th>
                <th className="py-2 pr-3 font-semibold">Files</th>
                <th className="py-2 pr-3 font-semibold">Size</th>
                <th className="py-2 pr-3 font-semibold">Created</th>
                <th className="py-2 pr-3 font-semibold">Elapsed</th>
                <th className="py-2 pr-3 font-semibold">Status</th>
                <th className="py-2 font-semibold text-right">Actions</th>
              </tr>
            </thead>
            <tbody>
              {jobs.map((job) => {
                const active = ACTIVE.includes(job.status);
                const busy = busyId === job.id;
                return (
                  <tr
                    key={job.id}
                    className={`border-b border-slate-800/70 align-top ${job.id === activeJobId ? 'bg-cyan-500/5' : ''}`}
                  >
                    <td className="py-2.5 pr-3">
                      <div className="font-semibold text-slate-100">{job.label}</div>
                      <div className="font-mono text-[10px] text-slate-500">{job.id.slice(0, 8)}</div>
                    </td>
                    <td className="py-2.5 pr-3 text-slate-300 max-w-[280px]">
                      {job.files.map((f) => (
                        <div key={f.slot} className="truncate" title={f.original_filename}>
                          <span className="text-slate-500">{f.slot === 'snp_indel' ? 'gVCF' : f.slot.toUpperCase()}:</span>{' '}
                          {f.original_filename}
                        </div>
                      ))}
                    </td>
                    <td className="py-2.5 pr-3 font-mono text-slate-300">{formatBytes(job.file_size)}</td>
                    <td className="py-2.5 pr-3 text-slate-300 whitespace-nowrap">{formatDate(job.created_at)}</td>
                    <td className="py-2.5 pr-3 font-mono text-slate-300">{formatElapsed(job.elapsed_seconds)}</td>
                    <td className="py-2.5 pr-3">
                      <StatusBadge status={job.status} progress={job.progress} />
                      {job.status === 'failed' && job.error_message && (
                        <div className="mt-1 text-[10px] text-rose-300 max-w-[220px] truncate" title={job.error_message}>
                          {job.error_message}
                        </div>
                      )}
                    </td>
                    <td className="py-2.5">
                      <div className="flex flex-wrap justify-end gap-1.5">
                        {active && (
                          <>
                            <button
                              type="button"
                              onClick={() => onContinue(job.id)}
                              className="px-2.5 py-1.5 rounded-md bg-cyan-500 text-slate-950 hover:bg-cyan-400 font-semibold flex items-center gap-1"
                            >
                              <PlayCircle className="w-3.5 h-3.5" /> Continue
                            </button>
                            <button
                              type="button"
                              disabled={busy}
                              onClick={() =>
                                act(job, discardJob, 'Discard this in-progress analysis? Its temporary files are deleted immediately.')
                              }
                              className="px-2.5 py-1.5 rounded-md border border-rose-500/50 text-rose-200 hover:bg-rose-500/10 font-semibold flex items-center gap-1 disabled:opacity-50"
                            >
                              <XCircle className="w-3.5 h-3.5" /> Discard
                            </button>
                          </>
                        )}
                        {job.can_open && (
                          <button
                            type="button"
                            onClick={() => onOpenExplorer(job)}
                            className="px-2.5 py-1.5 rounded-md bg-cyan-500 text-slate-950 hover:bg-cyan-400 font-semibold flex items-center gap-1"
                          >
                            <Sparkles className="w-3.5 h-3.5" /> Open
                          </button>
                        )}
                        {job.can_download && (
                          <a
                            href={downloadZipUrl(job.id)}
                            download
                            className="px-2.5 py-1.5 rounded-md border border-slate-700 bg-slate-950 hover:bg-slate-800 text-slate-100 font-semibold flex items-center gap-1"
                          >
                            <Download className="w-3.5 h-3.5 text-cyan-400" /> ZIP
                          </a>
                        )}
                        {!active && (
                          <button
                            type="button"
                            disabled={busy}
                            onClick={() =>
                              act(job, deleteJob, `Delete "${job.label}"? Its uploaded VCFs, catalog and results are removed permanently.`)
                            }
                            className="px-2.5 py-1.5 rounded-md border border-slate-700 text-slate-300 hover:text-rose-200 hover:border-rose-500/50 font-semibold flex items-center gap-1 disabled:opacity-50"
                            aria-label={`Delete ${job.label}`}
                          >
                            <Trash2 className="w-3.5 h-3.5" />
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
      )}
    </section>
  );
}
