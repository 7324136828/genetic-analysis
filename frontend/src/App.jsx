import React, { useCallback, useEffect, useState } from 'react';
import { Dna, History, PlusCircle, Sparkles, X } from 'lucide-react';
import { fetchConfig, getJob, listJobs, setActiveJob } from './api';
import ExplorerView from './views/ExplorerView';
import HistoryView from './views/HistoryView';
import NewAnalysisView from './views/NewAnalysisView';

const STORE_KEY = 'genomescope.ui.v1';

// Per-viewer convenience only (which screen / analysis was open); never genome data.
function loadUiState() {
  try {
    return JSON.parse(window.localStorage.getItem(STORE_KEY) || '{}') || {};
  } catch {
    return {};
  }
}

function saveUiState(state) {
  try {
    window.localStorage.setItem(STORE_KEY, JSON.stringify(state));
  } catch {
    /* storage unavailable (private mode): ignore */
  }
}

const TABS = [
  { key: 'new', label: 'New Analysis', icon: PlusCircle },
  { key: 'history', label: 'History', icon: History },
  { key: 'explorer', label: 'Explorer', icon: Sparkles },
];

export default function App() {
  const saved = loadUiState();
  const [view, setView] = useState(saved.view === 'explorer' ? 'history' : saved.view || 'new');
  const [monitoredJobId, setMonitoredJobId] = useState(saved.monitoredJobId || null);
  const [activeJob, setActiveJobState] = useState(null);
  const [config, setConfig] = useState(null);
  const [activeCount, setActiveCount] = useState(0);

  useEffect(() => {
    fetchConfig().then(setConfig).catch(() => setConfig(null));
  }, []);

  // Restore the previously opened analysis if it is still available
  useEffect(() => {
    if (!saved.activeJobId) return;
    getJob(saved.activeJobId)
      .then((job) => {
        if (job.can_open) {
          setActiveJob(job.id);
          setActiveJobState(job);
          if (saved.view === 'explorer') setView('explorer');
        }
      })
      .catch(() => {});
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  useEffect(() => {
    saveUiState({ view, monitoredJobId, activeJobId: activeJob?.id || null });
  }, [view, monitoredJobId, activeJob]);

  const refreshJobs = useCallback((data) => {
    const apply = (list) => setActiveCount(list.filter((j) => j.status === 'queued' || j.status === 'in_progress').length);
    if (Array.isArray(data)) apply(data);
    else listJobs().then(apply).catch(() => {});
  }, []);

  useEffect(() => {
    refreshJobs();
    const t = setInterval(refreshJobs, 5000);
    return () => clearInterval(t);
  }, [refreshJobs]);

  const openExplorer = useCallback((job) => {
    setActiveJob(job.id); // before the explorer mounts, so its first requests carry the id
    setActiveJobState(job);
    setView('explorer');
  }, []);

  const closeAnalysis = useCallback(() => {
    setActiveJob(null);
    setActiveJobState(null);
    setView((v) => (v === 'explorer' ? 'history' : v));
  }, []);

  const continueJob = useCallback((jobId) => {
    setMonitoredJobId(jobId);
    setView('new');
  }, []);

  const onDeleted = useCallback(
    (jobId) => {
      if (activeJob?.id === jobId) closeAnalysis();
      if (monitoredJobId === jobId) setMonitoredJobId(null);
    },
    [activeJob, monitoredJobId, closeAnalysis]
  );

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 flex flex-col">
      <nav className="sticky top-0 z-40 h-14 border-b border-slate-800 bg-slate-950/95 backdrop-blur-md">
        <div className="max-w-[1600px] h-full mx-auto px-4 sm:px-6 flex items-center gap-4">
          <div className="flex items-center gap-2.5 shrink-0">
            <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center shadow-lg shadow-cyan-500/20">
              <Dna className="w-5 h-5 text-white" />
            </div>
            <span className="font-bold tracking-tight text-white hidden sm:inline">GenomeScope</span>
          </div>
          <div className="flex items-center gap-1 overflow-x-auto" role="tablist">
            {TABS.map(({ key, label, icon: Icon }) => {
              const disabled = key === 'explorer' && !activeJob;
              return (
                <button
                  key={key}
                  type="button"
                  role="tab"
                  aria-selected={view === key}
                  disabled={disabled}
                  title={disabled ? 'Open a completed analysis from History first' : undefined}
                  onClick={() => setView(key)}
                  className={`px-3 py-1.5 rounded-md text-xs font-semibold flex items-center gap-1.5 whitespace-nowrap transition ${
                    view === key
                      ? 'bg-cyan-500 text-slate-950'
                      : disabled
                      ? 'text-slate-600 cursor-not-allowed'
                      : 'text-slate-300 hover:text-white hover:bg-slate-800'
                  }`}
                >
                  <Icon className="w-3.5 h-3.5" />
                  {label}
                  {key === 'history' && activeCount > 0 && (
                    <span className="ml-0.5 px-1.5 rounded-full bg-cyan-400/20 text-cyan-200 text-[10px] animate-pulse">{activeCount}</span>
                  )}
                </button>
              );
            })}
          </div>
          {activeJob && (
            <div className="ml-auto hidden md:flex items-center gap-1.5 px-2.5 py-1 rounded-lg border border-slate-700 bg-slate-900 text-[11px] text-slate-300 min-w-0">
              <span className="text-slate-500">Analysis:</span>
              <span className="font-semibold truncate max-w-[220px]" title={activeJob.label}>
                {activeJob.label}
              </span>
              <button type="button" onClick={closeAnalysis} className="p-0.5 rounded hover:bg-slate-800" aria-label="Close analysis">
                <X className="w-3.5 h-3.5" />
              </button>
            </div>
          )}
        </div>
      </nav>

      {view === 'explorer' && activeJob ? (
        <ExplorerView key={activeJob.id} job={activeJob} />
      ) : (
        <main className="flex-1 max-w-[1600px] w-full mx-auto px-4 sm:px-6 py-6">
          {view === 'history' ? (
            <HistoryView
              activeJobId={activeJob?.id}
              onContinue={continueJob}
              onOpenExplorer={openExplorer}
              onDeleted={onDeleted}
              onJobsChanged={refreshJobs}
            />
          ) : (
            <NewAnalysisView
              config={config}
              monitoredJobId={monitoredJobId}
              onMonitor={setMonitoredJobId}
              onOpenExplorer={openExplorer}
              onJobsChanged={refreshJobs}
            />
          )}
        </main>
      )}

      <footer className="border-t border-slate-900 py-3 text-center text-[11px] text-slate-600">
        Research-grade interpretation of sequencing data — not a medical diagnosis. Confirm clinically relevant findings with an
        accredited laboratory.
      </footer>
    </div>
  );
}
