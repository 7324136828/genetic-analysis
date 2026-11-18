import React, { useEffect, useRef, useState } from 'react';
import { ChevronDown, ChevronRight, Terminal } from 'lucide-react';
import { getJobLogs } from '../api';

/** Collapsible tail of the pipeline log; polls while the job is active. */
export default function LogViewer({ jobId, active, defaultOpen = false }) {
  const [open, setOpen] = useState(defaultOpen);
  const [lines, setLines] = useState([]);
  const [error, setError] = useState(null);
  const boxRef = useRef(null);

  useEffect(() => {
    if (!open || !jobId) return undefined;
    let alive = true;
    const load = () =>
      getJobLogs(jobId, 400)
        .then((b) => {
          if (!alive) return;
          setLines(b.lines || []);
          setError(null);
        })
        .catch((e) => alive && setError(e.message));
    load();
    const timer = active ? setInterval(load, 2500) : null;
    return () => {
      alive = false;
      if (timer) clearInterval(timer);
    };
  }, [open, jobId, active]);

  useEffect(() => {
    if (boxRef.current) boxRef.current.scrollTop = boxRef.current.scrollHeight;
  }, [lines]);

  return (
    <div className="rounded-xl border border-slate-800 bg-slate-950/70">
      <button
        type="button"
        onClick={() => setOpen((o) => !o)}
        className="w-full flex items-center gap-2 px-3 py-2 text-xs font-semibold text-slate-300 hover:text-white"
        aria-expanded={open}
      >
        {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
        <Terminal className="w-3.5 h-3.5 text-cyan-400" />
        Pipeline log
        {active && open && <span className="ml-auto text-[10px] text-cyan-400 animate-pulse">live</span>}
      </button>
      {open && (
        <pre
          ref={boxRef}
          className="max-h-72 overflow-auto px-3 pb-3 text-[11px] leading-relaxed font-mono text-slate-300 whitespace-pre-wrap break-words"
        >
          {error ? `Could not load the log: ${error}` : lines.length ? lines.join('\n') : 'No log output yet…'}
        </pre>
      )}
    </div>
  );
}
