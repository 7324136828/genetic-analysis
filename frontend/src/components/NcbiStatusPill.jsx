import React, { useEffect, useState } from 'react';
import { Globe } from 'lucide-react';
import { fetchNcbiStatus } from '../api';
import { asArray, formatFetchedAt } from './ncbiShared';

/** Small header pill showing NCBI E-utilities reachability (GET /api/ncbi/status, once on load). */
export default function NcbiStatusPill() {
  const [status, setStatus] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let active = true;
    fetchNcbiStatus()
      .then((s) => {
        if (active) setStatus(s);
      })
      .catch((err) => {
        if (active) setError(err.message || 'NCBI status unavailable');
      });
    return () => {
      active = false;
    };
  }, []);

  let tone = 'bg-slate-800 text-slate-400 border-slate-700';
  let dot = 'bg-slate-500 animate-pulse';
  let label = 'NCBI…';
  if (error) {
    tone = 'bg-amber-500/10 text-amber-300 border-amber-500/40';
    dot = 'bg-amber-400';
    label = 'NCBI ?';
  } else if (status) {
    if (status.reachable) {
      tone = 'bg-emerald-500/10 text-emerald-300 border-emerald-500/40';
      dot = 'bg-emerald-400';
      label = 'NCBI live';
    } else {
      tone = 'bg-amber-500/10 text-amber-300 border-amber-500/40';
      dot = 'bg-amber-400';
      label = 'NCBI offline';
    }
  }

  const tooltipLines = error
    ? [`Status check failed: ${error}`]
    : status
    ? [
        status.reachable ? 'NCBI E-utilities reachable' : 'NCBI unreachable — cached data only',
        `Rate limit: ${status.rate_limit_per_sec ?? '?'} req/s`,
        `API key: ${status.api_key_configured ? 'configured' : 'not configured'}`,
        `Cache: ${status.cache_entries ?? 0} entries, TTL ${status.cache_ttl_hours ?? '?'} h`,
        ...asArray(status.errors).map((e) => `Note: ${e}`),
        status.fetched_at ? `Checked ${formatFetchedAt(status.fetched_at)}` : null,
      ].filter(Boolean)
    : ['Checking NCBI connectivity…'];

  return (
    <div className="relative group">
      <span
        tabIndex={0}
        aria-label={tooltipLines.join('. ')}
        className={`inline-flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg border text-[11px] font-semibold cursor-default focus:outline-none focus:ring-1 focus:ring-cyan-500 ${tone}`}
      >
        <span className={`w-1.5 h-1.5 rounded-full ${dot}`} />
        <Globe className="w-3.5 h-3.5" />
        {label}
      </span>
      <div
        role="tooltip"
        className="pointer-events-none absolute right-0 top-full mt-2 w-64 rounded-xl border border-slate-700 bg-slate-900 p-3 text-[11px] text-slate-300 shadow-2xl opacity-0 group-hover:opacity-100 group-focus-within:opacity-100 transition-opacity z-40 space-y-0.5"
      >
        {tooltipLines.map((line, i) => (
          <div key={i} className={i === 0 ? 'font-semibold text-slate-100' : ''}>
            {line}
          </div>
        ))}
      </div>
    </div>
  );
}
