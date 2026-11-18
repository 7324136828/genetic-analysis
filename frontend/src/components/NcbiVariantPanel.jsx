import React, { useCallback, useEffect, useRef, useState } from 'react';
import { AlertTriangle, CheckCircle2, GitCompare, RefreshCw } from 'lucide-react';
import { fetchVariantNcbi } from '../api';
import { ErrorBanner, NcbiVariantDetail, SkeletonCards, asArray, formatFrequency } from './ncbiShared';

const FIELD_LABELS = {
  allele: 'Allele (ClinVar ALT)',
  classification: 'Classification',
  review_status: 'Review status',
  condition: 'Condition',
  rcv_id: 'RCV ID',
  risk_frequency: 'Risk allele frequency',
};

function formatCompareValue(field, value) {
  if (value === null || value === undefined || value === '') return '—';
  if (field === 'risk_frequency' && typeof value === 'number') return `${formatFrequency(value)} (${value})`;
  if (Array.isArray(value)) return value.length ? value.map((v) => (typeof v === 'object' ? JSON.stringify(v) : String(v))).join(', ') : '—';
  if (typeof value === 'object') {
    if (field === 'risk_frequency' && value.frequency !== undefined) {
      return `${formatFrequency(value.frequency)}${value.allele ? ` (${value.allele})` : ''}${value.study ? ` · ${value.study}` : ''}`;
    }
    return JSON.stringify(value);
  }
  if (typeof value === 'boolean') return value ? 'yes' : 'no';
  return String(value);
}

function CatalogVsNcbiTable({ rows }) {
  const list = asArray(rows);
  if (list.length === 0) return null;
  const mismatches = list.filter((r) => r && r.match === false).length;
  return (
    <div className="rounded-xl border border-slate-800 bg-slate-950 overflow-hidden">
      <div className="px-4 py-2.5 flex items-center justify-between border-b border-slate-800">
        <h4 className="text-xs font-semibold uppercase tracking-wider text-slate-400 flex items-center gap-1.5">
          <GitCompare className="w-3.5 h-3.5 text-cyan-400" />
          Catalog vs NCBI
        </h4>
        <span
          className={`text-[11px] font-semibold px-2 py-0.5 rounded-md border ${
            mismatches > 0
              ? 'bg-amber-500/15 text-amber-200 border-amber-500/40'
              : 'bg-emerald-500/15 text-emerald-300 border-emerald-500/40'
          }`}
        >
          {mismatches > 0 ? `${mismatches} difference${mismatches === 1 ? '' : 's'}` : 'All fields match'}
        </span>
      </div>
      <div className="overflow-x-auto">
        <table className="w-full text-left text-xs">
          <thead className="bg-slate-900 text-[10px] uppercase tracking-wider text-slate-500">
            <tr>
              <th className="px-3 py-2">Field</th>
              <th className="px-3 py-2">Local catalog</th>
              <th className="px-3 py-2">NCBI (live)</th>
              <th className="px-3 py-2 w-16">Match</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-800/70">
            {list.map((r, i) => {
              const mismatch = r.match === false;
              return (
                <tr key={`${r.field}-${i}`} className={mismatch ? 'bg-amber-500/10' : ''}>
                  <td className="px-3 py-2 text-slate-300 font-medium whitespace-nowrap">{FIELD_LABELS[r.field] || r.field}</td>
                  <td className="px-3 py-2 text-slate-300 break-words max-w-[220px]">{formatCompareValue(r.field, r.catalog)}</td>
                  <td className={`px-3 py-2 break-words max-w-[220px] ${mismatch ? 'text-amber-200 font-semibold' : 'text-slate-200'}`}>
                    {formatCompareValue(r.field, r.ncbi)}
                  </td>
                  <td className="px-3 py-2">
                    {r.match === true ? (
                      <span className="text-emerald-300">yes</span>
                    ) : r.match === false ? (
                      <span className="text-amber-300 font-semibold">no</span>
                    ) : (
                      <span className="text-slate-500">n/a</span>
                    )}
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

/**
 * NCBI tab inside the variant inspector drawer. Lazily loads
 * GET /api/variants/{recordId}/ncbi when mounted, and can sync the catalog row.
 */
export default function NcbiVariantPanel({ recordId, onSynced }) {
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [syncing, setSyncing] = useState(false);
  const [syncNotice, setSyncNotice] = useState(null);
  const requestSeq = useRef(0);
  // The panel unmounts when the drawer tab changes; a sync that already reached the server
  // must still refresh the parent views (table, summary, selected row).
  const onSyncedRef = useRef(onSynced);
  useEffect(() => {
    onSyncedRef.current = onSynced;
  }, [onSynced]);

  const load = useCallback(
    async (sync = false) => {
      if (recordId === null || recordId === undefined) return;
      const seq = ++requestSeq.current;
      if (sync) setSyncing(true);
      else setLoading(true);
      setError(null);
      try {
        const body = await fetchVariantNcbi(recordId, { sync });
        if (sync && body?.synced && onSyncedRef.current) onSyncedRef.current(body);
        if (seq !== requestSeq.current) return;
        setData(body);
        if (sync) {
          if (body?.synced) {
            setSyncNotice({ ok: true, text: 'Catalog row updated from NCBI ClinVar / dbSNP.' });
          } else {
            const reason = (Array.isArray(body?.errors) ? body.errors : []).find((e) => String(e).startsWith('Sync skipped'));
            setSyncNotice({ ok: false, text: reason || 'NCBI returned no data to sync; catalog row unchanged.' });
          }
        }
      } catch (err) {
        if (seq !== requestSeq.current) return;
        setError(err.message || 'Failed to load NCBI data');
      } finally {
        if (seq === requestSeq.current) {
          setLoading(false);
          setSyncing(false);
        }
      }
    },
    [recordId]
  );

  useEffect(() => {
    setData(null);
    setSyncNotice(null);
    load(false);
    return () => {
      // invalidate in-flight responses when the record changes / panel unmounts
      requestSeq.current += 1;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [recordId]);

  const comparisons = asArray(data?.catalog_vs_ncbi);
  const hasMismatch = comparisons.some((r) => r && r.match === false);
  const canSync = asArray(data?.clinvar).length > 0;

  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <p className="text-xs text-slate-400 max-w-md">
          Live clinical annotation from NCBI ClinVar, dbSNP, Gene and PubMed, joined with your own VCF read at the
          same GRCh38 locus.
        </p>
        <div className="flex items-center gap-2">
          <button
            type="button"
            onClick={() => load(false)}
            disabled={loading || syncing}
            className="px-3 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-xs font-medium flex items-center gap-1.5 disabled:opacity-50"
          >
            <RefreshCw className={`w-3.5 h-3.5 ${loading ? 'animate-spin' : ''}`} />
            Reload
          </button>
          <button
            type="button"
            onClick={() => load(true)}
            disabled={loading || syncing || !data || !canSync}
            title={
              canSync
                ? "Overwrite this catalog row's classification, review status, evidence, condition, RCV ids and risk frequency with NCBI values"
                : 'No NCBI ClinVar record to sync from'
            }
            className={`px-3 py-1.5 rounded-lg text-xs font-semibold flex items-center gap-1.5 disabled:opacity-50 ${
              hasMismatch
                ? 'bg-amber-500 hover:bg-amber-400 text-slate-950'
                : 'bg-cyan-500 hover:bg-cyan-400 text-slate-950'
            }`}
          >
            <RefreshCw className={`w-3.5 h-3.5 ${syncing ? 'animate-spin' : ''}`} />
            {syncing ? 'Syncing…' : 'Sync catalog from NCBI'}
          </button>
        </div>
      </div>

      {syncNotice && (
        <div
          className={`p-2.5 rounded-xl border text-xs flex items-center gap-2 ${
            syncNotice.ok
              ? 'bg-emerald-950/40 border-emerald-500/30 text-emerald-200'
              : 'bg-amber-950/40 border-amber-500/30 text-amber-200'
          }`}
        >
          {syncNotice.ok ? <CheckCircle2 className="w-3.5 h-3.5 shrink-0" /> : <AlertTriangle className="w-3.5 h-3.5 shrink-0" />}
          {syncNotice.text}
        </div>
      )}

      <ErrorBanner message={error} onRetry={() => load(false)} />

      {loading && !data ? (
        <SkeletonCards count={4} />
      ) : data ? (
        <>
          <CatalogVsNcbiTable rows={comparisons} />
          <NcbiVariantDetail data={data} />
        </>
      ) : null}
    </div>
  );
}
