import React, { useEffect, useState } from 'react';
import { ShieldCheck } from 'lucide-react';
import { fetchClinvarStatus } from '../api';

const TONES = {
  verified: ['bg-emerald-500/10 text-emerald-300 border-emerald-500/40', 'bg-emerald-400'],
  stale: ['bg-amber-500/10 text-amber-300 border-amber-500/40', 'bg-amber-400'],
  never_run: ['bg-amber-500/10 text-amber-300 border-amber-500/40', 'bg-amber-400'],
  index_missing: ['bg-amber-500/10 text-amber-300 border-amber-500/40', 'bg-amber-400'],
  no_analysis: ['bg-slate-800 text-slate-300 border-slate-700', 'bg-emerald-400'],
};

const LABELS = {
  verified: 'ClinVar verified',
  stale: 'ClinVar stale',
  never_run: 'ClinVar unscanned',
  index_missing: 'ClinVar not loaded',
  no_analysis: 'ClinVar ready',
};

/** Header pill: is the catalog verified against a local ClinVar release (GET /api/clinvar/status, once on load). */
export default function ClinvarStatusPill() {
  const [status, setStatus] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let active = true;
    fetchClinvarStatus()
      .then((s) => {
        if (active) setStatus(s);
      })
      .catch((err) => {
        if (active) setError(err.message || 'ClinVar status unavailable');
      });
    return () => {
      active = false;
    };
  }, []);

  let tone = 'bg-slate-800 text-slate-400 border-slate-700';
  let dot = 'bg-slate-500 animate-pulse';
  let label = 'ClinVar…';
  if (error) {
    [tone, dot] = TONES.stale;
    label = 'ClinVar ?';
  } else if (status) {
    [tone, dot] = TONES[status.state] || TONES.stale;
    label = LABELS[status.state] || 'ClinVar';
  }

  const scan = status?.last_scan;
  const lines = error
    ? [`Status check failed: ${error}`]
    : status
    ? [
        status.message || label,
        status.index ? `ClinVar release ${status.index.release_date} (${(status.index.rows_indexed || 0).toLocaleString()} alleles)` : null,
        scan?.finished_at && !scan.dry_run ? `Last scan ${scan.finished_at.replace('T', ' ').replace('+00:00', ' UTC')}` : null,
        scan?.coverage
          ? `Pathogenic/likely-pathogenic loci: ${scan.coverage.confirmed_reference.toLocaleString()} confirmed reference, ${scan.coverage.variant_called.toLocaleString()} variant rows, ${scan.coverage.no_gvcf_record.toLocaleString()} not covered`
          : null,
        status.catalog?.carried_rows_without_verified_allele
          ? `${status.catalog.carried_rows_without_verified_allele} carried row(s) lack a verified allele`
          : null,
        scan?.unresolved_high_impact_indels ? `${scan.unresolved_high_impact_indels} indel(s) need manual review (see scan_report.json)` : null,
      ].filter(Boolean)
    : ['Checking ClinVar scan status…'];

  return (
    <div className="relative group">
      <span
        tabIndex={0}
        aria-label={lines.join('. ')}
        className={`inline-flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg border text-[11px] font-semibold cursor-default focus:outline-none focus:ring-1 focus:ring-cyan-500 ${tone}`}
      >
        <span className={`w-1.5 h-1.5 rounded-full ${dot}`} />
        <ShieldCheck className="w-3.5 h-3.5" />
        {label}
      </span>
      <div
        role="tooltip"
        className="pointer-events-none absolute right-0 top-full mt-2 w-72 rounded-xl border border-slate-700 bg-slate-900 p-3 text-[11px] text-slate-300 shadow-2xl opacity-0 group-hover:opacity-100 group-focus-within:opacity-100 transition-opacity z-40 space-y-0.5"
      >
        {lines.map((line, i) => (
          <div key={i} className={i === 0 ? 'font-semibold text-slate-100' : ''}>
            {line}
          </div>
        ))}
      </div>
    </div>
  );
}
