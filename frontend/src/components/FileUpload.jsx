import React, { useRef, useState } from 'react';
import { FileUp, FileCheck2, X, Layers, Dna, Activity } from 'lucide-react';

export const SLOTS = [
  {
    key: 'snp_indel',
    label: 'SNP & Indel gVCF',
    hint: 'Required • genome VCF with hom-ref blocks (.vcf or .vcf.gz)',
    icon: Dna,
    required: true,
  },
  {
    key: 'sv',
    label: 'Structural Variants VCF',
    hint: 'Optional • e.g. Manta deletions / insertions / duplications',
    icon: Layers,
    required: false,
  },
  {
    key: 'cnv',
    label: 'Copy Number Variants VCF',
    hint: 'Optional • e.g. Canvas gains / losses / LOH',
    icon: Activity,
    required: false,
  },
];

const ACCEPT = '.vcf,.vcf.gz,.gz,.bgz';

export function formatBytes(n) {
  const v = Number(n) || 0;
  if (v >= 1024 ** 3) return `${(v / 1024 ** 3).toFixed(2)} GB`;
  if (v >= 1024 ** 2) return `${(v / 1024 ** 2).toFixed(1)} MB`;
  if (v >= 1024) return `${(v / 1024).toFixed(0)} KB`;
  return `${v} B`;
}

export function isVcfName(name) {
  return /\.(vcf|vcf\.gz|vcf\.bgz|gz|bgz)$/i.test(String(name || ''));
}

/** Guess the slot of a file from its name (used for multi-file drops and Ctrl+V paste). */
export function guessSlot(name) {
  const n = String(name || '').toLowerCase();
  if (/(^|[._-])(cnv|canvas|copy[-_]?number)([._-]|$)/.test(n)) return 'cnv';
  if (/(^|[._-])(sv|manta|structural)([._-]|$)/.test(n)) return 'sv';
  return 'snp_indel';
}

/** Assign several files to slots: by name first, then into any free slot. */
export function assignFiles(fileList, current = {}) {
  const next = { ...current };
  const rest = [];
  Array.from(fileList || []).forEach((f) => {
    const slot = guessSlot(f.name);
    if (!next[slot] || fileList.length === 1) next[slot] = f;
    else rest.push(f);
  });
  rest.forEach((f) => {
    const free = SLOTS.find((s) => !next[s.key]);
    if (free) next[free.key] = f;
  });
  return next;
}

function SlotDropzone({ slot, file, error, disabled, onFile, onClear }) {
  const inputRef = useRef(null);
  const [over, setOver] = useState(false);
  const Icon = slot.icon;

  const pick = (files) => {
    const f = files && files[0];
    if (f) onFile(slot.key, f);
  };

  return (
    <div
      onDragOver={(e) => {
        e.preventDefault();
        e.stopPropagation();
        if (!disabled) setOver(true);
      }}
      onDragLeave={() => setOver(false)}
      onDrop={(e) => {
        e.preventDefault();
        e.stopPropagation();
        setOver(false);
        if (!disabled) pick(e.dataTransfer.files);
      }}
      onClick={() => !disabled && !file && inputRef.current?.click()}
      onKeyDown={(e) => {
        if ((e.key === 'Enter' || e.key === ' ') && !disabled && !file) inputRef.current?.click();
      }}
      role="button"
      tabIndex={disabled ? -1 : 0}
      aria-label={`${slot.label}: ${file ? file.name : 'choose or drop a file'}`}
      className={`relative rounded-xl border-2 border-dashed p-4 transition-all ${
        disabled ? 'opacity-60 cursor-not-allowed' : file ? 'cursor-default' : 'cursor-pointer'
      } ${
        over
          ? 'border-cyan-400 bg-cyan-500/10'
          : error
          ? 'border-rose-500/60 bg-rose-500/5'
          : file
          ? 'border-emerald-500/50 bg-emerald-500/5'
          : 'border-slate-700 bg-slate-950/60 hover:border-slate-500'
      }`}
    >
      <input
        ref={inputRef}
        type="file"
        accept={ACCEPT}
        className="hidden"
        onChange={(e) => {
          pick(e.target.files);
          e.target.value = '';
        }}
      />
      <div className="flex items-start gap-3">
        <div
          className={`w-9 h-9 rounded-lg flex items-center justify-center shrink-0 ${
            file ? 'bg-emerald-500/15 text-emerald-300' : 'bg-slate-800 text-cyan-300'
          }`}
        >
          {file ? <FileCheck2 className="w-5 h-5" /> : <Icon className="w-5 h-5" />}
        </div>
        <div className="min-w-0 flex-1">
          <div className="flex items-center gap-2">
            <span className="text-sm font-semibold text-white">{slot.label}</span>
            {slot.required && (
              <span className="text-[10px] px-1.5 py-0.5 rounded bg-cyan-500/15 text-cyan-300 border border-cyan-500/30">
                required
              </span>
            )}
          </div>
          {file ? (
            <div className="mt-1 text-xs text-slate-300 font-mono truncate" title={file.name}>
              {file.name} <span className="text-slate-500">• {formatBytes(file.size)}</span>
            </div>
          ) : (
            <div className="mt-1 text-xs text-slate-400">
              {slot.hint}
              <div className="mt-1 flex items-center gap-1 text-slate-500">
                <FileUp className="w-3.5 h-3.5" /> Click, drop a file here, or paste with Ctrl+V
              </div>
            </div>
          )}
          {error && <div className="mt-1 text-xs text-rose-300">{error}</div>}
        </div>
        {file && !disabled && (
          <button
            type="button"
            onClick={(e) => {
              e.stopPropagation();
              onClear(slot.key);
            }}
            className="p-1 rounded-md text-slate-400 hover:text-white hover:bg-slate-800"
            aria-label={`Remove ${slot.label}`}
          >
            <X className="w-4 h-4" />
          </button>
        )}
      </div>
    </div>
  );
}

/**
 * Three VCF slots. Supports the file picker, drag-and-drop (per slot, or several files on
 * the whole area) and clipboard paste (handled by the parent view, see assignFiles).
 */
export default function FileUpload({ files, errors = {}, disabled = false, onChange }) {
  const [overAll, setOverAll] = useState(false);
  return (
    <div
      onDragOver={(e) => {
        e.preventDefault();
        if (!disabled) setOverAll(true);
      }}
      onDragLeave={(e) => {
        if (!e.currentTarget.contains(e.relatedTarget)) setOverAll(false);
      }}
      onDrop={(e) => {
        e.preventDefault();
        setOverAll(false);
        if (!disabled && e.dataTransfer.files?.length) onChange(assignFiles(e.dataTransfer.files, files));
      }}
      className={`grid grid-cols-1 lg:grid-cols-3 gap-3 rounded-2xl p-1 transition ${overAll ? 'ring-2 ring-cyan-500/40' : ''}`}
    >
      {SLOTS.map((slot) => (
        <SlotDropzone
          key={slot.key}
          slot={slot}
          file={files[slot.key]}
          error={errors[slot.key]}
          disabled={disabled}
          onFile={(key, f) => onChange({ ...files, [key]: f })}
          onClear={(key) => {
            const next = { ...files };
            delete next[key];
            onChange(next);
          }}
        />
      ))}
    </div>
  );
}
