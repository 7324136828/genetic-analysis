const API_BASE = import.meta.env.VITE_API_BASE || '/api';

// State-changing backend routes require this header (blocks cross-site form/img requests).
const APP_HEADERS = { 'X-Requested-With': 'GenomeScope' };
// Every explorer request names the analysis it works on (the backend activates that
// analysis' isolated temp workspace for the request).
const JOB_HEADER = 'X-GenomeScope-Job';

let activeJobId = null;

/** Select the analysis all explorer / NCBI calls operate on (null = none). */
export function setActiveJob(jobId) {
  activeJobId = jobId || null;
}

export function getActiveJob() {
  return activeJobId;
}

function apiFetch(url, options = {}) {
  const headers = { ...(options.headers || {}) };
  if (activeJobId && !headers[JOB_HEADER]) headers[JOB_HEADER] = activeJobId;
  return fetch(url, { ...options, headers });
}

export async function fetchSummary() {
  const res = await apiFetch(`${API_BASE}/summary`);
  if (!res.ok) throw new Error(`Failed to fetch summary: ${res.status}`);
  return res.json();
}

export async function fetchFilterOptions() {
  const res = await apiFetch(`${API_BASE}/filters`);
  if (!res.ok) throw new Error(`Failed to fetch filters: ${res.status}`);
  return res.json();
}

export async function fetchVariants(params = {}) {
  const searchParams = new URLSearchParams();
  Object.entries(params).forEach(([k, v]) => {
    if (v !== undefined && v !== null && String(v).trim() !== '') {
      searchParams.set(k, String(v).trim());
    }
  });
  const res = await apiFetch(`${API_BASE}/variants?${searchParams.toString()}`);
  if (!res.ok) throw new Error(`Failed to fetch variants: ${res.status}`);
  return res.json();
}

export async function fetchVariantDetail(id) {
  const res = await apiFetch(`${API_BASE}/variants/${id}`);
  if (!res.ok) throw new Error(`Failed to fetch variant detail: ${res.status}`);
  return res.json();
}

export async function fetchDirectVcfReads({ chr, position, variants_only = false, limit = 50 }) {
  const searchParams = new URLSearchParams({
    chr: String(chr),
    position: String(position),
    variants_only: String(variants_only),
    limit: String(limit),
  });
  const res = await apiFetch(`${API_BASE}/vcf/read?${searchParams.toString()}`);
  if (!res.ok) {
    const err = await res.json().catch(() => ({}));
    throw new Error(err.detail || `Failed to read VCF: ${res.status}`);
  }
  return res.json();
}

// ---------------------------------------------------------------------------
// NCBI (ClinVar / dbSNP / Gene / PubMed) live endpoints
// ---------------------------------------------------------------------------

function describeDetail(detail) {
  if (!detail) return null;
  if (typeof detail === 'string') return detail;
  if (Array.isArray(detail)) {
    // FastAPI validation errors: [{loc, msg, type}]
    return detail
      .map((d) => (d && typeof d === 'object' ? `${(d.loc || []).slice(1).join('.') || 'request'}: ${d.msg}` : String(d)))
      .join('; ');
  }
  if (typeof detail === 'object') return detail.message || JSON.stringify(detail);
  return String(detail);
}

async function requestJson(url, options = {}, label = 'Request') {
  let res;
  try {
    res = await apiFetch(url, options);
  } catch (networkErr) {
    throw new Error(`${label} failed: backend unreachable (${networkErr.message || 'network error'})`);
  }
  if (!res.ok) {
    const body = await res.json().catch(() => ({}));
    const msg = describeDetail(body.detail) || `${label} failed: HTTP ${res.status}`;
    const err = new Error(msg);
    err.status = res.status;
    throw err;
  }
  return res.json();
}

function buildQuery(params = {}) {
  const sp = new URLSearchParams();
  Object.entries(params).forEach(([k, v]) => {
    if (v !== undefined && v !== null && String(v).trim() !== '') {
      sp.set(k, String(v).trim());
    }
  });
  const qs = sp.toString();
  return qs ? `?${qs}` : '';
}

export function fetchClinvarStatus() {
  return requestJson(`${API_BASE}/clinvar/status`, {}, 'ClinVar scan status');
}

export function fetchNcbiStatus() {
  return requestJson(`${API_BASE}/ncbi/status`, {}, 'NCBI status');
}

/** params: { rsid, rcv_id, vcv_id, chr, position, include_publications } */
export function fetchNcbiVariant(params = {}) {
  return requestJson(`${API_BASE}/ncbi/variant${buildQuery(params)}`, {}, 'NCBI variant lookup');
}

export function fetchNcbiGene(symbol, { significance = 'pathogenic', limit = 40 } = {}) {
  const sym = encodeURIComponent(String(symbol || '').trim());
  return requestJson(
    `${API_BASE}/ncbi/gene/${sym}${buildQuery({ significance, limit })}`,
    {},
    'NCBI gene scan'
  );
}

export function searchNcbi(term, limit = 25) {
  return requestJson(`${API_BASE}/ncbi/search${buildQuery({ term, limit })}`, {}, 'NCBI ClinVar search');
}

export function importNcbiVariants(variationIds = []) {
  return requestJson(
    `${API_BASE}/ncbi/import`,
    {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', ...APP_HEADERS },
      body: JSON.stringify({ variation_ids: variationIds.map(String) }),
    },
    'NCBI import'
  );
}

/** Read-only NCBI comparison for a catalog row; { sync: true } overwrites the row (POST). */
export function fetchVariantNcbi(recordId, { sync = false } = {}) {
  const base = `${API_BASE}/variants/${encodeURIComponent(recordId)}/ncbi`;
  if (sync) {
    return requestJson(`${base}/sync`, { method: 'POST', headers: APP_HEADERS }, 'NCBI catalog sync');
  }
  return requestJson(base, {}, 'NCBI comparison');
}

// ---------------------------------------------------------------------------
// Analysis jobs (upload, progress, history, download)
// ---------------------------------------------------------------------------

export function fetchConfig() {
  return requestJson(`${API_BASE}/config`, {}, 'Configuration');
}

export function listJobs() {
  return requestJson(`${API_BASE}/jobs`, {}, 'Analysis history');
}

export function getJob(jobId) {
  return requestJson(`${API_BASE}/jobs/${encodeURIComponent(jobId)}`, {}, 'Analysis status');
}

export function getJobLogs(jobId, tail = 300) {
  return requestJson(`${API_BASE}/jobs/${encodeURIComponent(jobId)}/logs${buildQuery({ tail })}`, {}, 'Analysis log');
}

export function discardJob(jobId) {
  return requestJson(
    `${API_BASE}/jobs/${encodeURIComponent(jobId)}/discard`,
    { method: 'POST', headers: APP_HEADERS },
    'Discard analysis'
  );
}

export async function deleteJob(jobId) {
  let res;
  try {
    res = await fetch(`${API_BASE}/jobs/${encodeURIComponent(jobId)}`, { method: 'DELETE', headers: APP_HEADERS });
  } catch (networkErr) {
    throw new Error(`Delete failed: backend unreachable (${networkErr.message || 'network error'})`);
  }
  if (!res.ok) {
    const body = await res.json().catch(() => ({}));
    throw new Error(describeDetail(body.detail) || `Delete failed: HTTP ${res.status}`);
  }
  return true;
}

export function downloadZipUrl(jobId) {
  return `${API_BASE}/jobs/${encodeURIComponent(jobId)}/download-zip`;
}

/**
 * Upload VCFs and start an analysis. Uses XMLHttpRequest for upload progress (fetch has none).
 * files: { snp_indel: File, sv?: File, cnv?: File }; options: { label, clinvar_mode }
 * Returns { promise, abort }.
 */
export function createJob(files, { label = '', clinvar_mode = 'cached' } = {}, onProgress = () => {}) {
  const form = new FormData();
  Object.entries(files).forEach(([slot, file]) => {
    if (file) form.append(slot, file, file.name);
  });
  form.append('label', label);
  form.append('clinvar_mode', clinvar_mode);
  const xhr = new XMLHttpRequest();
  const promise = new Promise((resolve, reject) => {
    xhr.open('POST', `${API_BASE}/jobs`);
    xhr.setRequestHeader('X-Requested-With', APP_HEADERS['X-Requested-With']);
    xhr.responseType = 'json';
    xhr.upload.onprogress = (e) => {
      if (e.lengthComputable) onProgress(Math.round((e.loaded / e.total) * 100));
    };
    xhr.onload = () => {
      const body = xhr.response || {};
      if (xhr.status >= 200 && xhr.status < 300) resolve(body);
      else reject(new Error(describeDetail(body.detail) || `Upload failed: HTTP ${xhr.status}`));
    };
    xhr.onerror = () => reject(new Error('Upload failed: backend unreachable'));
    xhr.onabort = () => reject(new Error('Upload cancelled'));
    xhr.send(form);
  });
  return { promise, abort: () => xhr.abort() };
}
