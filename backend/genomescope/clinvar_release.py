"""
Local copy of a ClinVar release, used to scan an analysis' gVCF completely and
allele-specifically (see clinvar_scan.py).

The release is public reference data shared by every analysis. It is cached outside the
repository and outside the per-job temp folders, in GENOMESCOPE_CLINVAR_DIR (default:
%LOCALAPPDATA%/genomescope/clinvar on Windows, $XDG_CACHE_HOME/genomescope/clinvar or
~/.cache/genomescope/clinvar elsewhere).

Source: NCBI's monthly `variant_summary.txt.gz` (one line per variant and assembly, with
VCF-style alleles for GRCh38, RCV accessions, germline classification and review status).

  download_release()  fetch over verified TLS, check the MD5 NCBI publishes, keep a
                      sidecar JSON with the release timestamp
  build_index()       stream the GRCh38 rows into a compact SQLite index
  ClinVarIndex        read-only queries used by the scan

Rows without a precise VCF allele (large CNVs, repeat expansions, translocations ...)
cannot be matched against SNV / indel calls; they are counted in the index metadata and
surfaced in the scan report instead of being dropped silently.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sqlite3
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import httpx

from genomescope.alleles import is_plain_allele, normalize_vcf_allele
from genomescope.indexer import normalize_classification
from genomescope.ncbi_client import review_status_to_stars
from genomescope.vcf_engine import normalize_chrom


def default_clinvar_dir() -> Path:
    env = (os.environ.get("GENOMESCOPE_CLINVAR_DIR") or "").strip()
    if env:
        return Path(env).expanduser()
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA") or (Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache"))
    return base / "genomescope" / "clinvar"


CLINVAR_DIR = default_clinvar_dir()
SUMMARY_NAME = "variant_summary.txt.gz"
SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
SUMMARY_MD5_URL = SUMMARY_URL + ".md5"
INDEX_NAME = "clinvar_grch38.db"
SCHEMA_VERSION = 2

# Interpretation tiers (higher = more clinically significant). Tiers 0-2 (no assertion for
# the single variant, benign, uncertain significance) are stored as allele keys only (for
# adjudication, linking and counting); full annotation is stored from tier 3 up.
# Annotation for a VUS can be read on demand with annotations_from_summary().
TIER_NONE, TIER_BENIGN, TIER_VUS, TIER_CONFLICTING, TIER_OTHER, TIER_LIKELY_PATH, TIER_PATH = 0, 1, 2, 3, 4, 5, 6
TIER_NAMES = {
    TIER_NONE: "no clinical assertion",
    TIER_BENIGN: "benign / likely benign",
    TIER_VUS: "uncertain significance",
    TIER_CONFLICTING: "conflicting interpretations",
    TIER_OTHER: "risk factor / drug response / association / protective",
    TIER_LIKELY_PATH: "likely pathogenic",
    TIER_PATH: "pathogenic",
}
_NO_ASSERTION = {
    "", "-", "not provided", "other", "no classification provided",
    "no classification for the single variant", "no classifications from unflagged records",
}
_PLACEHOLDER_CONDITIONS = {"not provided", "not specified", "none provided", "see cases", "-", ""}
REQUIRED_COLUMNS = (
    "VariationID", "Type", "Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF",
    "AlternateAlleleVCF", "ClinicalSignificance", "ReviewStatus", "RCVaccession",
)
OPTIONAL_COLUMNS = (
    "#AlleleID", "Name", "GeneID", "GeneSymbol", "LastEvaluated", "RS# (dbSNP)", "PhenotypeList",
    "Origin", "NumberSubmitters",
)
_VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


class ClinVarError(RuntimeError):
    """Raised for problems with the ClinVar release (bad download, unexpected format)."""


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def classify_significance(sig: str) -> Tuple[int, str]:
    """ClinVar germline classification string -> (tier, display label)."""
    raw = (sig or "").strip()
    if raw.lower() in _NO_ASSERTION:
        return (TIER_NONE, "")
    label = normalize_classification(raw)
    low = label.lower()
    if low.startswith("pathogenic"):
        return (TIER_PATH, label)
    if low.startswith("likely pathogenic"):
        return (TIER_LIKELY_PATH, label)
    if low.startswith("conflicting"):
        return (TIER_CONFLICTING, label)
    if low in ("risk factor", "drug response", "association", "protective", "likely risk allele"):
        return (TIER_OTHER, label)
    if low.startswith("uncertain") or low.startswith("vus"):
        return (TIER_VUS, "Uncertain Significance")
    if "benign" in low:
        return (TIER_BENIGN, label)
    return (TIER_VUS, label)  # unknown term: keep it visible rather than hide it


def parse_clinvar_date(raw: str) -> Optional[str]:
    """'Sep 05, 2013' -> '2013-09-05'; '-' / unparseable -> None."""
    raw = (raw or "").strip()
    if not raw or raw == "-":
        return None
    for fmt in ("%b %d, %Y", "%Y-%m-%d", "%Y/%m/%d"):
        try:
            return datetime.strptime(raw, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    return None


def clean_conditions(phenotype_list: str) -> str:
    seen: List[str] = []
    for part in (phenotype_list or "").split("|"):
        p = part.strip()
        if p.lower() in _PLACEHOLDER_CONDITIONS or p in seen:
            continue
        seen.append(p)
    return " | ".join(seen)


def md5_of(path: Path, chunk: int = 4 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def _sidecar(dest_dir: Path) -> Path:
    return dest_dir / (SUMMARY_NAME + ".meta.json")


def read_release_info(dest_dir: Path = CLINVAR_DIR) -> Dict[str, Any]:
    try:
        return json.loads(_sidecar(dest_dir).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def download_release(
    dest_dir: Path = CLINVAR_DIR,
    force: bool = False,
    progress: Optional[Callable[[str], None]] = None,
    client: Optional[httpx.Client] = None,
) -> Dict[str, Any]:
    """
    Download variant_summary.txt.gz (verified TLS) and check it against the MD5 NCBI
    publishes next to it. An already-downloaded file whose MD5 equals the published one is
    kept. Returns the sidecar info {url, md5, size, last_modified, downloaded_at}.
    """
    say = progress or (lambda _m: None)
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / SUMMARY_NAME
    own = client is None
    http = client or httpx.Client(timeout=httpx.Timeout(60.0, connect=15.0), follow_redirects=True)
    try:
        r = http.get(SUMMARY_MD5_URL)
        r.raise_for_status()
        tokens = r.text.split()
        expected = tokens[0].lower() if tokens else ""
        if len(expected) != 32:
            raise ClinVarError(f"unexpected MD5 file from NCBI: {r.text[:80]!r}")
        if target.exists() and not force and md5_of(target) == expected:
            info = read_release_info(dest_dir)
            if info.get("md5") != expected:
                # verified file without a download record (e.g. fetched by hand): record it
                head = http.head(SUMMARY_URL)
                info = {
                    "url": SUMMARY_URL, "md5": expected, "size": target.stat().st_size,
                    "last_modified": head.headers.get("Last-Modified", ""),
                    "downloaded_at": datetime.fromtimestamp(target.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"),
                }
                _sidecar(dest_dir).write_text(json.dumps(info, indent=2), encoding="utf-8")
            say(f"ClinVar release already present and verified (md5 {expected}).")
            return info
        part = target.with_suffix(target.suffix + ".part")
        h = hashlib.md5()
        size = 0
        with http.stream("GET", SUMMARY_URL) as resp:
            resp.raise_for_status()
            last_modified = resp.headers.get("Last-Modified", "")
            total = int(resp.headers.get("Content-Length", "0") or 0)
            next_report = 0
            with open(part, "wb") as out:
                for block in resp.iter_bytes(1024 * 1024):
                    out.write(block)
                    h.update(block)
                    size += len(block)
                    if size >= next_report:
                        say(f"  downloaded {size / 1e6:,.0f} / {total / 1e6:,.0f} MB")
                        next_report = size + 50_000_000
        if h.hexdigest() != expected:
            part.unlink(missing_ok=True)
            raise ClinVarError(f"MD5 mismatch for {SUMMARY_NAME}: got {h.hexdigest()}, NCBI publishes {expected}")
        part.replace(target)
        info = {
            "url": SUMMARY_URL,
            "md5": expected,
            "size": size,
            "last_modified": last_modified,
            "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        _sidecar(dest_dir).write_text(json.dumps(info, indent=2), encoding="utf-8")
        say(f"ClinVar release downloaded and verified (md5 {expected}).")
        return info
    finally:
        if own:
            http.close()


def release_label(info: Dict[str, Any], summary_path: Path) -> str:
    """ISO date of the release: the server's Last-Modified, else the file's mtime."""
    lm = (info or {}).get("last_modified") or ""
    for fmt in ("%a, %d %b %Y %H:%M:%S %Z", "%a, %d %b %Y %H:%M:%S GMT"):
        try:
            return datetime.strptime(lm, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    try:
        return datetime.fromtimestamp(summary_path.stat().st_mtime, timezone.utc).strftime("%Y-%m-%d")
    except OSError:
        return "unknown"


# ---------------------------------------------------------------------------
# Index build
# ---------------------------------------------------------------------------

_SCHEMA = """
CREATE TABLE alleles (
    variation_id INTEGER NOT NULL,
    chr TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    npos INTEGER NOT NULL,
    nref TEXT NOT NULL,
    nalt TEXT NOT NULL,
    indel INTEGER NOT NULL,
    tier INTEGER NOT NULL,
    stars INTEGER NOT NULL,
    rsid INTEGER
);
CREATE TABLE annot (
    variation_id INTEGER PRIMARY KEY,
    allele_id INTEGER,
    type TEXT,
    name TEXT,
    gene_symbol TEXT,
    gene_id TEXT,
    rsid TEXT,
    significance TEXT,
    classification TEXT,
    review_status TEXT,
    last_evaluated TEXT,
    rcvs TEXT,
    conditions TEXT,
    origin TEXT,
    n_submitters INTEGER
);
CREATE TABLE rcv (rcv INTEGER NOT NULL, variation_id INTEGER NOT NULL);  -- accession number without the 'RCV' prefix
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
"""


def _open_gz(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")


def build_index(
    summary_gz: Path,
    db_path: Path,
    info: Optional[Dict[str, Any]] = None,
    progress: Optional[Callable[[str], None]] = None,
) -> Dict[str, Any]:
    """
    Stream the GRCh38 rows of `summary_gz` into a fresh SQLite index at `db_path`
    (replaced atomically). Returns the statistics stored in the `meta` table.
    """
    say = progress or (lambda _m: None)
    info = info if info is not None else read_release_info(summary_gz.parent)
    tmp = db_path.with_name(db_path.name + ".building")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(str(tmp))
    t0 = time.time()
    tiers: Counter = Counter()
    skipped_types: Counter = Counter()
    rows_total = rows_grch38 = rows_vcf = 0
    try:
        conn.executescript(_SCHEMA)
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        with _open_gz(summary_gz) as f:
            header = f.readline().rstrip("\r\n").split("\t")
            col = {h: i for i, h in enumerate(header)}
            missing = [c for c in REQUIRED_COLUMNS if c not in col]
            if missing:
                raise ClinVarError(
                    f"{summary_gz.name} is missing expected column(s) {missing}; ClinVar may have changed the "
                    f"file layout. Columns found: {header}"
                )
            ix = {c: col[c] for c in REQUIRED_COLUMNS}
            opt = {c: col.get(c, -1) for c in OPTIONAL_COLUMNS}

            def get(parts: Sequence[str], name: str) -> str:
                i = opt[name]
                return parts[i] if 0 <= i < len(parts) else ""

            allele_rows: List[tuple] = []
            annot_rows: List[tuple] = []
            rcv_rows: List[tuple] = []

            def flush() -> None:
                conn.executemany("INSERT INTO alleles VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", allele_rows)
                conn.executemany("INSERT OR REPLACE INTO annot VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", annot_rows)
                conn.executemany("INSERT INTO rcv VALUES (?,?)", rcv_rows)
                allele_rows.clear(); annot_rows.clear(); rcv_rows.clear()

            need = max(ix.values()) + 1
            for line in f:
                rows_total += 1
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < need or parts[ix["Assembly"]] != "GRCh38":
                    continue
                rows_grch38 += 1
                chrom = normalize_chrom(parts[ix["Chromosome"]])
                pos_s, ref, alt = parts[ix["PositionVCF"]], parts[ix["ReferenceAlleleVCF"]], parts[ix["AlternateAlleleVCF"]]
                if chrom not in _VALID_CHROMS or not pos_s.isdigit() or not (is_plain_allele(ref) and is_plain_allele(alt)):
                    skipped_types[parts[ix["Type"]]] += 1
                    continue
                tier, label = classify_significance(parts[ix["ClinicalSignificance"]])
                rows_vcf += 1
                tiers[tier] += 1
                vid = int(parts[ix["VariationID"]])
                pos, ref, alt = int(pos_s), ref.upper(), alt.upper()
                npos, nref, nalt = normalize_vcf_allele(pos, ref, alt)
                review = parts[ix["ReviewStatus"]]
                stars = review_status_to_stars(review)
                rs_s = get(parts, "RS# (dbSNP)").strip()
                allele_rows.append((vid, chrom, pos, ref, alt, npos, nref, nalt, 1 if len(nref) != len(nalt) else 0, tier, stars,
                                    int(rs_s) if rs_s.isdigit() and int(rs_s) > 0 else None))
                rcvs = [r for r in parts[ix["RCVaccession"]].split("|") if r.startswith("RCV")]
                rcv_rows.extend((int(r[3:]), vid) for r in rcvs if r[3:].isdigit())
                if tier >= TIER_CONFLICTING:
                    rs = get(parts, "RS# (dbSNP)").strip()
                    n_sub = get(parts, "NumberSubmitters")
                    annot_rows.append((
                        vid,
                        int(get(parts, "#AlleleID")) if get(parts, "#AlleleID").isdigit() else None,
                        parts[ix["Type"]],
                        get(parts, "Name"),
                        get(parts, "GeneSymbol"),
                        get(parts, "GeneID"),
                        f"rs{rs}" if rs.isdigit() and int(rs) > 0 else None,
                        parts[ix["ClinicalSignificance"]],
                        label,
                        review,
                        parse_clinvar_date(get(parts, "LastEvaluated")),
                        ", ".join(rcvs),
                        clean_conditions(get(parts, "PhenotypeList")),
                        get(parts, "Origin"),
                        int(n_sub) if n_sub.isdigit() else None,
                    ))
                if len(allele_rows) >= 50_000:
                    flush()
                    if rows_total % 500_000 < 50_000:
                        say(f"  indexed {rows_total:,} ClinVar rows ...")
            flush()
        if rows_vcf == 0:
            raise ClinVarError(f"no usable GRCh38 SNV/indel rows found in {summary_gz.name}")
        say("  building indexes ...")
        conn.executescript(
            """
            CREATE INDEX idx_alleles_key ON alleles(chr, npos, nref, nalt);
            CREATE INDEX idx_alleles_indel ON alleles(chr, pos) WHERE indel = 1;
            CREATE INDEX idx_alleles_var ON alleles(variation_id);
            CREATE INDEX idx_rcv ON rcv(rcv);
            CREATE INDEX idx_alleles_rsid ON alleles(rsid) WHERE rsid IS NOT NULL;
            CREATE INDEX idx_annot_rsid ON annot(rsid);
            """
        )
        stats = {
            "schema_version": SCHEMA_VERSION,
            "source_file": summary_gz.name,
            "source_md5": (info or {}).get("md5"),
            "source_size": summary_gz.stat().st_size,
            "release_date": release_label(info or {}, summary_gz),
            "release_date_source": "NCBI Last-Modified header" if (info or {}).get("last_modified") else "local file modification time (no download record)",
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "rows_total": rows_total,
            "rows_grch38": rows_grch38,
            "rows_indexed": rows_vcf,
            "rows_no_clinical_assertion": tiers[TIER_NONE],
            "rows_not_vcf_representable": sum(skipped_types.values()),
            "not_vcf_representable_by_type": dict(skipped_types.most_common()),
            "indexed_by_tier": {TIER_NAMES[t]: tiers[t] for t in sorted(TIER_NAMES)},
            "build_seconds": round(time.time() - t0, 1),
        }
        conn.executemany("INSERT INTO meta VALUES (?,?)", [(k, json.dumps(v)) for k, v in stats.items()])
        conn.commit()
    except BaseException:
        conn.close()
        tmp.unlink(missing_ok=True)
        raise
    conn.close()
    tmp.replace(db_path)
    return stats


def annotations_from_summary(summary_gz: Path, variation_ids: Iterable[int]) -> Dict[int, Dict[str, Any]]:
    """Annotation dicts (same keys as the `annot` table) for specific GRCh38 variation ids,
    read straight from the release file. Used for tiers whose annotation is not indexed."""
    want = {int(v) for v in variation_ids}
    out: Dict[int, Dict[str, Any]] = {}
    if not want:
        return out
    with _open_gz(summary_gz) as f:
        header = f.readline().rstrip("\r\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        if any(c not in col for c in REQUIRED_COLUMNS):
            raise ClinVarError(f"{summary_gz.name} does not have the expected columns")
        vi, asm = col["VariationID"], col["Assembly"]

        def g(parts: Sequence[str], name: str) -> str:
            i = col.get(name, -1)
            return parts[i] if 0 <= i < len(parts) else ""

        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= max(vi, asm) or parts[asm] != "GRCh38" or not parts[vi].isdigit() or int(parts[vi]) not in want:
                continue
            vid = int(parts[vi])
            tier, label = classify_significance(g(parts, "ClinicalSignificance"))
            rs = g(parts, "RS# (dbSNP)").strip()
            rcvs = [r for r in g(parts, "RCVaccession").split("|") if r.startswith("RCV")]
            out[vid] = {
                "variation_id": vid, "type": g(parts, "Type"), "name": g(parts, "Name"),
                "gene_symbol": g(parts, "GeneSymbol"), "gene_id": g(parts, "GeneID"),
                "rsid": f"rs{rs}" if rs.isdigit() and int(rs) > 0 else None,
                "significance": g(parts, "ClinicalSignificance"), "classification": label,
                "review_status": g(parts, "ReviewStatus"), "last_evaluated": parse_clinvar_date(g(parts, "LastEvaluated")),
                "rcvs": ", ".join(rcvs), "conditions": clean_conditions(g(parts, "PhenotypeList")), "origin": g(parts, "Origin"),
            }
            if len(out) == len(want):
                break
    return out


# ---------------------------------------------------------------------------
# Read-only access
# ---------------------------------------------------------------------------

INDEX_PATH = CLINVAR_DIR / INDEX_NAME


class ClinVarIndex:
    """Read-only view of the local ClinVar index."""

    def __init__(self, path: Path = INDEX_PATH):
        self.path = Path(path)
        if not self.path.exists():
            raise ClinVarError(
                f"ClinVar index not found at {self.path}. Build it with: python -m genomescope.clinvar_scan --workspace <job_dir> --update (or tick 'download ClinVar release' when starting an analysis)"
            )
        uri = self.path.resolve().as_uri() + "?mode=ro"
        self.conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
        self.conn.execute("PRAGMA query_only=ON")
        self.conn.execute("PRAGMA mmap_size=1073741824")
        self._meta: Optional[Dict[str, Any]] = None
        v = self.meta().get("schema_version")
        if v != SCHEMA_VERSION:
            raise ClinVarError(f"ClinVar index schema {v!r} != {SCHEMA_VERSION}; rebuild with --update")

    @staticmethod
    def available(path: Path = INDEX_PATH) -> bool:
        return Path(path).exists()

    def close(self) -> None:
        self.conn.close()

    def meta(self) -> Dict[str, Any]:
        if self._meta is None:
            self._meta = {k: json.loads(v) for k, v in self.conn.execute("SELECT key, value FROM meta")}
        return self._meta

    def exact(self, chrom: str, npos: int, nref: str, nalt: str) -> List[Tuple[int, int, str, str, int, int]]:
        """Alleles with this normalised (chr, pos, del, ins): (variation_id, pos, ref, alt, tier, stars)."""
        return self.conn.execute(
            "SELECT variation_id, pos, ref, alt, tier, stars FROM alleles WHERE chr = ? AND npos = ? AND nref = ? AND nalt = ?",
            (chrom, npos, nref, nalt),
        ).fetchall()

    def indels_between(self, chrom: str, lo: int, hi: int) -> List[Tuple[int, int, str, str, int, int]]:
        return self.conn.execute(
            "SELECT variation_id, pos, ref, alt, tier, stars FROM alleles WHERE indel = 1 AND chr = ? AND pos BETWEEN ? AND ?",
            (chrom, lo, hi),
        ).fetchall()

    def allele(self, variation_id: int) -> Optional[Dict[str, Any]]:
        rows = self.alleles_for([variation_id])
        return rows.get(int(variation_id))

    def alleles_for(self, variation_ids: Iterable[int]) -> Dict[int, Dict[str, Any]]:
        ids = sorted({int(v) for v in variation_ids})
        out: Dict[int, Dict[str, Any]] = {}
        for i in range(0, len(ids), 500):
            chunk = ids[i:i + 500]
            q = ",".join("?" * len(chunk))
            for r in self.conn.execute(
                f"SELECT variation_id, chr, pos, ref, alt, tier, stars FROM alleles WHERE variation_id IN ({q})", chunk
            ):
                out[r[0]] = {"variation_id": r[0], "chr": r[1], "pos": r[2], "ref": r[3], "alt": r[4], "tier": r[5], "stars": r[6]}
        return out

    def alleles_by_rsid(self, rsids: Iterable[str]) -> Dict[str, List[Dict[str, Any]]]:
        """All indexed alleles of the given dbSNP rsIDs (strings like 'rs671')."""
        nums = sorted({int(r[2:]) for r in rsids if r and r.lower().startswith("rs") and r[2:].isdigit()})
        out: Dict[str, List[Dict[str, Any]]] = {}
        for i in range(0, len(nums), 500):
            chunk = nums[i:i + 500]
            q = ",".join("?" * len(chunk))
            for r in self.conn.execute(
                f"SELECT variation_id, chr, pos, ref, alt, tier, stars, rsid FROM alleles WHERE rsid IN ({q})", chunk
            ):
                out.setdefault(f"rs{r[7]}", []).append(
                    {"variation_id": r[0], "chr": r[1], "pos": r[2], "ref": r[3], "alt": r[4], "tier": r[5], "stars": r[6]}
                )
        return out

    def annotations_for(self, variation_ids: Iterable[int]) -> Dict[int, Dict[str, Any]]:
        ids = sorted({int(v) for v in variation_ids})
        out: Dict[int, Dict[str, Any]] = {}
        for i in range(0, len(ids), 500):
            chunk = ids[i:i + 500]
            q = ",".join("?" * len(chunk))
            cur = self.conn.execute(f"SELECT * FROM annot WHERE variation_id IN ({q})", chunk)
            names = [d[0] for d in cur.description]
            for r in cur:
                d = dict(zip(names, r))
                out[d["variation_id"]] = d
        return out

    def variations_for_rcvs(self, rcvs: Iterable[str]) -> Dict[str, int]:
        nums = sorted({int(r.strip()[3:]) for r in rcvs if r and r.strip().startswith("RCV") and r.strip()[3:].isdigit()})
        out: Dict[str, int] = {}
        for i in range(0, len(nums), 500):
            chunk = nums[i:i + 500]
            q = ",".join("?" * len(chunk))
            for r in self.conn.execute(f"SELECT rcv, variation_id FROM rcv WHERE rcv IN ({q})", chunk):
                out[f"RCV{r[0]:09d}"] = r[1]
        return out

    def high_impact_loci(self) -> Dict[str, List[int]]:
        """Sorted distinct positions of Pathogenic / Likely pathogenic alleles, per chromosome."""
        out: Dict[str, List[int]] = {}
        for chrom, pos in self.conn.execute(
            "SELECT DISTINCT chr, pos FROM alleles WHERE tier >= ? ORDER BY chr, pos", (TIER_LIKELY_PATH,)
        ):
            out.setdefault(chrom, []).append(pos)
        return out
