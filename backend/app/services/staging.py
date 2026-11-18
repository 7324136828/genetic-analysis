"""
Staging of uploaded VCF files into a job's ``inputs/`` folder.

- Accepts plain ``.vcf`` or gzip / bgzip compressed ``.vcf.gz`` (decompressed on the fly:
  the gVCF reader needs byte offsets into an uncompressed file).
- Files are stored under canonical names (genomescope.workspace.SNP_VCF_NAME, ...), so a
  user-supplied file name never reaches the catalog, reports or logs. The original name is
  kept only in ``inputs/manifest.json`` (shown in the History screen).
- The header is validated: VCF signature, a ``#CHROM`` line with a sample column, and the
  reference assembly (GRCh37 / hg19 inputs are rejected; annotations are GRCh38).
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
from pathlib import Path
from typing import Any, BinaryIO, Dict, List, Optional

from genomescope.workspace import CNV_VCF_NAME, SNP_VCF_NAME, SV_VCF_NAME, Workspace

SLOTS: Dict[str, Dict[str, str]] = {
    "snp_indel": {"canonical_name": SNP_VCF_NAME, "label": "SNP & Indel gVCF"},
    "sv": {"canonical_name": SV_VCF_NAME, "label": "Structural Variants VCF"},
    "cnv": {"canonical_name": CNV_VCF_NAME, "label": "Copy Number Variants VCF"},
}
REQUIRED_SLOTS = ("snp_indel",)

GZIP_MAGIC = b"\x1f\x8b"
CHUNK = 4 * 1024 * 1024
# Zip-bomb guard: a decompressed VCF larger than this is rejected.
MAX_DECOMPRESSED_BYTES = 64 * 1024 ** 3
MAX_HEADER_LINES = 200_000

# chr1 length identifies the assembly reliably when ##reference is vague.
_CHR1_LENGTHS = {248956422: "GRCh38", 249250621: "GRCh37"}
_GRCH37_HINTS = ("grch37", "hg19", "b37", "hs37", "human_g1k_v37", "ncbi37")
_GRCH38_HINTS = ("grch38", "hg38", "hs38")


class StagingError(ValueError):
    """The upload is not a usable VCF; the message is safe to show to the user."""


def safe_display_name(name: Optional[str]) -> str:
    base = Path(str(name or "")).name.replace("\\", "/").split("/")[-1]
    base = re.sub(r"[\x00-\x1f\x7f]", "", base).strip()
    return base[:200] or "upload.vcf"


def _copy_stream(src: BinaryIO, dest: Path) -> Dict[str, Any]:
    head = src.read(2)
    compressed = head == GZIP_MAGIC
    try:
        src.seek(0)
    except (AttributeError, OSError):
        raise StagingError("upload stream is not seekable")
    reader: BinaryIO = gzip.GzipFile(fileobj=src, mode="rb") if compressed else src  # type: ignore[assignment]
    written = 0
    try:
        with open(dest, "wb") as out:
            while True:
                block = reader.read(CHUNK)
                if not block:
                    break
                written += len(block)
                if written > MAX_DECOMPRESSED_BYTES:
                    raise StagingError("decompressed VCF is larger than 64 GB")
                out.write(block)
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        dest.unlink(missing_ok=True)
        raise StagingError(f"could not read the upload ({'corrupt gzip' if compressed else 'read error'}: {e.__class__.__name__})")
    except StagingError:
        dest.unlink(missing_ok=True)
        raise
    return {"compressed": compressed, "size_bytes": written}


def inspect_header(path: Path) -> Dict[str, Any]:
    """Validate a staged VCF and return {reference, samples, header_lines}."""
    reference: Optional[str] = None
    reference_line = ""
    columns: List[str] = []
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
        if not first.startswith("##fileformat=VCF"):
            raise StagingError("not a VCF file (missing '##fileformat=VCF' on the first line)")
        n = 1
        for line in f:
            n += 1
            if n > MAX_HEADER_LINES:
                raise StagingError("VCF header is too long or the '#CHROM' line is missing")
            if line.startswith("##reference="):
                reference_line = line.strip().lower()
            elif line.startswith("##contig=") and reference is None:
                m = re.search(r"ID=(?:chr)?1,.*?length=(\d+)", line) or re.search(r"length=(\d+).*?ID=(?:chr)?1[,>]", line)
                if m:
                    reference = _CHR1_LENGTHS.get(int(m.group(1)))
            elif line.startswith("#CHROM"):
                columns = line.rstrip("\r\n").split("\t")
                break
            elif not line.startswith("#"):
                break
    if not columns:
        raise StagingError("VCF has no '#CHROM' header line")
    if len(columns) < 10:
        raise StagingError("VCF has no sample (genotype) column; a single-sample VCF is required")
    if reference is None and reference_line:
        if any(h in reference_line for h in _GRCH37_HINTS):
            reference = "GRCh37"
        elif any(h in reference_line for h in _GRCH38_HINTS):
            reference = "GRCh38"
    if reference == "GRCh37":
        raise StagingError("this VCF is aligned to GRCh37/hg19; GenomeScope annotates GRCh38 coordinates only")
    return {"reference": reference or "GRCh38", "reference_detected": reference is not None, "sample_columns": len(columns) - 9, "header_lines": n}


def stage_file(ws: Workspace, slot: str, stream: BinaryIO, original_name: Optional[str]) -> Dict[str, Any]:
    if slot not in SLOTS:
        raise StagingError(f"unknown upload slot {slot!r}")
    canonical = SLOTS[slot]["canonical_name"]
    dest = ws.inputs_dir / canonical
    info = _copy_stream(stream, dest)
    if info["size_bytes"] == 0:
        dest.unlink(missing_ok=True)
        raise StagingError(f"{SLOTS[slot]['label']}: the file is empty")
    try:
        header = inspect_header(dest)
    except StagingError as e:
        dest.unlink(missing_ok=True)
        raise StagingError(f"{SLOTS[slot]['label']}: {e}")
    return {
        "slot": slot,
        "label": SLOTS[slot]["label"],
        "canonical_name": canonical,
        "original_filename": safe_display_name(original_name),
        "size_bytes": info["size_bytes"],
        "compressed_upload": info["compressed"],
        "reference": header["reference"],
        "reference_detected": header["reference_detected"],
        "sample_columns": header["sample_columns"],
    }


def write_manifest(ws: Workspace, files: List[Dict[str, Any]]) -> Dict[str, Any]:
    refs = {f["reference"] for f in files if f.get("reference_detected")}
    manifest = {"reference": (refs.pop() if len(refs) == 1 else "GRCh38"), "files": files}
    ws.manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def clear_inputs(ws: Workspace) -> None:
    if ws.inputs_dir.exists():
        shutil.rmtree(ws.inputs_dir, ignore_errors=True)
    ws.inputs_dir.mkdir(parents=True, exist_ok=True)
