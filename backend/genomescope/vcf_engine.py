import json
import os
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from genomescope import workspace
from genomescope.workspace import CNV_VCF_NAME, SNP_VCF_NAME, SV_VCF_NAME  # noqa: F401  (re-exported)

# Bump when the sparse-index format / sampling changes so cached indexes are rebuilt.
SPARSE_INDEX_VERSION = 2


def normalize_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.upper() == "M":
        return "MT"
    return c.upper() if c.upper() in ("X", "Y", "MT") else c


def parse_info_field(info_str: str) -> Dict[str, Any]:
    if not info_str or info_str == ".":
        return {}
    result: Dict[str, Any] = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            result[k] = v
        elif item:
            result[item] = True
    return result


def parse_sample_field(format_str: str, sample_str: str) -> Dict[str, str]:
    if not format_str or not sample_str:
        return {}
    keys = format_str.split(":")
    vals = sample_str.split(":")
    return {k: (vals[i] if i < len(vals) else ".") for i, k in enumerate(keys)}


def determine_variant_type(ref: str, alt: str, info_map: Dict[str, Any], source_type: str = "snp-indel") -> str:
    if source_type == "cnv":
        svtype = info_map.get("SVTYPE", "CNV")
        if svtype == "LOH":
            return "LOH (Loss of Heterozygosity)"
        if "CN0" in alt and "CN2" not in alt and "CN3" not in alt:
            return "CNV Loss"
        if any(tag in alt for tag in ("<CN2>", "<CN3>", "<CN4>", "<CN5>", "<DUP>")):
            return "CNV Gain"
        if alt == ".":
            return "CNV Reference"
        return "CNV"
    if source_type == "sv":
        svtype = info_map.get("SVTYPE", "")
        if svtype == "DEL":
            return "SV Deletion"
        if svtype == "INS":
            return "SV Insertion"
        if svtype == "DUP":
            return "SV Duplication"
        if svtype == "BND":
            return "SV Breakend (BND)"
        return f"SV ({svtype})" if svtype else "Structural Variant"
    if alt == "." or alt == ref:
        return "SNV"
    first_alt = alt.split(",")[0]
    if info_map.get("INDEL"):
        if len(first_alt) > len(ref):
            return "Insertion"
        if len(first_alt) < len(ref):
            return "Deletion"
        return "Indel"
    if len(ref) == 1 and len(first_alt) == 1:
        return "SNV"
    if len(first_alt) > len(ref):
        return "Insertion"
    if len(first_alt) < len(ref):
        return "Deletion"
    return "Indel"


def format_vcf_genotype_and_status(
    ref: str,
    alt: str,
    qual: str,
    filt: str,
    info_map: Dict[str, Any],
    sample_map: Dict[str, str],
    source_type: str = "snp-indel",
    expected_risk_allele: Optional[str] = None,
) -> Dict[str, Any]:
    gt = sample_map.get("GT", "./.")
    dp_str = sample_map.get("DP") or info_map.get("DP") or "."
    gq_str = sample_map.get("GQ") or sample_map.get("QS") or "."
    min_dp_str = sample_map.get("MIN_DP", ".")
    ad_str = sample_map.get("AD", ".")

    dp_val = int(dp_str) if str(dp_str).isdigit() else None
    gq_val = None
    try:
        if gq_str != ".":
            gq_val = int(float(gq_str))
    except ValueError:
        pass

    qual_val = None
    try:
        if qual != ".":
            qual_val = round(float(qual), 2)
    except ValueError:
        pass

    if source_type == "cnv":
        cn = sample_map.get("CN", ".")
        rc = sample_map.get("RC", ".")
        bc = sample_map.get("BC", ".")
        qs = sample_map.get("QS", qual)
        svtype = info_map.get("SVTYPE", "CNV")
        cnvlen = info_map.get("CNVLEN", ".")
        if alt == ".":
            your_status = "Reference (Not Detected)"
        elif filt != "PASS":
            your_status = "Low Quality Call (Unconfirmed)"
        elif svtype == "LOH":
            your_status = "Loss of Heterozygosity (LOH)"
        elif cn == "0":
            your_status = "Detected (Homozygous)"
        elif cn == "1":
            your_status = "Carrier (Heterozygous)"
        else:
            your_status = "Detected (Copy Number Gain)"
        your_data = f"CN={cn} ({gt}) | RC: {rc} | Bins: {bc} | QS: {qs} | {filt}"
        risk_version = alt if alt != "." else ref
        return {
            "gt": gt,
            "alleles_called": f"CN={cn}",
            "your_data": your_data,
            "your_status": your_status,
            "risk_version": risk_version,
            "dp": dp_val,
            "gq": gq_val,
            "ad": ad_str,
            "qual": qual_val,
            "filter": filt,
            "is_carried": alt != "." and filt == "PASS",
        }

    if source_type == "sv":
        pr = sample_map.get("PR", ".")
        sr = sample_map.get("SR", ".")
        svtype = info_map.get("SVTYPE", "SV")
        svlen = info_map.get("SVLEN", ".")
        if filt != "PASS":
            your_status = "Low Quality Call (Unconfirmed)"
        elif gt == "1/1":
            your_status = "Detected (Homozygous)"
        elif gt in ("0/1", "1/0", "0/2", "1/2"):
            your_status = "Carrier (Heterozygous)"
        else:
            your_status = "Reference (Not Detected)"
        alt_short = alt if len(alt) <= 24 else f"{alt[:21]}..."
        your_data = f"{gt} ({svtype}) | PR: {pr} | SR: {sr} | GQ: {gq_str} | {filt}"
        return {
            "gt": gt,
            "alleles_called": f"{ref[:12]} / {alt_short}",
            "your_data": your_data,
            "your_status": your_status,
            "risk_version": alt_short if expected_risk_allele is None else expected_risk_allele,
            "dp": dp_val,
            "gq": gq_val,
            "ad": f"PR:{pr} SR:{sr}",
            "qual": qual_val,
            "filter": filt,
            "is_carried": gt in ("0/1", "1/1", "1/2") and filt == "PASS",
        }

    # SNP / Indel gVCF handling
    alleles = [ref] + ([a for a in alt.split(",") if a != "."] if alt != "." else [])
    sep = "/" if "/" in gt else ("|" if "|" in gt else "/")
    gt_parts = gt.replace("|", "/").split("/")
    called_bases: List[str] = []
    for idx_str in gt_parts:
        if idx_str.isdigit():
            idx = int(idx_str)
            if idx < len(alleles):
                called_bases.append(alleles[idx])
            else:
                called_bases.append(ref)
        else:
            called_bases.append("?")
    alleles_str = f"{called_bases[0]} {sep} {called_bases[1]}" if len(called_bases) == 2 else "/".join(called_bases)

    is_ref_block = alt == "." and "END" in info_map
    if gt == "./." or gt == ".":
        your_status = "No Call (Low Coverage)"
        is_carried = False
    elif alt == "." or gt == "0/0" or gt == "0|0" or gt == "0":
        your_status = "Reference (Not Detected)"
        is_carried = False
    elif filt != "PASS":
        your_status = "Low Quality Call (Unconfirmed)"
        is_carried = False
    elif gt in ("1/1", "1|1", "2/2", "1"):
        your_status = "Detected (Homozygous)"
        is_carried = True
    elif gt in ("0/1", "1/0", "0|1", "1|0", "0/2", "1/2"):
        your_status = "Carrier (Heterozygous)"
        is_carried = True
    else:
        your_status = "Variant Call"
        is_carried = True

    metrics = [f"{alleles_str} ({gt})"]
    if dp_val is not None:
        metrics.append(f"DP: {dp_val}")
    if min_dp_str != "." and is_ref_block:
        metrics.append(f"MIN_DP: {min_dp_str}")
    if ad_str != "." and not is_ref_block:
        metrics.append(f"AD: {ad_str}")
    if gq_val is not None and gq_val > 0:
        metrics.append(f"GQ: {gq_val}")
    metrics.append(filt if not is_ref_block else f"{filt} [Ref Block]")

    risk_ver = expected_risk_allele if expected_risk_allele else (alt if alt != "." else ref)

    return {
        "gt": gt,
        "alleles_called": alleles_str,
        "your_data": " | ".join(metrics),
        "your_status": your_status,
        "risk_version": risk_ver,
        "dp": dp_val,
        "gq": gq_val,
        "ad": ad_str,
        "qual": qual_val,
        "filter": filt,
        "is_carried": is_carried,
    }


class VCFGenomeEngine:
    """Reader for one analysis' VCF files (gVCF binary search + SV / CNV scans)."""

    def __init__(self, genome_dir: Path, sparse_index_path: Optional[Path] = None, manifest: Optional[Dict[str, Any]] = None):
        self.genome_dir = Path(genome_dir)
        self.snp_vcf_path = self.genome_dir / SNP_VCF_NAME
        self.sv_vcf_path = self.genome_dir / SV_VCF_NAME
        self.cnv_vcf_path = self.genome_dir / CNV_VCF_NAME
        self.sparse_index_path = Path(sparse_index_path) if sparse_index_path else self.genome_dir / "vcf_sparse_index.json"
        self.manifest: Dict[str, Any] = dict(manifest or {})
        self.header_end: int = 0
        self.blocks: List[Tuple[int, str, int]] = []
        self.chrom_ranges: Dict[str, Tuple[int, int]] = {}
        self._loaded = False
        self._load_lock = threading.Lock()

    def ensure_loaded(self) -> None:
        if self._loaded:
            return
        with self._load_lock:
            if not self._loaded:
                self._load()

    def _load(self) -> None:
        if not self.snp_vcf_path.exists():
            self._loaded = True
            return

        file_size = os.path.getsize(self.snp_vcf_path)
        if self.sparse_index_path.exists():
            try:
                cached = json.loads(self.sparse_index_path.read_text(encoding="utf-8"))
                if cached.get("file_size") == file_size and cached.get("version") == SPARSE_INDEX_VERSION:
                    self.header_end = cached["header_end"]
                    self.blocks = [tuple(x) for x in cached["blocks"]]
                    self.chrom_ranges = {k: tuple(v) for k, v in cached["chrom_ranges"].items()}
                    self._loaded = True
                    return
            except Exception:
                pass

        # Build sparse index every 256 KB
        step = 256 * 1024
        blocks: List[Tuple[int, str, int]] = []
        with open(self.snp_vcf_path, "rb") as f:
            pos = 0
            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith(b"#"):
                    self.header_end = pos
                    break
                pos = f.tell()

            for offset in range(self.header_end, file_size, step):
                f.seek(offset)
                if offset > self.header_end:
                    f.readline()
                cur = f.tell()
                line = f.readline()
                if not line:
                    break
                parts = line.split(b"\t", 3)
                if len(parts) >= 2:
                    blocks.append((cur, normalize_chrom(parts[0].decode("ascii", errors="replace")), int(parts[1])))

            # Sampling every 256 KB misses contigs whose records fit between two samples
            # (chrMT, chrY, small or filtered gVCFs). Scan only the gaps where the contig
            # changes (and the tail) and add an exact block at every contig's first record.
            blocks = self._add_contig_boundaries(f, blocks, file_size)

        chrom_ranges: Dict[str, List[int]] = {}
        for i, (off, ch, _p) in enumerate(blocks):
            if ch not in chrom_ranges:
                prev_off = blocks[i - 1][0] if i > 0 else self.header_end
                chrom_ranges[ch] = [prev_off, file_size]
            if i + 1 < len(blocks) and blocks[i + 1][1] != ch:
                chrom_ranges[ch][1] = blocks[i + 1][0]

        self.blocks = blocks
        self.chrom_ranges = {k: (v[0], v[1]) for k, v in chrom_ranges.items()}
        try:
            self.sparse_index_path.parent.mkdir(parents=True, exist_ok=True)
            self.sparse_index_path.write_text(
                json.dumps(
                    {
                        "version": SPARSE_INDEX_VERSION,
                        "file_size": file_size,
                        "header_end": self.header_end,
                        "blocks": self.blocks,
                        "chrom_ranges": self.chrom_ranges,
                    }
                ),
                encoding="utf-8",
            )
        except Exception:
            pass
        self._loaded = True

    @staticmethod
    def _add_contig_boundaries(f, blocks: List[Tuple[int, str, int]], file_size: int) -> List[Tuple[int, str, int]]:
        if not blocks:
            return blocks
        extra: List[Tuple[int, str, int]] = []
        gaps = [(blocks[i][0], blocks[i + 1][0], blocks[i][1]) for i in range(len(blocks) - 1) if blocks[i][1] != blocks[i + 1][1]]
        gaps.append((blocks[-1][0], file_size, blocks[-1][1]))
        for start, end, chrom in gaps:
            f.seek(start)
            current = chrom
            while f.tell() < end:
                off = f.tell()
                line = f.readline()
                if not line:
                    break
                parts = line.split(b"\t", 3)
                if len(parts) < 2 or line.startswith(b"#"):
                    continue
                ch = normalize_chrom(parts[0].decode("ascii", errors="replace"))
                if ch != current:
                    try:
                        extra.append((off, ch, int(parts[1])))
                    except ValueError:
                        pass
                    current = ch
        if not extra:
            return blocks
        merged = {b[0]: b for b in blocks}
        for b in extra:
            merged[b[0]] = b
        return [merged[k] for k in sorted(merged)]

    def _seek_near_position(self, f, chrom: str, target_pos: int) -> Tuple[int, int]:
        self.ensure_loaded()
        if chrom not in self.chrom_ranges:
            return (0, 0)
        lo, hi = self.chrom_ranges[chrom]
        for i, (off, ch, p) in enumerate(self.blocks):
            if ch == chrom:
                if p <= target_pos:
                    lo = self.blocks[i - 1][0] if (i > 0 and self.blocks[i - 1][1] == chrom) else self.chrom_ranges[chrom][0]
                if p > target_pos:
                    hi = off
                    break

        left, right = lo, hi
        while right - left > 16384:
            mid = (left + right) // 2
            f.seek(mid)
            f.readline()
            line = f.readline()
            if not line:
                right = mid
                continue
            parts = line.split(b"\t", 3)
            if len(parts) < 2:
                right = mid
                continue
            ch = normalize_chrom(parts[0].decode("ascii", errors="replace"))
            p = int(parts[1])
            if ch != chrom:
                if mid < (lo + hi) // 2:
                    left = mid
                else:
                    right = mid
            elif p < target_pos - 4000:
                left = mid
            else:
                right = mid

        start_seek = max(self.header_end, left - 2048)
        f.seek(start_seek)
        if start_seek > self.header_end:
            f.readline()
        return (start_seek, hi + 16384)

    def lookup_snp_indel_locus(
        self, chrom: str, target_pos: int, f_handle=None, include_context: int = 0
    ) -> Tuple[Optional[List[str]], List[List[str]]]:
        """
        Look up a single (chrom, target_pos) coordinate in snp-indel.genome.vcf.
        Returns (matching_parts, context_reads_parts).
        If multiple records overlap target_pos, prioritizes actual variant call (ALT != '.') over hom-ref block.
        """
        self.ensure_loaded()
        chrom = normalize_chrom(chrom)
        if not self.snp_vcf_path.exists() or chrom not in self.chrom_ranges:
            return (None, [])

        own_file = f_handle is None
        f = open(self.snp_vcf_path, "rb") if own_file else f_handle
        try:
            _start_seek, max_offset = self._seek_near_position(f, chrom, target_pos)
            best_match: Optional[List[str]] = None
            context_window: List[List[str]] = []

            while f.tell() <= max_offset:
                line = f.readline()
                if not line:
                    break
                parts = line.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")
                if len(parts) < 10:
                    continue
                ch = normalize_chrom(parts[0])
                if ch != chrom:
                    if best_match is not None:
                        break
                    continue
                p = int(parts[1])
                info = parts[7]
                end_p = p
                if info.startswith("END="):
                    try:
                        end_p = int(info.split(";")[0].split("=")[1])
                    except ValueError:
                        pass

                if include_context > 0 and abs(p - target_pos) <= include_context:
                    context_window.append(parts)

                if p == target_pos and parts[4] != ".":
                    best_match = parts
                    if include_context == 0:
                        break
                elif p <= target_pos <= end_p:
                    if best_match is None or best_match[4] == ".":
                        best_match = parts

                if p > target_pos + include_context:
                    break

            return (best_match, context_window[:200])
        finally:
            if own_file:
                f.close()

    def read_vcf_region(
        self,
        chrom: str,
        start_pos: int,
        end_pos: Optional[int] = None,
        variants_only: bool = False,
        limit: int = 50,
    ) -> List[Dict[str, Any]]:
        """
        Reads a series of VCF rows across snp-indel.genome.vcf, sv.vcf, and cnv.vcf
        within [start_pos, end_pos] on chrom.
        """
        self.ensure_loaded()
        chrom = normalize_chrom(chrom)
        if end_pos is None:
            end_pos = start_pos

        results: List[Dict[str, Any]] = []

        # 1. Read from snp-indel.genome.vcf via binary search
        if self.snp_vcf_path.exists() and chrom in self.chrom_ranges:
            with open(self.snp_vcf_path, "rb") as f:
                _s, max_offset = self._seek_near_position(f, chrom, start_pos)
                while f.tell() <= max_offset and len(results) < limit:
                    line = f.readline()
                    if not line:
                        break
                    raw_line = line.decode("utf-8", errors="replace").rstrip("\r\n")
                    parts = raw_line.split("\t")
                    if len(parts) < 10:
                        continue
                    ch = normalize_chrom(parts[0])
                    if ch != chrom:
                        if results:
                            break
                        continue
                    p = int(parts[1])
                    if p > end_pos:
                        break
                    info_map = parse_info_field(parts[7])
                    rec_end = int(info_map.get("END", p))
                    if rec_end < start_pos:
                        continue
                    if variants_only and parts[4] == ".":
                        continue
                    sample_map = parse_sample_field(parts[8], parts[9])
                    fmt = format_vcf_genotype_and_status(parts[3], parts[4], parts[5], parts[6], info_map, sample_map, "snp-indel")
                    vtype = determine_variant_type(parts[3], parts[4], info_map, "snp-indel")
                    results.append(
                        {
                            "source_file": SNP_VCF_NAME,
                            "chr": ch,
                            "position": p,
                            "end_position": rec_end,
                            "variant_id": parts[2] if parts[2] != "." else f"chr{ch}:{p}{parts[3]}>{parts[4]}",
                            "ref": parts[3],
                            "alt": parts[4],
                            "qual": fmt["qual"],
                            "filter": parts[6],
                            "variant_type": vtype,
                            "your_data": fmt["your_data"],
                            "risk_version": fmt["risk_version"],
                            "your_status": fmt["your_status"],
                            "dp": fmt["dp"],
                            "gq": fmt["gq"],
                            "ad": fmt["ad"],
                            "is_carried": fmt["is_carried"],
                            "raw_vcf_line": raw_line,
                        }
                    )

        # 2. Check sv.vcf for overlapping structural variants
        if self.sv_vcf_path.exists():
            with open(self.sv_vcf_path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    raw_line = line.rstrip("\r\n")
                    parts = raw_line.split("\t")
                    if len(parts) < 10 or normalize_chrom(parts[0]) != chrom:
                        continue
                    p = int(parts[1])
                    info_map = parse_info_field(parts[7])
                    rec_end = int(info_map.get("END", p))
                    if p <= end_pos and rec_end >= start_pos:
                        sample_map = parse_sample_field(parts[8], parts[9])
                        fmt = format_vcf_genotype_and_status(parts[3], parts[4], parts[5], parts[6], info_map, sample_map, "sv")
                        vtype = determine_variant_type(parts[3], parts[4], info_map, "sv")
                        results.append(
                            {
                                "source_file": SV_VCF_NAME,
                                "chr": chrom,
                                "position": p,
                                "end_position": rec_end,
                                "variant_id": parts[2],
                                "ref": parts[3],
                                "alt": parts[4],
                                "qual": fmt["qual"],
                                "filter": parts[6],
                                "variant_type": vtype,
                                "your_data": fmt["your_data"],
                                "risk_version": fmt["risk_version"],
                                "your_status": fmt["your_status"],
                                "dp": fmt["dp"],
                                "gq": fmt["gq"],
                                "ad": fmt["ad"],
                                "is_carried": fmt["is_carried"],
                                "raw_vcf_line": raw_line,
                            }
                        )

        # 3. Check cnv.vcf for overlapping CNVs
        if self.cnv_vcf_path.exists():
            with open(self.cnv_vcf_path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    raw_line = line.rstrip("\r\n")
                    parts = raw_line.split("\t")
                    if len(parts) < 10 or normalize_chrom(parts[0]) != chrom:
                        continue
                    p = int(parts[1])
                    info_map = parse_info_field(parts[7])
                    rec_end = int(info_map.get("END", p))
                    if variants_only and parts[4] == ".":
                        continue
                    if p <= end_pos and rec_end >= start_pos:
                        sample_map = parse_sample_field(parts[8], parts[9])
                        fmt = format_vcf_genotype_and_status(parts[3], parts[4], parts[5], parts[6], info_map, sample_map, "cnv")
                        vtype = determine_variant_type(parts[3], parts[4], info_map, "cnv")
                        results.append(
                            {
                                "source_file": CNV_VCF_NAME,
                                "chr": chrom,
                                "position": p,
                                "end_position": rec_end,
                                "variant_id": parts[2],
                                "ref": parts[3],
                                "alt": parts[4],
                                "qual": fmt["qual"],
                                "filter": parts[6],
                                "variant_type": vtype,
                                "your_data": fmt["your_data"],
                                "risk_version": fmt["risk_version"],
                                "your_status": fmt["your_status"],
                                "dp": fmt["dp"],
                                "gq": fmt["gq"],
                                "ad": fmt["ad"],
                                "is_carried": fmt["is_carried"],
                                "raw_vcf_line": raw_line,
                            }
                        )

        return results[:limit]

    def get_vcf_files_metadata(self) -> List[Dict[str, Any]]:
        uploaded = {f.get("canonical_name"): f for f in self.manifest.get("files", []) if isinstance(f, dict)}
        reference = str(self.manifest.get("reference") or workspace.DEFAULT_REFERENCE)
        files = []
        for name, label, desc in [
            (SNP_VCF_NAME, "SNP & Indel Genome VCF (gVCF)", "Single nucleotide variants, indels and homozygous-reference blocks"),
            (SV_VCF_NAME, "Structural Variants VCF (e.g. Manta)", "Deletions, insertions, tandem duplications and breakend translocations"),
            (CNV_VCF_NAME, "Copy Number Variants VCF (e.g. Canvas)", "Copy number gains, losses and regions of loss of heterozygosity"),
        ]:
            p = self.genome_dir / name
            exists = p.exists()
            size_mb = round(os.path.getsize(p) / (1024 * 1024), 2) if exists else 0
            files.append(
                {
                    "filename": name,
                    "original_filename": (uploaded.get(name) or {}).get("original_filename"),
                    "label": label,
                    "description": desc,
                    "exists": exists,
                    "size_mb": size_mb,
                    "reference": reference,
                }
            )
        return files


# ---------------------------------------------------------------------------
# One engine per workspace (its sparse index is loaded once and reused)
# ---------------------------------------------------------------------------

_engines: Dict[Path, VCFGenomeEngine] = {}
_engines_lock = threading.Lock()


def engine_for(ws: "workspace.Workspace") -> VCFGenomeEngine:
    key = Path(ws.root).resolve()
    with _engines_lock:
        eng = _engines.get(key)
        if eng is None:
            eng = VCFGenomeEngine(ws.inputs_dir, ws.sparse_index_path, ws.manifest())
            _engines[key] = eng
        return eng


def forget_engine(root: Path) -> None:
    """Drop the cached engine of a workspace (called when an analysis is deleted)."""
    with _engines_lock:
        _engines.pop(Path(root).resolve(), None)


def active_engine() -> VCFGenomeEngine:
    return engine_for(workspace.current())


class _ActiveEngineProxy:
    """`vcf_engine` module attribute: forwards to the engine of the active workspace.

    Attributes set on the proxy itself (e.g. by tests with monkeypatch) take precedence.
    """

    def __getattr__(self, name: str) -> Any:
        return getattr(active_engine(), name)


vcf_engine = _ActiveEngineProxy()
