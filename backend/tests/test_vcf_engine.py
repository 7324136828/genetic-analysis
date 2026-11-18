"""gVCF reader: sparse index coverage of every contig and chr-prefixed contig names."""

from genomescope.vcf_engine import VCFGenomeEngine
from genomescope.workspace import SNP_VCF_NAME


def _write(tmp_path, lines):
    d = tmp_path / "inputs"
    d.mkdir()
    (d / SNP_VCF_NAME).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return VCFGenomeEngine(d, tmp_path / "sparse.json")


HEAD = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]


def test_small_contigs_between_samples_are_indexed(tmp_path):
    # chr1 spans several 256 KB samples; chrY and chrM are tiny and fall between samples
    body = [f"chr1\t{1000 + i * 10}\t.\tA\t.\t.\tPASS\tEND={1009 + i * 10}\tGT:DP\t0/0:30" for i in range(9000)]
    body += ["chrY\t2781480\trs000001\tC\tT\t50\tPASS\t.\tGT:DP\t1:20"]
    body += ["chrM\t73\trs000002\tA\tG\t50\tPASS\t.\tGT:DP\t1/1:900"]
    eng = _write(tmp_path, HEAD + body)
    eng.ensure_loaded()
    assert {"1", "Y", "MT"} <= set(eng.chrom_ranges)
    row, _ = eng.lookup_snp_indel_locus("Y", 2781480)
    assert row and row[2] == "rs000001"
    row, _ = eng.lookup_snp_indel_locus("chrM", 73)
    assert row and row[2] == "rs000002"
    reads = eng.read_vcf_region("1", 50000, 50020)
    assert reads and all(r["chr"] == "1" for r in reads)


def test_sparse_index_is_cached_and_versioned(tmp_path):
    eng = _write(tmp_path, HEAD + ["chr2\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"])
    eng.ensure_loaded()
    assert (tmp_path / "sparse.json").exists()
    again = VCFGenomeEngine(tmp_path / "inputs", tmp_path / "sparse.json")
    again.ensure_loaded()
    assert again.chrom_ranges == eng.chrom_ranges and "2" in again.chrom_ranges
