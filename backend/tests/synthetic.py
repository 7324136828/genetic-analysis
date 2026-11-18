"""
Synthetic test genome (no real person's data).

Every genotype below is invented to exercise a specific code path (hom-ref blocks, no-calls,
split rows at one site, FAIL calls, haploid chrX calls, left-aligned indels, SV / CNV rows).
Positions are public GRCh38 coordinates of well-known ClinVar / dbSNP records so that the
offline NCBI fixtures in fixtures/ncbi/ line up with them.
"""

from __future__ import annotations

from typing import Dict, List

from genomescope import indexer
from genomescope.vcf_engine import engine_for
from genomescope.workspace import Workspace

SAMPLE = "SYNTHETIC_SAMPLE"
FMT = "GT:DP:GQ:AD"
BLOCK_FMT = "GT:DP:GQ:MIN_DP"


def _block(chrom: str, start: int, end: int, ref: str = "A", dp: int = 35) -> str:
    return f"{chrom}\t{start}\t.\t{ref}\t.\t.\tPASS\tEND={end}\t{BLOCK_FMT}\t0/0:{dp}:0:{dp - 5}"


def _call(chrom: str, pos: int, vid: str, ref: str, alt: str, gt: str, filt: str = "PASS", dp: int = 32, gq: int = 99, ad: str = "16,16", qual: str = "60.0") -> str:
    return f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t.\t{FMT}\t{gt}:{dp}:{gq}:{ad}"


def _filler(chrom: str, start: int, count: int, step: int = 1000) -> List[str]:
    """Hom-ref blocks to make the file large enough for several sparse-index samples."""
    return [_block(chrom, start + i * step, start + i * step + step - 1) for i in range(count)]


def gvcf_lines() -> List[str]:
    head = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##reference=file:///reference/GRCh38.fa",
        "##contig=<ID=1,length=248956422>",
        "##contig=<ID=2,length=242193529>",
        "##contig=<ID=6,length=170805979>",
        "##contig=<ID=7,length=159345973>",
        "##contig=<ID=12,length=133275309>",
        "##contig=<ID=17,length=83257441>",
        "##contig=<ID=X,length=156040895>",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
    ]
    body: List[str] = []
    # chr1 ------------------------------------------------------------------
    body.append("1\t10050\t.\tA\t.\t0\tPASS\tEND=10060\tGT:DP:GQ:MIN_DP\t./.:0:0:0")         # no-call
    body += _filler("1", 20000, 4000)                                                       # ~4 Mb of hom-ref blocks
    body.append(_block("1", 11796000, 11796320, "C"))
    body.append(_call("1", 11796321, "rs1801133", "G", "A", "0/1", ad="15,17"))             # het carrier
    body.append(_block("1", 11796322, 11797000, "C"))
    body.append(_block("1", 100186000, 100186847, "T"))
    body.append(_call("1", 100186848, ".", "C", "CAT", "1/1", filt="FAIL", dp=12, gq=10, ad="2,10"))   # FAIL row first
    body.append(_call("1", 100186848, "rs4364901", "C", "T", "0/1"))                      # PASS SNV at same POS
    body.append(_block("1", 100186849, 100187000, "G"))
    body.append(_block("1", 236414000, 236414461, "A"))
    body.append(_call("1", 236414462, "rs61021118", "C", "CT", "0/1"))                    # insertion row ...
    body.append(_call("1", 236414462, "rs625443", "C", "T", "0/1"))                       # ... then SNV row
    body.append(_block("1", 236414463, 236415000, "G"))
    body += _filler("1", 240000000, 2000)
    # chr2 ------------------------------------------------------------------
    body.append(_block("2", 47414400, 47414419, "G"))
    body.append(_call("2", 47414420, "rs900000002", "T", "C", "0/1"))                     # call at the deleted base
    body.append(_block("2", 47414421, 47415000, "A"))
    # chr6 ------------------------------------------------------------------
    body.append(_block("6", 26092800, 26092912, "A"))
    body.append(_block("6", 26092913, 26093000, "G"))                                      # hom-ref at HFE locus
    # chr7 ------------------------------------------------------------------
    body.append(_block("7", 117559500, 117559700, "T"))                                    # hom-ref at CFTR locus
    body.append(_block("7", 157480000, 157480040, "A"))
    body.append(_call("7", 157480041, "rs59826218", "CTG", "C", "0/1", ad="14,15"))        # left-aligned deletion
    body.append(_block("7", 157480044, 157481000, "A"))
    # chr12 -----------------------------------------------------------------
    body.append(_block("12", 111803000, 111803961, "C"))
    body.append(_call("12", 111803962, "rs671", "G", "A", "1/1", dp=30, gq=90, ad="0,30"))  # homozygous ALT
    body.append(_block("12", 111803963, 111805000, "C"))
    # chr17 -----------------------------------------------------------------
    body.append(_block("17", 46337000, 46337317, "A"))
    body.append(_call("17", 46337318, ".", "T", "C", "1/1", filt="FAIL", dp=8, gq=5, ad="0,8"))
    body.append(_block("17", 46337319, 46338000, "A"))
    # chrX (haploid call) ----------------------------------------------------
    body.append(_block("X", 150596000, 150596741, "G"))
    body.append(f"X\t150596742\trs900000003\tT\tC\t60.0\tPASS\t.\t{FMT}\t1:30:99:0,30")
    body.append(_block("X", 150596743, 150597000, "G"))
    return head + body


def sv_lines() -> List[str]:
    return [
        "##fileformat=VCFv4.1",
        "##reference=file:///reference/GRCh38.fa",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
        "2\t5000000\tMantaDEL:1:0:1:0:0:0\tA\t<DEL>\t800\tPASS\tEND=5004000;SVTYPE=DEL;SVLEN=-4000\tGT:GQ:PR:SR\t0/1:500:20,10:15,8",
        "6\t3000000\tMantaINS:2:0:1:0:0:0\tC\t<INS>\t650\tPASS\tEND=3000001;SVTYPE=INS;SVLEN=300\tGT:GQ:PR:SR\t0/1:400:18,9:12,7",
        "7\t9000000\tMantaDEL:3:0:1:0:0:0\tG\t<DEL>\t100\tMinGQ\tEND=9001000;SVTYPE=DEL;SVLEN=-1000\tGT:GQ:PR:SR\t0/1:20:5,2:3,1",
    ]


def cnv_lines() -> List[str]:
    return [
        "##fileformat=VCFv4.1",
        "##reference=file:///reference/GRCh38.fa",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
        "1\t1000001\tCanvas:REF:1:1000001-2000000\tN\t.\t.\tPASS\tSVTYPE=CNV;END=2000000\tRC:BC:CN\t60:500:2",
        "1\t30000001\tCanvas:LOSS:1:30000001-30045000\tN\t<CN0>\t40\tPASS\tSVTYPE=CNV;END=30045000;CNVLEN=45000\tRC:BC:CN:MCC\t30:40:1:.",
        "12\t50000001\tCanvas:GAIN:12:50000001-50060000\tN\t<CN3>\t35\tPASS\tSVTYPE=CNV;END=50060000;CNVLEN=60000\tRC:BC:CN:MCC\t90:55:3:.",
    ]


# MyVariant-style ClinVar records joined to the gVCF by the legacy catalog builder.
SEED_CLINVAR: List[Dict] = [
    {
        "rcv_id": "RCV000020060", "all_rcv_ids": "RCV000020060", "chr": "12", "position": 111803962, "end_position": 111803962,
        "variant_id": "rs671", "gene": "ALDH2", "ref_allele": "G", "alt_allele": "A", "af": 0.0061,
        "classification": "Drug Response", "review_status": "reviewed by expert panel",
        "condition": "Alcohol flushing response (synthetic test record)", "hgvs": "c.1510G>A (p.Glu504Lys)", "clinvar_type": "single nucleotide variant",
    },
    {
        "rcv_id": "RCV900000001", "all_rcv_ids": "RCV900000001", "chr": "1", "position": 11796321, "end_position": 11796321,
        "variant_id": "rs1801133", "gene": "MTHFR", "ref_allele": "G", "alt_allele": "A", "af": 0.31,
        "classification": "Risk Factor", "review_status": "practice guideline",
        "condition": "Folate metabolism trait (synthetic test record)", "hgvs": "c.665C>T (p.Ala222Val)", "clinvar_type": "single nucleotide variant",
    },
]


def write_inputs(ws: Workspace) -> None:
    ws.ensure_dirs()
    ws.snp_vcf_path.write_text("\n".join(gvcf_lines()) + "\n", encoding="utf-8")
    ws.sv_vcf_path.write_text("\n".join(sv_lines()) + "\n", encoding="utf-8")
    ws.cnv_vcf_path.write_text("\n".join(cnv_lines()) + "\n", encoding="utf-8")


def build_catalog(ws: Workspace, monkeypatch=None) -> None:
    """Seed the catalog offline (requires `ws` to be the active workspace)."""
    indexer.init_db()
    eng = engine_for(ws)
    eng.ensure_loaded()
    rows = []
    for cv in SEED_CLINVAR:
        parts, _ctx = eng.lookup_snp_indel_locus(cv["chr"], cv["position"])
        row = indexer.build_record_from_clinvar_and_vcf(cv, parts)
        assert row is not None, cv["variant_id"]
        rows.append(row)
    indexer.upsert_records(rows)
    original = indexer.annotate_genomic_interval_genes
    indexer.annotate_genomic_interval_genes = lambda chrom, start, end: ("SYNGENE1", "synthetic gene")  # no network
    try:
        indexer.index_cnv_and_sv_files()
    finally:
        indexer.annotate_genomic_interval_genes = original

