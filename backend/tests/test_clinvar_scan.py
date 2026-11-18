"""
Offline tests for the ClinVar release scan (genomescope/clinvar_release.py, clinvar_scan.py,
genomescope/alleles.py). A small synthetic ClinVar summary and gVCF (fictitious genes,
variants and genotypes) make every behaviour deterministic; each test runs in its own
temporary analysis workspace.
"""

import gzip
import hashlib
import json
import shutil
import sqlite3
from pathlib import Path

import httpx
import pytest

from genomescope import clinvar_release as cr
from genomescope import clinvar_scan as cs
from genomescope import indexer, workspace
from genomescope import vcf_engine as vcf_mod
from genomescope.alleles import (
    DIFFERENT_ALT_STATUS,
    UNVERIFIED_ALLELE_STATUS,
    find_alt_index,
    normalize_vcf_allele,
)

HEADER = [
    "#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "HGNC_ID", "ClinicalSignificance", "ClinSigSimple",
    "LastEvaluated", "RS# (dbSNP)", "nsv/esv (dbVar)", "RCVaccession", "PhenotypeIDS", "PhenotypeList", "Origin",
    "OriginSimple", "Assembly", "ChromosomeAccession", "Chromosome", "Start", "Stop", "ReferenceAllele",
    "AlternateAllele", "Cytogenetic", "ReviewStatus", "NumberSubmitters", "Guidelines", "TestedInGTR", "OtherIDs",
    "SubmitterCategories", "VariationID", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
]
EXPERT = "reviewed by expert panel"
SINGLE = "criteria provided, single submitter"


def cv_row(vid, chrom, pos, ref, alt, sig, rcv, gene="GENE1", rs=-1, review=SINGLE, cond="Hereditary syndrome A",
           assembly="GRCh38", typ="single nucleotide variant", name=None, last="Sep 05, 2013"):
    d = {h: "-" for h in HEADER}
    d.update({
        "#AlleleID": str(vid + 10), "Type": typ, "Name": name or f"NM_000000.1({gene}):c.{vid}", "GeneID": "1",
        "GeneSymbol": gene, "ClinicalSignificance": sig, "LastEvaluated": last, "RS# (dbSNP)": str(rs),
        "RCVaccession": rcv, "PhenotypeList": cond, "Origin": "germline", "Assembly": assembly, "Chromosome": chrom,
        "ReviewStatus": review, "NumberSubmitters": "2", "VariationID": str(vid), "PositionVCF": str(pos),
        "ReferenceAlleleVCF": ref, "AlternateAlleleVCF": alt,
    })
    return "\t".join(d[h] for h in HEADER)


CLINVAR_ROWS = [
    cv_row(1001, "1", 500, "G", "A", "Pathogenic", "RCV000001001", gene="GENE1", rs=1001),
    cv_row(1002, "1", 600, "C", "T", "Likely pathogenic", "RCV000001002", gene="GENE2", rs=1002),     # sample has C>G
    cv_row(1003, "1", 700, "T", "A", "Pathogenic", "RCV000001003", gene="GENE3", rs=1003),            # multi-allelic 1/2
    cv_row(1004, "1", 800, "G", "T", "Likely pathogenic", "RCV000001004", gene="GENE4"),             # FAIL filter
    cv_row(1005, "1", 900, "A", "C", "Benign", "RCV000001005", gene="GENE5"),                        # benign
    cv_row(1006, "1", 950, "T", "G", "Uncertain significance", "RCV000001006", gene="GENE6"),        # VUS
    cv_row(1007, "1", 1200, "G", "C", "Pathogenic", "RCV000001007", gene="GENE7", rs=1007),          # sample hom-ref
    cv_row(2001, "1", 102, "GAG", "G", "Pathogenic", "RCV000002001", gene="GENEX", rs=900000001, typ="Deletion",
           review=EXPERT, name="NM_000000.1(GENEX):c.100_101del (p.Xaa34fs)"),                       # shifted indel
    cv_row(2002, "1", 2000, "C", "CAG", "Pathogenic", "RCV000002002", gene="GENEX", typ="Duplication"),  # not carried
    cv_row(3001, "1", -1, "na", "na", "Pathogenic", "RCV000003001", gene="BIGDEL", typ="copy number loss"),  # not representable
    cv_row(3002, "1", 500, "G", "A", "Pathogenic", "RCV000003002", gene="OLD", assembly="GRCh37"),   # other assembly
    cv_row(3003, "1", 3000, "A", "G", "not provided", "RCV000003003", gene="NOCLASS"),               # no assertion
    cv_row(4001, "1", 1100, "A", "G", "no classification for the single variant", "RCV000004001", gene="STAR", rs=4001),  # indexed for linking, never cataloged by the scan
]

FILLER_LINES = 5200  # keep the gVCF over 256 KB so the sparse byte-offset index has several blocks


def gvcf_text() -> str:
    head = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##contig=<ID=1>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_SAMPLE",
    ]
    f = "GT:DP:GQ:AD"
    body = [
        f"1\t100\t.\tC\t.\t50\tPASS\tEND=100\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t100\trs900000001\tCAG\tC\t60\tPASS\t.\t{f}\t0/1:30:99:15,15",
        f"1\t101\t.\tA\t.\t50\tPASS\tEND=499\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t500\trs1001\tG\tA\t50\tPASS\t.\t{f}\t0/1:30:99:15,15",
        f"1\t501\t.\tG\t.\t50\tPASS\tEND=599\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t600\t.\tC\tG\t50\tPASS\t.\t{f}\t0/1:30:99:15,15",
        f"1\t601\t.\tT\t.\t50\tPASS\tEND=699\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t700\t.\tT\tG,A\t50\tPASS\t.\t{f}\t1/2:30:99:0,15,15",
        f"1\t701\t.\tT\t.\t50\tPASS\tEND=799\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t800\t.\tG\tT\t50\tFAIL\t.\t{f}\t0/1:30:99:15,15",
        f"1\t801\t.\tG\t.\t50\tPASS\tEND=899\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t900\t.\tA\tC\t50\tPASS\t.\t{f}\t1/1:30:99:0,30",
        f"1\t901\t.\tA\t.\t50\tPASS\tEND=949\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t950\t.\tT\tG\t50\tPASS\t.\t{f}\t0/1:30:99:15,15",
        f"1\t951\t.\tT\t.\t50\tPASS\tEND=1300\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99",
        f"1\t1100\trs4001\tA\tG\t50\tPASS\t.\t{f}\t0/1:30:99:15,15",
    ]
    filler = [f"1\t{3000 + i * 10}\t.\tA\t.\t50\tPASS\tEND={3000 + i * 10 + 9}\tGT:DP:GQ:MIN_DP:AD\t0/0:40:0:30:99" for i in range(FILLER_LINES)]
    return "\n".join(head + body + filler) + "\n"


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    root = tmp_path_factory.mktemp("clinvar_world")
    summary = root / "variant_summary.txt.gz"
    with gzip.open(summary, "wt", encoding="utf-8") as f:
        f.write("\t".join(HEADER) + "\n" + "\n".join(CLINVAR_ROWS) + "\n")
    genome = root / "genome"
    genome.mkdir()
    gvcf = genome / vcf_mod.SNP_VCF_NAME
    gvcf.write_text(gvcf_text(), encoding="utf-8")
    index_path = root / "clinvar_grch38.db"
    stats = cr.build_index(summary, index_path, info={"md5": "abc", "last_modified": "Sun, 01 Nov 2026 14:17:51 GMT"})
    return {"root": root, "summary": summary, "genome": genome, "gvcf": gvcf, "index_path": index_path, "stats": stats}


@pytest.fixture
def index(world):
    ix = cr.ClinVarIndex(world["index_path"])
    yield ix
    ix.close()


@pytest.fixture
def ws(world, tmp_path):
    """A temporary analysis workspace holding the synthetic gVCF (active for the test)."""
    w = workspace.Workspace(tmp_path / "analysis").ensure_dirs()
    shutil.copy2(world["gvcf"], w.snp_vcf_path)  # copy2 keeps the mtime the scan records
    with workspace.activate(w):
        yield w
    vcf_mod.forget_engine(w.root)


@pytest.fixture
def catalog(ws):
    """Fresh catalog DB of the active workspace."""
    indexer.init_db()
    return ws.db_path


def legacy_row(**kw):
    base = dict(
        unique_key="SNP:1:0:rs0:RCV0", rcv_id="RCV0", all_rcv_ids="RCV0", chr="1", position=0, end_position=0, variant_id="rs0",
        gene="G", gene_name="", ref_allele="A", alt_allele="G", risk_version="G", risk_frequency=None, risk_frequency_text="Rare / Not Reported",
        variant_type="SNV", your_data="x", your_status="Carrier (Heterozygous)", genotype="0/1", alleles_called="A / G", condition="c",
        classification="Pathogenic", evidence="e", evidence_stars=1, review_status="r", hgvs="", molecular_consequence="", dp=30, gq=99,
        ad="15,15", qual=50.0, filter_status="PASS", is_carried=1, source_file=vcf_mod.SNP_VCF_NAME, raw_vcf_line="raw",
    )
    base.update(kw)
    return base


def rows(db, where="1=1", params=()):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    try:
        return [dict(r) for r in conn.execute(f"SELECT * FROM variants WHERE {where} ORDER BY id", params)]
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_classify_significance_tiers():
    assert cr.classify_significance("Pathogenic")[0] == cr.TIER_PATH
    assert cr.classify_significance("Pathogenic/Likely pathogenic")[0] == cr.TIER_PATH
    assert cr.classify_significance("Likely pathogenic")[0] == cr.TIER_LIKELY_PATH
    assert cr.classify_significance("Conflicting classifications of pathogenicity")[0] == cr.TIER_CONFLICTING
    assert cr.classify_significance("Uncertain significance")[0] == cr.TIER_VUS
    assert cr.classify_significance("VUS-mid") == (cr.TIER_VUS, "Uncertain Significance")
    assert cr.classify_significance("Benign/Likely benign")[0] == cr.TIER_BENIGN
    assert cr.classify_significance("drug response")[0] == cr.TIER_OTHER
    for none in ("-", "not provided", "", "no classification for the single variant"):
        assert cr.classify_significance(none)[0] == cr.TIER_NONE


def test_clinvar_date_and_conditions():
    assert cr.parse_clinvar_date("Sep 05, 2013") == "2013-09-05"
    assert cr.parse_clinvar_date("-") is None
    assert cr.clean_conditions("Syndrome A|not provided|Syndrome A|Neoplasm") == "Syndrome A | Neoplasm"


def test_normalized_key_is_anchor_independent():
    # two representations of one repeat deletion (gVCF anchor vs. annotation-source anchor)
    assert normalize_vcf_allele(1000, "CAG", "C") == normalize_vcf_allele(1001, "AGA", "A")
    # a different event does not collapse onto it
    assert normalize_vcf_allele(1000, "CAG", "C") != normalize_vcf_allele(1000, "C", "CAG")


def test_find_alt_index_exact_equivalent_and_different():
    assert find_alt_index(500, "G", ["A"], 500, "G", "A") == (1, "exact")
    assert find_alt_index(700, "T", ["G", "A"], 700, "T", "A") == (2, "exact")
    # same position, different substitution: must NOT match
    assert find_alt_index(600, "C", ["G"], 600, "C", "T")[0] is None
    # deletion of AG written at a different offset of the AGAG repeat
    assert find_alt_index(100, "CAG", ["C"], 102, "GAG", "G") == (1, "equivalent")
    # unusable ClinVar allele cannot be attributed
    assert find_alt_index(500, "G", ["A"], 500, "N", "N") == (None, "undecidable")


# ---------------------------------------------------------------------------
# Release download + index
# ---------------------------------------------------------------------------

def test_build_index_filters_and_counts(world, index):
    st = world["stats"]
    assert st["rows_grch38"] == len(CLINVAR_ROWS) - 1                       # GRCh37 row ignored
    assert st["rows_not_vcf_representable"] == 1 and st["not_vcf_representable_by_type"] == {"copy number loss": 1}
    assert st["rows_no_clinical_assertion"] == 2
    assert st["rows_indexed"] == 11                 # incl. the unclassified alleles (kept for linking)
    assert st["release_date"] == "2026-11-01" and st["release_date_source"] == "NCBI Last-Modified header"
    assert index.allele(1001)["alt"] == "A"
    # benign / VUS keep their allele keys; only tier >= conflicting has full annotation
    assert index.allele(1005)["tier"] == cr.TIER_BENIGN
    assert set(index.annotations_for([1001, 1005, 1006])) == {1001}
    assert [a["variation_id"] for a in index.alleles_by_rsid(["rs4001"])["rs4001"]] == [4001]
    assert index.variations_for_rcvs(["RCV000001001", "RCV999999999"]) == {"RCV000001001": 1001}
    assert cr.annotations_from_summary(world["summary"], [1006])[1006]["classification"] == "Uncertain Significance"


def test_build_index_rejects_unexpected_layout(tmp_path):
    bad = tmp_path / "variant_summary.txt.gz"
    with gzip.open(bad, "wt") as f:
        f.write("#AlleleID\tType\tAssembly\n1\tx\tGRCh38\n")
    with pytest.raises(cr.ClinVarError, match="missing expected column"):
        cr.build_index(bad, tmp_path / "x.db", info={})
    assert not (tmp_path / "x.db").exists()


def _mock_client(body: bytes, published_md5: str) -> httpx.Client:
    def handler(request: httpx.Request) -> httpx.Response:
        if str(request.url).endswith(".md5"):
            return httpx.Response(200, text=f"{published_md5}  /ftp/{cr.SUMMARY_NAME}\n")
        return httpx.Response(200, content=body, headers={"Last-Modified": "Sun, 01 Nov 2026 14:17:51 GMT", "Content-Length": str(len(body))})
    return httpx.Client(transport=httpx.MockTransport(handler), follow_redirects=True)


def test_download_verifies_md5_and_writes_sidecar(tmp_path):
    body = b"pretend-gzip-bytes" * 1000
    info = cr.download_release(tmp_path, client=_mock_client(body, hashlib.md5(body).hexdigest()))
    assert info["md5"] == hashlib.md5(body).hexdigest() and (tmp_path / cr.SUMMARY_NAME).read_bytes() == body
    assert cr.read_release_info(tmp_path)["last_modified"].startswith("Sun, 01 Nov 2026")
    # an already verified file is kept (no second download)
    again = cr.download_release(tmp_path, client=_mock_client(b"different", hashlib.md5(body).hexdigest()))
    assert again["md5"] == info["md5"] and (tmp_path / cr.SUMMARY_NAME).read_bytes() == body


def test_download_rejects_md5_mismatch(tmp_path):
    with pytest.raises(cr.ClinVarError, match="MD5 mismatch"):
        cr.download_release(tmp_path, client=_mock_client(b"corrupted download", "0" * 32))
    assert not (tmp_path / cr.SUMMARY_NAME).exists()
    assert not list(tmp_path.glob("*.part"))


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

def test_scan_matches_alleles_not_positions(world, index):
    res = cs.scan_gvcf(world["gvcf"], index)
    by_vid = {f.variation_id: f for f in res.findings}

    assert 1002 not in by_vid                       # sample has C>G, ClinVar says C>T: not carried
    assert 1007 not in by_vid and 2002 not in by_vid  # hom-ref / absent
    assert 1005 not in by_vid                       # benign is counted, never a finding
    assert {1001, 1003, 1004, 1006, 2001} == set(by_vid)

    assert by_vid[1001].match == "exact" and by_vid[1001].copies == 1 and by_vid[1001].ploidy == 2
    assert by_vid[1003].alt == "A" and by_vid[1003].alt_index == 2          # second ALT of a 1/2 call
    assert not by_vid[1004].passes                                           # FAIL filter kept, flagged
    # same deletion written at another offset of the repeat
    assert by_vid[2001].match == "equivalent" and by_vid[2001].pos == 100 and by_vid[2001].cv_pos == 102
    assert res.stats["equivalent_representation_matches"] == 1
    assert res.stats["clinvar_matches_by_tier"]["benign / likely benign"]["PASS"] == 1
    assert res.stats["unresolved_high_impact_indels"] == 0


def test_scan_coverage_distinguishes_not_carried_from_not_callable(world, index):
    cov = cs.scan_gvcf(world["gvcf"], index).coverage
    assert cov["high_impact_loci"] == 7             # 1001,1002,1003,1004,1007,2001 + 2002
    assert cov["confirmed_reference"] == 1          # 1007: PASS block, depth 30
    assert cov["variant_called"] == 5
    assert cov["no_gvcf_record"] == 1               # 2002 lies in no gVCF row
    assert cov["confirmed_reference"] + cov["reference_low_confidence"] + cov["variant_called"] + cov["no_gvcf_record"] == 7


# ---------------------------------------------------------------------------
# Catalog reconciliation
# ---------------------------------------------------------------------------

def _seed_legacy(db):
    indexer.upsert_records([
        # the sample has C>G at 600 but this row claims the ClinVar C>T allele is carried
        legacy_row(unique_key="SNP:1:600:rs1002:RCV000001002", rcv_id="RCV000001002", all_rcv_ids="RCV000001002", position=600,
                   variant_id="rs1002", ref_allele="C", alt_allele="G", risk_version="T", classification="Likely Pathogenic"),
        # a correct carried row that must be updated in place, not duplicated
        legacy_row(unique_key="SNP:1:500:rs1001:RCV000001001", rcv_id="RCV000001001", all_rcv_ids="RCV000001001", position=500,
                   variant_id="rs1001", ref_allele="G", alt_allele="A", risk_version="A"),
        # a non-carried marker with the old ref == alt artefact
        legacy_row(unique_key="SNP:1:1200:rs1007:RCV000001007", rcv_id="RCV000001007", all_rcv_ids="RCV000001007", position=1200,
                   variant_id="rs1007", ref_allele="G", alt_allele="G", risk_version="G", your_status="Reference (Not Detected)", is_carried=0),
        # claims carried, but no ClinVar allele can be tied to it
        legacy_row(unique_key="SNP:1:777:rs777:RCV000009999", rcv_id="RCV000009999", all_rcv_ids="RCV000009999", position=777, variant_id="rs777"),
        # retired RCV, but the allele still exists in ClinVar (now without a single-variant classification)
        legacy_row(unique_key="SNP:1:1100:rs4001:RCV000999999", rcv_id="RCV000999999", all_rcv_ids="RCV000999999", position=1100,
                   variant_id="rs4001", ref_allele="A", alt_allele="A", risk_version="G", classification="Drug Response", evidence="Moderate"),
        # benign allele the sample carries, wrongly stored as not detected
        legacy_row(unique_key="SNP:1:900:rs1005:RCV000001005", rcv_id="RCV000001005", all_rcv_ids="RCV000001005", position=900,
                   variant_id="rs1005", ref_allele="A", alt_allele="A", risk_version="C", your_status="Reference (Not Detected)", is_carried=0, classification="Benign"),
        # no ClinVar record at all, but the row's own substitution is verifiably carried
        legacy_row(unique_key="SNP:1:500:rs9:RCV000888888", rcv_id="RCV000888888", all_rcv_ids="RCV000888888", position=500,
                   variant_id="rs9", ref_allele="G", alt_allele="A", risk_version="A", classification="Risk Factor"),
        # direct VCF read: makes no ClinVar claim and must be left alone
        legacy_row(unique_key="RAW:snp:1:600:x", rcv_id="VCF-DIRECT-READ", all_rcv_ids="VCF-DIRECT-READ", position=600, variant_id="chr1:600CG", classification="Variant Call"),
    ])


def test_reconcile_adds_missing_carriers_and_corrects_false_ones(world, index, catalog, tmp_path):
    _seed_legacy(catalog)
    summary = cs.run_scan(world["gvcf"], index, progress=None, summary_gz=world["summary"],
                          report_path=tmp_path / "report.json", tsv_path=tmp_path / "carried.tsv")
    assert summary["postcondition_violations"] == []

    # the shifted GENEX indel is now in the catalog even though no legacy row existed for it
    genex = rows(catalog, "clinvar_variation_id = '2001'")
    assert len(genex) == 1 and genex[0]["unique_key"].startswith("CLINVAR:") and genex[0]["variant_id"] == "rs900000001"
    assert genex[0]["is_carried"] == 1 and genex[0]["your_status"] == "Carrier (Heterozygous)"
    assert genex[0]["clinvar_match"] == "equivalent" and genex[0]["clinvar_release"] == "2026-11-01"
    assert genex[0]["classification"] == "Pathogenic" and "Expert Panel" in genex[0]["evidence"] and genex[0]["gene"] == "GENEX"
    assert genex[0]["condition"] == "Hereditary syndrome A" and genex[0]["rcv_id"] == "RCV000002001"
    assert (genex[0]["ref_allele"], genex[0]["alt_allele"], genex[0]["risk_version"]) == ("CAG", "C", "G")

    # multi-allelic 1/2 call: the matched ALT is recorded, not the whole list
    gt3 = rows(catalog, "clinvar_variation_id = '1003'")[0]
    assert gt3["alt_allele"] == "A" and gt3["is_carried"] == 1

    # FAIL call is recorded for review but never counted as carried
    f4 = rows(catalog, "clinvar_variation_id = '1004'")[0]
    assert f4["is_carried"] == 0 and f4["your_status"] == cs.LOW_QUALITY_STATUS

    # benign and (by default) VUS are not cataloged
    assert rows(catalog, "clinvar_variation_id = '1006'") == []
    assert all(not r["unique_key"].startswith("CLINVAR:") for r in rows(catalog, "clinvar_variation_id = '1005'"))  # only the pre-existing legacy row

    # legacy false carrier corrected, not deleted
    fc = rows(catalog, "unique_key = 'SNP:1:600:rs1002:RCV000001002'")[0]
    assert fc["is_carried"] == 0 and fc["your_status"] == DIFFERENT_ALT_STATUS and fc["clinvar_match"] == "adjudicated_not_carried"
    assert any(x["variant_id"] == "rs1002" for x in json.loads((tmp_path / "report.json").read_text())["catalog_changes_detail"]["false_carriers_corrected"])

    # correct legacy row updated in place; no duplicate scan row for the same variation
    ok = rows(catalog, "clinvar_variation_id = '1001'")
    assert len(ok) == 1 and ok[0]["unique_key"] == "SNP:1:500:rs1001:RCV000001001" and ok[0]["is_carried"] == 1

    # non-carried marker: ref == alt artefact repaired, status stays reference
    ref_row = rows(catalog, "unique_key = 'SNP:1:1200:rs1007:RCV000001007'")[0]
    assert (ref_row["ref_allele"], ref_row["alt_allele"]) == ("G", "C") and ref_row["is_carried"] == 0
    assert "Reference" in ref_row["your_status"]

    # unverifiable carried row is demoted, with a trace in the report
    unv = rows(catalog, "unique_key = 'SNP:1:777:rs777:RCV000009999'")[0]
    assert unv["is_carried"] == 0 and unv["your_status"] == UNVERIFIED_ALLELE_STATUS and unv["clinvar_match"] == "unverified"

    # retired RCV -> re-linked by rsID; genotype verified; stale classification flagged, not silently trusted
    st4 = rows(catalog, "unique_key = 'SNP:1:1100:rs4001:RCV000999999'")[0]
    assert st4["clinvar_variation_id"] == "4001" and st4["is_carried"] == 1 and st4["clinvar_match"] == "adjudicated_carried"
    assert (st4["ref_allele"], st4["alt_allele"]) == ("A", "G") and "no germline classification" in st4["evidence"]
    detail = json.loads((tmp_path / "report.json").read_text())["catalog_changes_detail"]
    assert any(x["variant_id"] == "rs4001" for x in detail["stale_classification_rows"])

    # benign false negative fixed; it is not reported as a scan miss
    bn = rows(catalog, "unique_key = 'SNP:1:900:rs1005:RCV000001005'")[0]
    assert bn["is_carried"] == 1 and bn["your_status"] == "Detected (Homozygous)"
    assert any(x["variant_id"] == "rs1005" for x in detail["false_negatives_corrected"])
    assert detail["scan_vs_adjudication_disagreements"] == []

    # no ClinVar record, but the row's own allele is verified in the gVCF: kept, and labelled as such
    own = rows(catalog, "unique_key = 'SNP:1:500:rs9:RCV000888888'")[0]
    assert own["is_carried"] == 1 and own["clinvar_match"] == "legacy_allele_verified"

    # direct reads untouched
    raw = rows(catalog, "unique_key LIKE 'RAW:%'")[0]
    assert raw["is_carried"] == 1 and raw["clinvar_match"] is None

    assert (tmp_path / "carried.tsv").read_text().count("\n") == len(cs.scan_gvcf(world["gvcf"], index).findings) + 1


def test_reconcile_is_idempotent_and_backs_up(world, index, catalog, ws, tmp_path):
    _seed_legacy(catalog)
    kw = dict(progress=None, summary_gz=world["summary"], report_path=tmp_path / "r.json", tsv_path=tmp_path / "c.tsv")
    first = cs.run_scan(world["gvcf"], index, **kw)
    snapshot = [(r["unique_key"], r["is_carried"], r["your_status"], r["clinvar_match"]) for r in rows(catalog)]
    second = cs.run_scan(world["gvcf"], index, **kw)
    assert [(r["unique_key"], r["is_carried"], r["your_status"], r["clinvar_match"]) for r in rows(catalog)] == snapshot
    assert first["catalog"]["rows_added"] > 0 and second["catalog"]["rows_added"] == 0
    assert second["postcondition_violations"] == []
    assert len(list(ws.backup_dir.glob("genetic_variants.db.pre-clinvar-scan-*"))) >= 1


def test_include_vus_reads_annotation_from_release_file(world, index, catalog, tmp_path):
    kw = dict(progress=None, report_path=tmp_path / "r.json", tsv_path=tmp_path / "c.tsv")
    cs.run_scan(world["gvcf"], index, min_tier=cr.TIER_VUS, summary_gz=world["summary"], **kw)
    vus = rows(catalog, "clinvar_variation_id = '1006'")
    assert len(vus) == 1 and vus[0]["classification"] == "Uncertain Significance" and vus[0]["is_carried"] == 1


def test_dry_run_changes_nothing(world, index, catalog, tmp_path):
    _seed_legacy(catalog)
    before = rows(catalog)
    summary = cs.run_scan(world["gvcf"], index, dry_run=True, progress=None, report_path=tmp_path / "r.json", tsv_path=tmp_path / "c.tsv")
    assert rows(catalog) == before and summary["catalog"] is None
    assert (tmp_path / "r.json").exists()


def test_verify_catalog_flags_unverified_carriers(world, index, catalog):
    indexer.upsert_records([legacy_row(unique_key="SNP:1:5:rs5:RCV5", position=5)])
    problems = cs.verify_catalog(cs.scan_gvcf(world["gvcf"], index))
    assert any("marked carried without a verified ClinVar allele match" in p for p in problems)
    assert any("VCV000001001" in p for p in problems)   # carried per scan, absent from the catalog


# ---------------------------------------------------------------------------
# Legacy (position-keyed) builder must no longer invent carriers
# ---------------------------------------------------------------------------

def test_legacy_builder_requires_the_clinvar_allele():
    cv = {"chr": "1", "position": 600, "end_position": 600, "variant_id": "rs1002", "rcv_id": "RCV1", "all_rcv_ids": "RCV1",
          "gene": "G", "ref_allele": "C", "alt_allele": "T", "af": None, "classification": "Pathogenic", "review_status": SINGLE,
          "condition": "c", "hgvs": "", "clinvar_type": "single nucleotide variant"}
    line = "1\t600\t.\tC\tG\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15".split("\t")
    row = indexer.build_record_from_clinvar_and_vcf(cv, line)
    assert row["is_carried"] == 0 and row["your_status"] == DIFFERENT_ALT_STATUS

    line_ok = "1\t600\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15".split("\t")
    assert indexer.build_record_from_clinvar_and_vcf(cv, line_ok)["is_carried"] == 1

    # ClinVar alt on the second ALT but the call only has the first: not carried
    line_multi = "1\t600\t.\tC\tG,T\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15,0".split("\t")
    assert indexer.build_record_from_clinvar_and_vcf(cv, line_multi)["is_carried"] == 0

    # degenerate / unusable ClinVar allele (ref == alt) can never be "carried"
    bad = dict(cv, alt_allele="C")
    assert indexer.build_record_from_clinvar_and_vcf(bad, line_ok)["your_status"] == UNVERIFIED_ALLELE_STATUS


# ---------------------------------------------------------------------------
# Status (GET /api/clinvar/status)
# ---------------------------------------------------------------------------

def test_status_reports_each_state(world, catalog, ws, tmp_path):
    assert cs.status(tmp_path / "missing.db")["state"] == "index_missing"

    idx = world["index_path"]
    assert cs.status(idx)["state"] == "never_run"

    index = cr.ClinVarIndex(idx)
    try:
        cs.run_scan(world["gvcf"], index, dry_run=True, progress=None, report_path=tmp_path / "r.json", tsv_path=tmp_path / "c.tsv")
        assert cs.status(idx)["state"] == "never_run"            # a dry run does not verify the catalog
        cs.run_scan(world["gvcf"], index, progress=None, summary_gz=world["summary"], report_path=tmp_path / "r.json", tsv_path=tmp_path / "c.tsv")
    finally:
        index.close()
    st = cs.status(idx)
    assert st["state"] == "verified" and st["index"]["release_date"] == "2026-11-01"
    assert st["catalog"]["carried_rows_without_verified_allele"] == 0
    assert st["last_scan"]["coverage"]["high_impact_loci"] == 7

    # the gVCF changed after the scan -> the verification no longer applies
    extra = "\t".join(["1", "999999", ".", "A", ".", "50", "PASS", "END=999999", "GT", "0/0"]) + "\n"
    ws.snp_vcf_path.write_text(gvcf_text() + extra, encoding="utf-8")
    assert cs.status(idx)["state"] == "stale"


def test_status_without_selected_analysis(world):
    # no active workspace: the index is described, the per-analysis state is not
    assert workspace.get() is None
    st = cs.status(world["index_path"])
    assert st["state"] == "no_analysis" and st["index"]["release_date"] == "2026-11-01" and st["catalog"] is None


def test_status_endpoint_is_served(monkeypatch, bare_api):
    monkeypatch.setattr(cs, "status", lambda *a, **k: {"state": "verified", "message": "ok"})
    r = bare_api.get("/api/clinvar/status")
    assert r.status_code == 200 and r.json()["state"] == "verified"
