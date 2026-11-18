"""
Tests for the NCBI integration (genomescope/ncbi_client.py, genomescope/ncbi_service.py and
the /api/ncbi/* endpoints).

- Parser tests use public NCBI responses saved under tests/fixtures/ncbi/ (trimmed).
- Endpoint tests monkeypatch NcbiClient._transport with a fixture router (no network) and
  run against a fresh copy of the synthetic analysis (tests/synthetic.py) for every test.
- Live tests run only with RUN_LIVE_NCBI=1.
"""

import json
import os
import threading
import time
import urllib.parse
from pathlib import Path

import httpx
import pytest

from genomescope import indexer, ncbi_client, ncbi_service
from genomescope import vcf_engine as vcf_engine_mod
from genomescope.ncbi_client import (
    NcbiClient,
    NcbiError,
    RateLimiter,
    alleles_equivalent,
    allele_names,
    dbsnp_spdis,
    extract_clinvar_spdi,
    normalize_ncbi_date,
    parse_clinvar_summary,
    parse_dbsnp_summary,
    parse_esearch,
    parse_esummary_result,
    parse_gene_summary,
    parse_global_maf,
    parse_pubmed_summary,
    parse_refsnp_frequencies,
    review_status_to_stars,
    risk_allele_frequency_from_dbsnp,
    risk_allele_frequency_from_refsnp,
    sanitize_phrase,
    spdi_to_vcf,
    translate_search_term,
)

FIX = Path(__file__).parent / "fixtures" / "ncbi"
VCF_AVAILABLE = True  # the synthetic analysis always provides a gVCF
pytestmark = pytest.mark.usefixtures("job_ws")


def load(name):
    return json.loads((FIX / name).read_text(encoding="utf-8"))


CLINVAR = load("clinvar_esummary.json")
SNP = load("snp_esummary.json")
GENE = load("gene_esummary.json")
PUBMED = load("pubmed_esummary.json")
ELINK = load("clinvar_elink_pubmed.json")
REFSNP_334 = load("refsnp_334.json")
FASTA_START = 117559585
FASTA_SEQ = "".join(l for l in (FIX / "efetch_nc000007_117559585_117559595.fasta").read_text().splitlines() if not l.startswith(">"))


def cv(uid):
    return CLINVAR["result"][uid]


def snp(uid):
    return SNP["result"][uid]


# ---------------------------------------------------------------------------
# Fixture-backed fake transport
# ---------------------------------------------------------------------------

ESEARCH_MAP = {
    ("clinvar", "rs671[VRID]"): (1, ["18390"]),
    ("clinvar", "RCV000007523"): (1, ["7105"]),
    ("clinvar", "RCV000000019"): (1, ["9"]),
    ("clinvar", "VCV000017865[VACC]"): (1, ["17865"]),
    ("clinvar", "rs1800562[VRID]"): (1, ["9"]),
    ("clinvar", '"HFE"[gene] AND ("clinsig pathogenic"[prop] OR "clinsig likely pathogenic"[prop])'): (57, ["9"]),
    ("clinvar", '"hereditary hemochromatosis"[dis]'): (1220, ["9", "15333"]),
    ("clinvar", "12[chr] AND 111803962[chrpos38]"): (1, ["18390"]),
    ("clinvar", "17865[VID]"): (1, []),
    ("gene", '"HFE"[sym] AND "Homo sapiens"[orgn]'): (1, ["3077"]),
}


def _esummary_subset(src, ids):
    res = {"uids": []}
    for i in ids:
        if i in src["result"]:
            res["uids"].append(i)
            res[i] = src["result"][i]
    return {"header": {"type": "esummary"}, "result": res}


class FakeTransport:
    def __init__(self, fail=False):
        self.calls = []
        self.fail = fail

    def __call__(self, url, params, timeout):
        self.calls.append((url, dict(params)))
        if self.fail:
            raise httpx.ConnectError("simulated outage")
        assert params.get("tool") == "genomescope_wgs_explorer" or "eutils" not in url
        if url.endswith("esearch.fcgi"):
            count, ids = ESEARCH_MAP.get((params["db"], params["term"]), (0, []))
            body = {"esearchresult": {"count": str(count), "idlist": ids[: int(params.get("retmax", 20))], "querytranslation": params["term"]}}
            return (200, json.dumps(body), {})
        if url.endswith("esummary.fcgi"):
            ids = params["id"].split(",")
            src = {"clinvar": CLINVAR, "snp": SNP, "gene": GENE, "pubmed": PUBMED}[params["db"]]
            return (200, json.dumps(_esummary_subset(src, ids)), {})
        if url.endswith("elink.fcgi"):
            if params["dbfrom"] == "clinvar" and "18390" in params["id"].split(","):
                return (200, json.dumps(ELINK), {})
            return (200, json.dumps({"linksets": [{"dbfrom": params["dbfrom"], "ids": params["id"].split(",")}]}), {})
        if url.endswith("efetch.fcgi"):
            a, b = int(params["seq_start"]), int(params["seq_stop"])
            if params["id"] == "NC_000007.14" and FASTA_START <= a <= b < FASTA_START + len(FASTA_SEQ):
                return (200, ">NC_000007.14 slice\n" + FASTA_SEQ[a - FASTA_START: b - FASTA_START + 1] + "\n", {})
            return (400, "bad range", {})
        if "/refsnp/334" in url:
            return (200, json.dumps(REFSNP_334), {})
        return (404, "{}", {})


@pytest.fixture
def temp_db(job_ws):
    """The catalog of the fresh synthetic analysis (a per-test copy)."""
    import sqlite3

    db = job_ws.db_path
    conn = sqlite3.connect(str(db))
    conn.execute("DROP TABLE IF EXISTS ncbi_cache")  # start each test with an empty NCBI cache
    conn.commit()
    conn.close()
    # start from a migrated schema so row-immutability tests do not see the in-place migration
    indexer.ensure_ncbi_columns(force=True)
    return db


@pytest.fixture
def fake_ncbi(temp_db, monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    client = NcbiClient()
    client.limiter = RateLimiter(1000)
    client.backoff_base = 0.0
    fake = FakeTransport()
    monkeypatch.setattr(client, "_transport", fake)
    ncbi_client.set_client(client)
    yield fake
    ncbi_client.set_client(None)


APP_HEADERS = {"X-Requested-With": "GenomeScope"}


# ---------------------------------------------------------------------------
# Pure parser tests
# ---------------------------------------------------------------------------

def test_parse_esearch_and_esummary_envelopes():
    es = parse_esearch({"esearchresult": {"count": "2", "idlist": ["18390", "15077"], "querytranslation": "671[VRID]"}})
    assert es == {"count": 2, "ids": ["18390", "15077"], "translated": "671[VRID]"}
    with pytest.raises(NcbiError):
        parse_esearch({"esearchresult": {"ERROR": "Invalid query"}})
    recs = parse_esummary_result({"result": {"uids": ["1", "2"], "1": {"uid": "1"}, "2": {"uid": "2", "error": "cannot get document summary"}}})
    assert list(recs) == ["1"]


def test_parse_clinvar_rs671_aldh2():
    rec = parse_clinvar_summary(cv("18390"))
    assert rec["variation_id"] == "18390"
    assert rec["vcv_id"] == "VCV000018390"
    assert rec["title"].startswith("NM_000690.4(ALDH2):c.1510G>A")
    assert rec["chr"] == "12" and rec["start"] == 111803962 and rec["stop"] == 111803962
    assert rec["assembly"] == "GRCh38"
    assert (rec["ref"], rec["alt"]) == ("G", "A")
    assert rec["rsid"] == "rs671"
    assert rec["classification"] == "drug response"
    assert rec["classification_normalized"] == "Drug Response"
    assert rec["review_status"] == "reviewed by expert panel" and rec["review_stars"] == 3
    assert rec["last_evaluated"] == "2021-03-24"
    assert "RCV000020060" in rec["rcv_accessions"]
    assert rec["genes"][0] == {"symbol": "ALDH2", "gene_id": "217"}
    assert rec["molecular_consequences"] == ["missense variant"]
    assert rec["url"] == "https://www.ncbi.nlm.nih.gov/clinvar/variation/18390/"
    assert set(rec) == {
        "variation_id", "vcv_id", "title", "variant_type", "rcv_accessions", "classification",
        "classification_normalized", "review_status", "review_stars", "last_evaluated", "conditions",
        "genes", "chr", "start", "stop", "assembly", "ref", "alt", "rsid", "molecular_consequences",
        "protein_change", "url", "somatic_classification",
    }
    assert rec["somatic_classification"] is None


def test_parse_clinvar_conditions_and_stars_hfe():
    rec = parse_clinvar_summary(cv("9"))
    assert rec["classification_normalized"] == "Pathogenic / Risk Factor"
    assert rec["review_stars"] == 2
    hemo = rec["conditions"][0]
    assert hemo["name"] == "Hemochromatosis type 1"
    assert hemo["medgen_id"] == "C3469186"
    assert hemo["omim_ids"] == ["235200"]
    assert hemo["url"] == "https://www.ncbi.nlm.nih.gov/medgen/C3469186"
    assert rec["rsid"] == "rs1800562"


def test_review_stars_and_dates():
    assert review_status_to_stars("practice guideline") == 4
    assert review_status_to_stars("reviewed by expert panel") == 3
    assert review_status_to_stars("criteria provided, multiple submitters, no conflicts") == 2
    assert review_status_to_stars("criteria provided, conflicting classifications") == 1
    assert review_status_to_stars("criteria provided, single submitter") == 1
    assert review_status_to_stars("no assertion criteria provided") == 0
    assert normalize_ncbi_date("2021/03/24 00:00") == "2021-03-24"
    assert normalize_ncbi_date("1/01/01 00:00") is None
    assert normalize_ncbi_date("") is None


def test_spdi_to_vcf_snv_and_indel_anchor():
    assert spdi_to_vcf(111803961, "G", "A") == (111803962, "G", "A", None)
    # CFTR F508del: NC_000007.14:117559590:TCTT:T -> left-aligned VCF chr7:117559590 ATCT>A
    sp = extract_clinvar_spdi(cv("7105"))
    assert sp == ("NC_000007.14", 117559590, "TCTT", "T")
    pos, ref, alt, need = spdi_to_vcf(sp[1], sp[2], sp[3])
    assert need == 117559590 and (pos, ref, alt) == (117559591, "TCTT", "T")  # un-normalised fallback
    anchor = FASTA_SEQ[117559590 - FASTA_START]
    assert spdi_to_vcf(sp[1], sp[2], sp[3], anchor_base=anchor) == (117559590, "ATCT", "A", None)
    # Pure insertion without shared base and without anchor -> unknown
    assert spdi_to_vcf(100, "", "AG")[:3] == (None, None, None)
    # parse_clinvar_summary leaves indel ref/alt for the service to anchor
    rec = parse_clinvar_summary(cv("7105"))
    assert rec["ref"] is None and rec["alt"] is None and rec["start"] == 117559591


def test_alleles_equivalent_indel_representations():
    # SPDI (fully justified) vs left-aligned VCF with anchor base
    assert alleles_equivalent(117559590, "TCTT", "T", 117559590, "ATCT", "A") is True
    assert alleles_equivalent(117559590, "TCTT", "T", 117559590, "ATCT", "AT") is False
    # rs59826218: dbSNP 157480041:TGT:T vs personal VCF 157480041 CTG>C
    assert alleles_equivalent(157480041, "TGT", "T", 157480041, "CTG", "C") is True
    # SNV
    assert alleles_equivalent(111803961, "G", "A", 111803962, "G", "A") is True
    assert alleles_equivalent(111803961, "G", "A", 111803962, "G", "T") is False
    # Non-overlapping with a gap: undecidable
    assert alleles_equivalent(1000, "A", "G", 990, "C", "T") is None
    # Symbolic / missing ALT
    assert alleles_equivalent(1000, "A", "G", 1001, "A", ".") is False


def test_parse_dbsnp_and_global_maf():
    assert parse_global_maf({"study": "GnomAD_genomes", "freq": "A=0.0079314/1183"}) == {"study": "GnomAD_genomes", "allele": "A", "frequency": 0.0079314, "count": 1183}
    assert parse_global_maf({"study": "TWINSUK", "freq": "A=0./0"}) == {"study": "TWINSUK", "allele": "A", "frequency": 0.0, "count": 0}
    assert parse_global_maf({"study": "X", "freq": "garbage"}) is None
    d = parse_dbsnp_summary(snp("671"))
    assert d["rsid"] == "rs671" and d["chr"] == "12" and d["position"] == 111803962
    assert d["ref"] == "G" and d["alts"] == ["A"]
    assert d["spdi"] == "NC_000012.12:111803961:G:A"
    assert d["variant_class"] == "snv"
    assert "missense_variant" in d["function_class"]
    assert "drug-response" in d["clinical_significance"]
    assert d["genes"] == [{"symbol": "ALDH2", "gene_id": "217"}]
    assert d["preferred_frequency"] == {"study": "ALFA", "allele": "A", "frequency": 0.0060531}
    assert d["url"] == "https://www.ncbi.nlm.nih.gov/snp/rs671"
    multi = parse_dbsnp_summary(snp("1800562"))
    assert multi["alts"] == ["A", "T"]
    indel = parse_dbsnp_summary(snp("113993960"))
    assert indel["ref"] == "TCTT" and indel["alts"] == ["T", "TCTTCTT"] and indel["position"] == 117559591


def test_risk_allele_frequency_conversions():
    rs671 = parse_dbsnp_summary(snp("671"))
    target = ("NC_000012.12", 111803961, "G", "A")
    assert risk_allele_frequency_from_dbsnp(rs671, target, "A") == {"allele": "A", "frequency": 0.0060531, "study": "ALFA"}

    # Reported allele is the REF on a biallelic site -> 1 - f
    ref_reported = dict(rs671, frequencies=[{"study": "ALFA", "allele": "G", "frequency": 0.9, "count": 10}])
    rf = risk_allele_frequency_from_dbsnp(ref_reported, target, "A")
    assert rf["study"] == "ALFA" and abs(rf["frequency"] - 0.1) < 1e-9

    # Multi-allelic site reporting the REF: cannot convert -> skipped
    hfe = parse_dbsnp_summary(snp("1800562"))
    hfe_ref = dict(hfe, frequencies=[{"study": "ALFA", "allele": "G", "frequency": 0.9, "count": 10}])
    assert risk_allele_frequency_from_dbsnp(hfe_ref, ("NC_000006.12", 26092912, "G", "A"), "A") is None
    # ... but the real record reports the ClinVar ALT 'A'
    assert risk_allele_frequency_from_dbsnp(hfe, ("NC_000006.12", 26092912, "G", "A"), "A")["frequency"] == 0.052719
    # ClinVar ALT not among the rs alleles -> None
    assert risk_allele_frequency_from_dbsnp(hfe, ("NC_000006.12", 26092912, "G", "C"), "C") is None

    # rs334 (HBB, multi-allelic): KOREAN reports 'G' (another ALT) and must not be used for T>A
    hbb = parse_dbsnp_summary(snp("334"))
    rf = risk_allele_frequency_from_dbsnp(hbb, ("NC_000011.10", 5227001, "T", "A"), "A")
    assert rf == {"allele": "A", "frequency": 0.000785, "study": "ALFA"}
    hbb_g = dict(hbb, frequencies=[f for f in hbb["frequencies"] if f["study"] in ("KOREAN",)])
    assert risk_allele_frequency_from_dbsnp(hbb_g, ("NC_000011.10", 5227001, "T", "A"), "A") is None


def test_risk_allele_frequency_indel_naming():
    f508 = parse_dbsnp_summary(snp("113993960"))
    target = ("NC_000007.14", 117559590, "TCTT", "T")
    ref_names, alt_names = allele_names(*target[1:])
    assert "-" in alt_names and "T" in alt_names
    # ALFA reports the SPDI inserted sequence 'T'
    assert risk_allele_frequency_from_dbsnp(f508, target, "A") == {"allele": "A", "frequency": 0.003121, "study": "ALFA"}
    # other studies report '-' for the deletion
    no_alfa = dict(f508, frequencies=[f for f in f508["frequencies"] if f["study"] != "ALFA"])
    assert risk_allele_frequency_from_dbsnp(no_alfa, target, "A") == {"allele": "A", "frequency": 0.0079298, "study": "GnomAD_genomes"}


def test_refsnp_exact_frequencies():
    freqs = parse_refsnp_frequencies(REFSNP_334)
    assert any(f["study"] == "ALFA" for f in freqs)  # dbGaP_PopFreq is ALFA
    rf = risk_allele_frequency_from_refsnp(freqs, ("NC_000011.10", 5227001, "T", "A"), "A")
    assert rf["study"] == "ALFA" and abs(rf["frequency"] - 35 / 44606) < 1e-6
    assert risk_allele_frequency_from_refsnp(freqs, ("NC_000011.10", 5227001, "T", "C"), "C") is None


def test_parse_gene_and_pubmed():
    g = parse_gene_summary(GENE["result"]["217"])
    assert g["gene_id"] == "217" and g["symbol"] == "ALDH2"
    assert g["chr"] == "12" and g["map_location"] == "12q24.12"
    assert g["start"] == 111766933 and g["stop"] == 111817532
    assert "ALDM" in g["aliases"] and g["mim_ids"] == ["100650"]
    assert g["url"] == "https://www.ncbi.nlm.nih.gov/gene/217"
    assert g["summary"]
    uid = PUBMED["result"]["uids"][0]
    p = parse_pubmed_summary(PUBMED["result"][uid])
    assert p["pmid"] == uid and p["url"] == f"https://pubmed.ncbi.nlm.nih.gov/{uid}/"
    assert len(p["year"]) == 4 and p["title"] and p["journal"]


def test_translate_and_sanitize_search_terms():
    assert translate_search_term("RCV000019449") == [("rcv", "RCV000019449")]
    assert translate_search_term("vcv17865") == [("vcv", "VCV000017865[VACC]")]
    assert translate_search_term("rs671") == [("rsid", "rs671[VRID]")]
    assert translate_search_term("chr12:111,803,962") == [("locus", "12[chr] AND 111803962[chrpos38]")]
    assert translate_search_term("HFE")[0] == ("gene", '"HFE"[gene]')
    assert translate_search_term("CYP2C19")[0][0] == "gene"
    kinds = translate_search_term("hereditary hemochromatosis")
    assert kinds[0] == ("condition", '"hereditary hemochromatosis"[dis]')
    assert kinds[1] == ("text", '"hereditary"[All Fields] AND "hemochromatosis"[All Fields]')
    # Injection attempt: field tags / quotes / boolean syntax are neutralised
    evil = translate_search_term('BRCA1[gene] OR "x"[uid] #1 & 2')
    for _kind, q in evil:
        inner = q.replace('"[dis]', "").replace('"[All Fields]', "").replace('"', "")
        assert "[gene]" not in inner and "[uid]" not in inner and "#" not in inner and "&" not in inner
    assert sanitize_phrase('a"b[c]d') == "a b c d"
    assert translate_search_term("   ") == []


# ---------------------------------------------------------------------------
# Rate limiter, retry, cache
# ---------------------------------------------------------------------------

def test_rate_limiter_spacing_thread_safe():
    rl = RateLimiter(20)  # 50 ms interval
    stamps = []
    lock = threading.Lock()

    def worker():
        rl.acquire()
        with lock:
            stamps.append(time.monotonic())

    threads = [threading.Thread(target=worker) for _ in range(6)]
    t0 = time.monotonic()
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert time.monotonic() - t0 >= 0.24  # 5 intervals after the first slot
    stamps.sort()
    gaps = [b - a for a, b in zip(stamps, stamps[1:])]
    # per-gap floor is loose: time.monotonic() ticks in ~15.6 ms steps on Windows
    assert min(gaps) >= 0.02


def test_default_rate_limits(monkeypatch, temp_db):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    assert NcbiClient().rate_limit_per_sec == 3
    monkeypatch.setenv("NCBI_API_KEY", "dummy-key")
    c = NcbiClient()
    assert c.rate_limit_per_sec == 10
    sent = c._with_identity(ncbi_client.EUTILS_BASE + "esearch.fcgi", {"db": "clinvar"})
    assert sent["api_key"] == "dummy-key" and sent["tool"] == "genomescope_wgs_explorer"
    assert "email" not in sent  # only when NCBI_EMAIL is set
    monkeypatch.setenv("NCBI_EMAIL", "someone@example.org")
    assert NcbiClient()._with_identity(ncbi_client.EUTILS_BASE + "x", {})["email"] == "someone@example.org"
    # cache key never contains credentials
    assert "dummy-key" not in NcbiClient.cache_key("u", sent)


def test_retry_backoff_and_persistent_cache(temp_db, monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    c = NcbiClient()
    c.limiter = RateLimiter(1000)
    c.backoff_base = 0.0
    responses = [(503, "busy", {}), (429, "slow down", {"Retry-After": "0"}), (200, json.dumps({"esearchresult": {"count": "1", "idlist": ["9"]}}), {})]
    calls = []

    def transport(url, params, timeout):
        calls.append(url)
        return responses[min(len(calls) - 1, len(responses) - 1)]

    monkeypatch.setattr(c, "_transport", transport)
    assert c.esearch("clinvar", "x")["ids"] == ["9"]
    assert len(calls) == 3
    # second identical call is served from the SQLite cache
    assert c.esearch("clinvar", "x")["ids"] == ["9"]
    assert len(calls) == 3
    assert c.cache.count() >= 1
    # a fresh client (new process) still hits the persistent cache
    c2 = NcbiClient()
    monkeypatch.setattr(c2, "_transport", transport)
    assert c2.esearch("clinvar", "x")["ids"] == ["9"] and len(calls) == 3
    # TTL expiry forces a refetch
    c2.cache.ttl_seconds = -1
    c2.limiter = RateLimiter(1000)
    c2.esearch("clinvar", "x")
    assert len(calls) == 4


def test_client_errors_are_not_retried_and_failures_raise(temp_db, monkeypatch):
    c = NcbiClient()
    c.limiter = RateLimiter(1000)
    c.backoff_base = 0.0
    calls = []

    def bad_request(url, params, timeout):
        calls.append(1)
        return (400, "bad", {})

    monkeypatch.setattr(c, "_transport", bad_request)
    with pytest.raises(NcbiError):
        c.esearch("clinvar", "y")
    assert len(calls) == 1

    def down(url, params, timeout):
        calls.append(1)
        raise httpx.ConnectTimeout("timeout")

    monkeypatch.setattr(c, "_transport", down)
    with pytest.raises(NcbiError):
        c.esearch("clinvar", "z")
    assert len(calls) == 1 + ncbi_client.MAX_ATTEMPTS


# ---------------------------------------------------------------------------
# YourGenotype against the synthetic gVCF
# ---------------------------------------------------------------------------

def test_your_genotype_explicit_alt_hom_rs671():
    gt = ncbi_service.build_your_genotype("12", 111803962, target_spdi=("NC_000012.12", 111803961, "G", "A"))
    assert gt["found"] and gt["genotype"] == "1/1"
    assert gt["carries_ncbi_alt"] is True and gt["zygosity_for_ncbi_alt"] == "homozygous"
    assert gt["is_reference_block"] is False and gt["dp"] == 30 and gt["gq"] == 90 and gt["ad"] == "0,30"
    assert gt["raw_vcf_line"].startswith("12\t111803962\trs671")
    # A different ALT at the same position is not the NCBI allele
    other = ncbi_service.build_your_genotype("12", 111803962, target_spdi=("NC_000012.12", 111803961, "G", "T"))
    assert other["carries_ncbi_alt"] is False and other["zygosity_for_ncbi_alt"] == "absent"
    assert other["is_carried"] is True and other["your_status"] == ncbi_service.DIFFERENT_ALT_STATUS


def test_your_genotype_reference_block_and_no_call():
    hfe = ncbi_service.build_your_genotype("chr6", 26092913, target_spdi=("NC_000006.12", 26092912, "G", "A"))
    assert hfe["found"] and hfe["is_reference_block"] is True
    assert hfe["carries_ncbi_alt"] is False and hfe["zygosity_for_ncbi_alt"] == "absent"
    nc = ncbi_service.build_your_genotype("1", 10050, target_spdi=(None, 10049, "A", "G"))
    assert nc["found"] and nc["genotype"] == "./."
    assert nc["zygosity_for_ncbi_alt"] == "no_call" and nc["carries_ncbi_alt"] is None


def test_your_genotype_indel_spdi_vs_left_aligned_vcf():
    # dbSNP rs59826218 SPDI 157480041:TGT:T; the synthetic gVCF has 157480041 CTG>C 0/1
    gt = ncbi_service.build_your_genotype("7", 157480042, target_spdi=("NC_000007.14", 157480041, "TGT", "T"))
    assert gt["found"] and gt["position"] == 157480041 and gt["ref"] == "CTG" and gt["alt"] == "C"
    assert gt["carries_ncbi_alt"] is True and gt["zygosity_for_ncbi_alt"] == "heterozygous"
    # The other dbSNP allele (insertion) is not carried
    ins = ncbi_service.build_your_genotype("7", 157480042, target_spdi=("NC_000007.14", 157480041, "TGT", "TGTGT"))
    assert ins["carries_ncbi_alt"] is False


def test_your_genotype_not_found():
    gt = ncbi_service.build_your_genotype("Y", 999_999_999, target_spdi=(None, 999_999_998, "A", "G"))
    assert gt["found"] is False and gt["zygosity_for_ncbi_alt"] == "unknown" and gt["carries_ncbi_alt"] is None


# ---------------------------------------------------------------------------
# Endpoints (fixture transport, temp DB)
# ---------------------------------------------------------------------------

VARIANT_KEYS = {"query", "resolved", "clinvar", "dbsnp", "gene", "publications", "your_genotype", "risk_allele_frequency", "links", "errors", "fetched_at"}
GENOTYPE_KEYS = {
    "found", "chr", "position", "ref", "alt", "genotype", "alleles_called", "your_status", "is_carried", "dp", "gq", "ad",
    "filter", "is_reference_block", "raw_vcf_line", "source_file", "carries_ncbi_alt", "zygosity_for_ncbi_alt",
}
DBSNP_KEYS = {"rsid", "chr", "position", "ref", "alts", "spdi", "variant_class", "function_class", "clinical_significance", "genes", "frequencies", "preferred_frequency", "url"}
GENE_KEYS = {"gene_id", "symbol", "name", "summary", "aliases", "chr", "start", "stop", "map_location", "mim_ids", "url"}


def test_endpoint_variant_rs671(api, fake_ncbi):
    r = api.get("/api/ncbi/variant", params={"rsid": "rs671"})
    assert r.status_code == 200
    b = r.json()
    assert set(b) == VARIANT_KEYS
    assert b["errors"] == [] and b["fetched_at"].endswith("Z")
    assert b["resolved"] == {"chr": "12", "position": 111803962, "rsid": "rs671"}
    assert b["clinvar"][0]["variation_id"] == "18390"
    assert set(b["dbsnp"]) == DBSNP_KEYS and set(b["gene"]) == GENE_KEYS
    assert set(b["your_genotype"]) == GENOTYPE_KEYS
    assert b["your_genotype"]["carries_ncbi_alt"] is True and b["your_genotype"]["zygosity_for_ncbi_alt"] == "homozygous"
    assert b["risk_allele_frequency"] == {"allele": "A", "frequency": 0.0060531, "study": "ALFA"}
    assert 1 <= len(b["publications"]) <= 10
    assert set(b["publications"][0]) == {"pmid", "title", "journal", "year", "authors", "url"}
    assert b["links"]["clinvar"].endswith("/clinvar/variation/18390/")
    assert b["links"]["dbsnp"].endswith("/snp/rs671") and b["links"]["gene"].endswith("/gene/217")
    assert b["links"]["pubmed"]
    # locus query resolves to the same record
    r2 = api.get("/api/ncbi/variant", params={"chr": "chr12", "position": "111803962", "include_publications": "false"}).json()
    assert r2["clinvar"][0]["variation_id"] == "18390" and r2["publications"] == []


def test_endpoint_variant_rcv_indel_anchor(api, fake_ncbi):
    b = api.get("/api/ncbi/variant", params={"rcv_id": "RCV000007523"}).json()
    rec = b["clinvar"][0]
    assert rec["variation_id"] == "7105" and (rec["ref"], rec["alt"]) == ("ATCT", "A")
    assert b["your_genotype"]["zygosity_for_ncbi_alt"] == "absent"
    assert b["risk_allele_frequency"]["frequency"] == 0.003121
    assert any(u.endswith("efetch.fcgi") for u, _p in fake_ncbi.calls)


def test_endpoint_variant_validation(api, fake_ncbi):
    assert api.get("/api/ncbi/variant").status_code == 400
    assert api.get("/api/ncbi/variant", params={"chr": "12"}).status_code == 400
    assert api.get("/api/ncbi/variant", params={"rsid": "rs671; DROP"}).status_code == 400
    assert api.get("/api/ncbi/gene/HFE", params={"significance": "bogus"}).status_code == 400
    assert api.get("/api/ncbi/gene/HFE", params={"limit": 0}).status_code == 422
    assert api.get("/api/ncbi/search", params={"term": "x", "limit": 101}).status_code == 422
    assert api.post("/api/ncbi/import", headers=APP_HEADERS, json={"variation_ids": []}).status_code == 400
    assert api.post("/api/ncbi/import", headers=APP_HEADERS, json={"variation_ids": [str(i) for i in range(51)]}).status_code == 400


def test_endpoint_gene_and_search(api, fake_ncbi):
    g = api.get("/api/ncbi/gene/HFE", params={"significance": "pathogenic", "limit": 5}).json()
    assert set(g) == {"gene", "significance", "clinvar_total", "variants", "carried_count", "errors", "fetched_at"}
    assert g["gene"]["symbol"] == "HFE" and g["clinvar_total"] == 57 and g["significance"] == "pathogenic"
    v = g["variants"][0]
    assert set(v) == {"clinvar", "your_genotype"}
    assert v["clinvar"]["variation_id"] == "9" and v["your_genotype"]["zygosity_for_ncbi_alt"] == "absent"
    assert g["carried_count"] == 0

    s = api.get("/api/ncbi/search", params={"term": "VCV000017865"}).json()
    assert set(s) == {"term", "translated_query", "total", "items", "errors", "fetched_at"}
    assert s["total"] == 1 and s["items"][0]["clinvar"]["variation_id"] == "17865"
    assert s["translated_query"] == "VCV000017865[VACC]"

    s2 = api.get("/api/ncbi/search", params={"term": "hereditary hemochromatosis", "limit": 1}).json()
    assert s2["total"] == 1220 and len(s2["items"]) == 1

    s3 = api.get("/api/ncbi/search", params={"term": 'zz"[uid] OR 1'}).json()
    assert s3["total"] == 0 and s3["items"] == []
    sent_terms = [p.get("term", "") for _u, p in fake_ncbi.calls if p.get("db") == "clinvar"]
    assert not any("[uid]" in t for t in sent_terms)


def test_endpoint_import_and_catalog_sync(api, fake_ncbi, temp_db):
    imp = api.post("/api/ncbi/import", headers=APP_HEADERS, json={"variation_ids": ["18390", "VCV000000009", "99999999"]}).json()
    assert set(imp) == {"imported", "record_ids", "skipped", "errors", "fetched_at"}
    assert imp["imported"] == 2 and len(imp["record_ids"]) == 2
    assert imp["skipped"] == [{"variation_id": "99999999", "reason": "not found in NCBI ClinVar"}]

    found = api.get("/api/variants", params={"variant_id": "rs1800562"}).json()
    assert found["total"] >= 1
    row = next(i for i in found["items"] if i["id"] in imp["record_ids"])
    assert row["ncbi_variation_id"] == "9" and "NCBI ClinVar" in row["evidence"]
    assert row["is_carried"] == 0 and row["frequency_source"] == "NCBI dbSNP ALFA"

    rs671_id = api.get("/api/variants", params={"variant_id": "rs671"}).json()["items"][0]["id"]
    assert rs671_id in imp["record_ids"]  # existing catalog row updated in place, id preserved

    r = api.get(f"/api/variants/{rs671_id}/ncbi").json()
    assert VARIANT_KEYS <= set(r) and r["record_id"] == rs671_id and r["synced"] is False
    fields = {c["field"]: c for c in r["catalog_vs_ncbi"]}
    assert set(fields) == {"allele", "classification", "review_status", "condition", "rcv_id", "risk_frequency"}
    assert fields["allele"]["match"] is True
    assert fields["classification"]["ncbi"] == "Drug Response"

    assert api.get(f"/api/variants/{rs671_id}/ncbi", params={"sync": "true"}).status_code == 405  # GET is read-only
    assert api.post(f"/api/variants/{rs671_id}/ncbi/sync").status_code == 403  # app header required
    s = api.post(f"/api/variants/{rs671_id}/ncbi/sync", headers=APP_HEADERS).json()
    assert s["synced"] is True
    fields = {c["field"]: c for c in s["catalog_vs_ncbi"]}
    assert fields["classification"]["match"] and fields["review_status"]["match"] and fields["risk_frequency"]["match"]
    detail = api.get(f"/api/variants/{rs671_id}").json()
    assert detail["classification"] == "Drug Response" and detail["ncbi_synced_at"]
    assert detail["risk_frequency"] == pytest.approx(0.0060531, rel=1e-3)

    assert api.get("/api/variants/999999999/ncbi").status_code == 404


def test_variants_search_falls_back_to_ncbi(api, fake_ncbi, monkeypatch):
    from app.api import explorer as main_mod

    monkeypatch.setattr(main_mod, "resolve_live_query", lambda **kw: 0)
    before = api.get("/api/variants", params={"rcv_id": "RCV000000019"}).json()
    assert before["total"] >= 1  # imported live from NCBI ClinVar then re-searched
    item = before["items"][0]
    assert "RCV000000019" in (item["all_rcv_ids"] or "") + (item["rcv_id"] or "")


def test_ncbi_outage_degrades_gracefully(api, temp_db, monkeypatch):
    c = NcbiClient()
    c.limiter = RateLimiter(1000)
    c.backoff_base = 0.0
    monkeypatch.setattr(c, "_transport", FakeTransport(fail=True))
    ncbi_client.set_client(c)
    try:
        r = api.get("/api/ncbi/variant", params={"rsid": "rs671"})
        assert r.status_code == 200
        b = r.json()
        assert b["errors"] and b["clinvar"] == [] and b["dbsnp"] is None
        if VCF_AVAILABLE:
            assert b["your_genotype"] is None or isinstance(b["your_genotype"], dict)
        g = api.get("/api/ncbi/gene/HFE").json()
        assert g["errors"] and g["variants"] == []
        st = api.get("/api/ncbi/status").json()
        assert st["reachable"] is False and st["errors"]
    finally:
        ncbi_client.set_client(None)


def test_status_and_health(api, fake_ncbi):
    st = api.get("/api/ncbi/status").json()
    assert set(st) == {"reachable", "api_key_configured", "rate_limit_per_sec", "cache_entries", "cache_ttl_hours", "services", "errors", "fetched_at"}
    assert st["reachable"] is True and st["rate_limit_per_sec"] == 3 and st["api_key_configured"] is False
    assert st["services"]["eutils"].startswith("https://eutils.ncbi.nlm.nih.gov/")
    h = api.get("/api/health").json()
    assert set(h["ncbi"]) == {"api_key_configured", "cache_entries"}


def test_schema_migration_adds_ncbi_columns(temp_db):
    indexer.init_db()
    conn = indexer.get_db_connection()
    cols = {r[1] for r in conn.execute("PRAGMA table_info(variants)").fetchall()}
    conn.close()
    assert {"ncbi_variation_id", "ncbi_last_evaluated", "ncbi_synced_at", "frequency_source"} <= cols
    assert {"clinvar_variation_id", "clinvar_release", "clinvar_match"} <= cols


# ---------------------------------------------------------------------------
# Live tests (opt-in): RUN_LIVE_NCBI=1 pytest tests/test_ncbi.py -k live  (from backend/)
# ---------------------------------------------------------------------------

live = pytest.mark.skipif(os.environ.get("RUN_LIVE_NCBI") != "1", reason="set RUN_LIVE_NCBI=1 to query real NCBI")


@pytest.fixture
def live_client(temp_db):
    ncbi_client.set_client(NcbiClient())
    yield
    ncbi_client.set_client(None)


@live
def test_live_rs671(api, live_client):
    b = api.get("/api/ncbi/variant", params={"rsid": "rs671"}).json()
    assert b["errors"] == []
    assert b["clinvar"] and b["clinvar"][0]["genes"][0]["symbol"] == "ALDH2"
    assert b["dbsnp"]["position"] == 111803962 and b["gene"]["symbol"] == "ALDH2"
    assert b["risk_allele_frequency"]["allele"] == "A"
    if VCF_AVAILABLE:
        assert b["your_genotype"]["zygosity_for_ncbi_alt"] == "homozygous"


@live
def test_live_rcv_and_gene(api, live_client):
    b = api.get("/api/ncbi/variant", params={"rcv_id": "RCV000019449", "include_publications": "false"}).json()
    assert b["clinvar"][0]["variation_id"] == "17865"
    g = api.get("/api/ncbi/gene/HFE", params={"significance": "pathogenic", "limit": 10}).json()
    assert g["gene"]["gene_id"] == "3077" and g["clinvar_total"] > 10 and len(g["variants"]) == 10
    s = api.get("/api/ncbi/status").json()
    assert s["reachable"] is True


def test_complex_records_are_not_confirmed_by_single_locus():
    import copy

    raw = copy.deepcopy(cv("18390"))
    raw["obj_type"] = "Haplotype"
    raw["variation_set"] = raw["variation_set"] * 2
    rec = parse_clinvar_summary(raw)
    gt = ncbi_service.genotype_for_clinvar(raw, rec, {"18390": (111803962, "G", "A")})
    assert gt["found"] and gt["genotype"] == "1/1"
    assert gt["carries_ncbi_alt"] is None and gt["zygosity_for_ncbi_alt"] == "unknown"
    # a simple record at the same locus is confirmed
    simple = ncbi_service.genotype_for_clinvar(cv("18390"), parse_clinvar_summary(cv("18390")), {"18390": (111803962, "G", "A")})
    assert simple["carries_ncbi_alt"] is True


def test_reference_allele_record_counts_reference_calls():
    # ClinVar "=" records (deleted == inserted) describe the reference allele itself
    gt = ncbi_service.build_your_genotype("6", 26092913, target_spdi=("NC_000006.12", 26092912, "G", "G"))
    # Reference copies are not "carrying an ALT": neutral 'reference_allele' state, never counted / ranked as carried
    assert gt["is_reference_block"] and gt["carries_ncbi_alt"] is None and gt["zygosity_for_ncbi_alt"] == "reference_allele"
    assert gt["ref"] == "G" and gt["alleles_called"] == "G / G"


def test_gene_scan_pulls_clinvar_records_at_users_own_loci(monkeypatch):
    """Large genes: the user's own non-reference loci are queried so carriers reach the ranking."""
    rows = [
        {"source_file": ncbi_service.SNP_VCF_NAME, "position": 100, "ref": "G", "alt": "A", "is_carried": True},
        {"source_file": ncbi_service.SNP_VCF_NAME, "position": 200, "ref": "ATCT", "alt": "A", "is_carried": True},
        {"source_file": ncbi_service.SNP_VCF_NAME, "position": 300, "ref": "C", "alt": "T", "is_carried": False},
        {"source_file": "sv.vcf", "position": 400, "ref": "N", "alt": "<DEL>", "is_carried": True},
    ]
    monkeypatch.setattr(vcf_engine_mod.active_engine(), "read_vcf_region", lambda *a, **k: rows)
    gene = {"chr": "13", "start": 50, "stop": 500}
    assert ncbi_service.user_carried_positions(gene) == ("13", [100, 200, 201])
    assert ncbi_service.user_carried_positions({"chr": "13", "start": None, "stop": None}) == (None, [])

    seen = []

    def fake_search(term, retmax, errors):
        seen.append(term)
        return (2, ["11", "22"], term)

    monkeypatch.setattr(ncbi_service, "search_clinvar_ids", fake_search)
    errors = []
    assert ncbi_service.search_user_loci('"BRCA2"[gene]', gene, errors) == ["11", "22"]
    assert seen == ['("BRCA2"[gene]) AND 13[chr] AND (100[chrpos38] OR 200[chrpos38] OR 201[chrpos38])']
    assert errors == []


# ---------------------------------------------------------------------------
# Regression tests for the reviewed findings
# ---------------------------------------------------------------------------

import copy as _copy  # noqa: E402

from app.api import explorer as main_mod  # noqa: E402
from genomescope.ncbi_client import (  # noqa: E402
    SIGNIFICANCE_FILTERS,
    classification_rank,
    normalize_clinvar_classification,
)


def _db_row(where, params=()):
    conn = indexer.get_db_connection()
    try:
        r = conn.execute(f"SELECT * FROM variants WHERE {where} ORDER BY id LIMIT 1", params).fetchone()
        return dict(r) if r else None
    finally:
        conn.close()


def _set_row(row_id, **fields):
    conn = indexer.get_db_connection()
    try:
        conn.execute(f"UPDATE variants SET {', '.join(k + ' = ?' for k in fields)} WHERE id = ?", list(fields.values()) + [row_id])
        conn.commit()
    finally:
        conn.close()


def _fake_item(uid, alt, rcvs, complex_=False, start=100, ref="G", cls="drug response", rsid_num="671"):
    raw = _copy.deepcopy(cv("18390"))
    raw["uid"] = uid
    raw["accession"] = f"VCV{int(uid):09d}"
    raw["supporting_submissions"] = {"rcv": list(rcvs), "scv": []}
    raw["germline_classification"] = dict(raw["germline_classification"], description=cls)
    vs = raw["variation_set"][0]
    vs["canonical_spdi"] = f"NC_000012.12:{start - 1}:{ref}:{alt}"
    vs["variation_xrefs"] = [{"db_source": "dbSNP", "db_id": rsid_num}]
    if complex_:
        raw["obj_type"] = "Diplotype"
        raw["variation_set"] = [vs, _copy.deepcopy(vs)]
    rec = parse_clinvar_summary(raw)
    return (raw, rec)


# -- upsert-replace-churns-ids-and-wipes-ncbi-sync / variants-fallback-imports-unrelated-records

def test_upsert_keeps_row_id_and_ncbi_sync(temp_db):
    row = _db_row("variant_id = 'rs671' AND source_file = ?", (ncbi_service.SNP_VCF_NAME,))
    _set_row(row["id"], ncbi_variation_id="18390", ncbi_synced_at="2026-01-01T00:00:00Z", classification="Drug Response", condition="Synced condition")
    rec = {c: row[c] for c in indexer.UPSERT_BASE_COLUMNS}
    rec.update(classification="Pathogenic", condition="MyVariant condition", your_status="Changed status")
    indexer.upsert_records([rec])
    after = _db_row("unique_key = ?", (row["unique_key"],))
    assert after["id"] == row["id"]  # no INSERT OR REPLACE id churn
    assert after["classification"] == "Drug Response" and after["condition"] == "Synced condition"
    assert after["ncbi_variation_id"] == "18390" and after["ncbi_synced_at"] == "2026-01-01T00:00:00Z"
    assert after["your_status"] == "Changed status"  # VCF-derived fields still refresh
    # a row never synced from NCBI is updated normally (still in place)
    _set_row(row["id"], ncbi_synced_at=None)
    indexer.upsert_records([rec])
    again = _db_row("unique_key = ?", (row["unique_key"],))
    assert again["id"] == row["id"] and again["classification"] == "Pathogenic" and again["ncbi_variation_id"] == "18390"


def test_variants_live_resolvers_only_for_absent_identifiers(api, temp_db, monkeypatch):
    calls = []
    monkeypatch.setattr(main_mod, "resolve_live_query", lambda **kw: calls.append(("mv", kw)) or 0)
    monkeypatch.setattr(ncbi_service, "fallback_import_for_search", lambda **kw: calls.append(("ncbi", kw)) or 0)
    before = _db_row("variant_id = 'rs671'")
    # identifier present, other filters empty the result: nothing resolved / written
    r = api.get("/api/variants", params={"variant_id": "rs671", "classification": "Benign"}).json()
    assert r["total"] == 0 and calls == []
    assert api.get("/api/variants", params={"gene": "ALDH2", "carried_only": "true", "classification": "Benign"}).json()["total"] == 0
    assert calls == []
    # free text never triggers live imports
    assert api.get("/api/variants", params={"q": "hypertrophic", "carried_only": "true"}).json()["total"] == 0
    assert calls == []
    assert _db_row("id = ?", (before["id"],)) == before
    # an identifier absent from the catalog is resolved (MyVariant first, then NCBI)
    api.get("/api/variants", params={"variant_id": "rs999999999991", "carried_only": "true"})
    assert [c[0] for c in calls] == ["mv", "ncbi"] and calls[1][1]["variant_id"] == "rs999999999991"
    assert "condition" not in calls[1][1] and "q" not in calls[1][1]


def test_fallback_import_ignores_free_text(fake_ncbi):
    assert ncbi_service.fallback_import_for_search() == 0
    assert ncbi_service.fallback_import_for_search(gene="not a gene!") == 0
    assert fake_ncbi.calls == []


# -- primary-allele-fallback-wrong-allele-sync / sync-writes-wrong-allele-record / primary-ignores-genotype

def test_select_primary_prefers_catalog_identity_then_carried_allele():
    simple_t = _fake_item("1703261", "T", ["RCV1"], cls="drug response")
    simple_a = _fake_item("16897", "A", ["RCV000018393"], cls="likely benign")
    diplo = _fake_item("634007", "A", ["RCV000782691"], complex_=True)
    items = sorted([simple_t, simple_a, diplo], key=ncbi_service._item_order)
    sel = ncbi_service._select_primary
    assert sel(items, prefer_rcvs=["RCV000018393"])[0][1]["variation_id"] == "16897"
    assert sel(items, prefer_rcvs=["RCV000018393"])[1] == "rcv"
    # catalog RCV belonging only to a diplotype -> that record (sync then refuses it)
    it, how = sel(items, prefer_rcvs=["RCV000782691"])
    assert it[1]["variation_id"] == "634007" and how == "rcv"
    # prefer_alt never picks a complex record
    assert sel([diplo, simple_a], prefer_alt="A")[0][1]["variation_id"] == "16897"
    # no catalog match: the user's carried allele beats severity ordering
    assert sel(items, carried_ids=["16897"]) == (simple_a, "carried")
    assert sel(items)[1] == "fallback" and sel(items)[0] is simple_t


def test_sync_refusal_rules():
    rec_ok = dict(parse_clinvar_summary(cv("18390")), variation_id="1", vcv_id="VCV000000001", ref="C", alt="A")
    body = {"clinvar": [rec_ok], "your_genotype": {"carries_ncbi_alt": False, "is_carried": False, "alt": "."}}
    row = {"rcv_id": "RCV1", "risk_version": "A"}
    assert ncbi_service._sync_refusal(row, body, {"how": "rcv"}) is None
    assert "haplotype" in ncbi_service._sync_refusal(row, body, {"how": "rcv", "primary_complex": True})
    assert "no NCBI ClinVar record matches" in ncbi_service._sync_refusal(row, body, {"how": "fallback"})
    assert "no NCBI ClinVar record matches" in ncbi_service._sync_refusal(row, body, {"how": "carried"})
    diff = dict(body, your_genotype={"carries_ncbi_alt": False, "is_carried": True, "alt": "T"})
    assert "different ALT" in ncbi_service._sync_refusal(row, diff, {"how": "allele"})
    cmp = {c["field"]: c for c in ncbi_service._comparisons({"risk_version": "A", "raw_vcf_line": "2\t1\t.\tC\tT\t.\tPASS\t.\tGT\t0/1"}, diff)}
    assert cmp["allele"]["match"] is False and "your VCF ALT T" in cmp["allele"]["catalog"]


def test_sync_refuses_catalog_row_for_another_allele(api, fake_ncbi, temp_db):
    row = _db_row("variant_id = 'rs671'")
    # catalog row describing a different (T) allele with an unrelated RCV
    _set_row(row["id"], rcv_id="RCV000000001", all_rcv_ids="RCV000000001", risk_version="T", ncbi_variation_id=None, classification="Benign")
    r = api.post(f"/api/variants/{row['id']}/ncbi/sync", headers=APP_HEADERS).json()
    assert r["synced"] is False and any(e.startswith("Sync skipped") for e in r["errors"])
    assert r["clinvar"][0]["variation_id"] == "18390"  # shown (the allele you carry) but not written
    assert _db_row("id = ?", (row["id"],))["classification"] == "Benign"


def test_rsid_lookup_fetches_all_records_and_picks_carried_simple(api, fake_ncbi, monkeypatch):
    ids = []
    for k in range(129):
        uid = str(900000 + k)
        raw, _rec = _fake_item(uid, "A", [f"RCV9{k:08d}"], complex_=True, start=111803962)
        monkeypatch.setitem(CLINVAR["result"], uid, raw)
        ids.append(uid)
    ids.append("18390")  # the simple record is the oldest -> last in esearch order
    monkeypatch.setitem(ESEARCH_MAP, ("clinvar", "rs671[VRID]"), (130, ids))
    b = api.get("/api/ncbi/variant", params={"rsid": "rs671", "include_publications": "false"}).json()
    assert b["clinvar"][0]["variation_id"] == "18390"
    assert b["your_genotype"]["carries_ncbi_alt"] is True and b["your_genotype"]["zygosity_for_ncbi_alt"] == "homozygous"
    assert len(b["clinvar"]) == 130 and not any("only the" in e for e in b["errors"])
    retmaxes = [int(p["retmax"]) for u, p in fake_ncbi.calls if u.endswith("esearch.fcgi") and p["term"] == "rs671[VRID]"]
    assert max(retmaxes) >= 130


# -- import-overwrites-different-clinvar-record-row

def test_import_never_overwrites_another_records_row(api, fake_ncbi, temp_db, monkeypatch):
    row = _db_row("variant_id = 'rs671'")
    _set_row(row["id"], ncbi_variation_id=None)
    hap, _r = _fake_item("441268", "A", ["RCV000999001"], complex_=True, start=111803962, cls="Pathogenic/Likely pathogenic")
    sib, _r2 = _fake_item("441269", "A", ["RCV000999002"], start=111803962, cls="Pathogenic")
    monkeypatch.setitem(CLINVAR["result"], "441268", hap)
    monkeypatch.setitem(CLINVAR["result"], "441269", sib)
    before = _db_row("id = ?", (row["id"],))
    imp = api.post("/api/ncbi/import", headers=APP_HEADERS, json={"variation_ids": ["441268", "441269"]}).json()
    assert imp["imported"] == 2 and row["id"] not in imp["record_ids"]
    assert _db_row("id = ?", (row["id"],)) == before  # the user's rs671 row is untouched
    hap_row = _db_row("ncbi_variation_id = '441268'")
    assert hap_row["is_carried"] == 0  # complex record: carriage unknown -> not carried
    assert hap_row["classification"] == "Pathogenic/Likely Pathogenic"
    # a record whose RCV matches the catalog row does update it in place
    sib_same = _copy.deepcopy(sib)
    sib_same["uid"] = "441270"
    sib_same["supporting_submissions"] = {"rcv": [row["rcv_id"]], "scv": []}
    monkeypatch.setitem(CLINVAR["result"], "441270", sib_same)
    imp2 = api.post("/api/ncbi/import", headers=APP_HEADERS, json={"variation_ids": ["441270"]}).json()
    assert imp2["record_ids"] == [row["id"]]


# -- split-row-snv-missed / filter-ignored-in-carries / hemizygous-called-two-copies

def test_genotype_reads_every_row_at_split_sites():
    # synthetic 1:236414462: rs61021118 C>CT 0/1 followed by rs625443 C>T PASS 0/1
    gt = ncbi_service.build_your_genotype("1", 236414462, target_spdi=("NC_000001.11", 236414461, "C", "T"))
    assert gt["carries_ncbi_alt"] is True and gt["zygosity_for_ncbi_alt"] == "heterozygous"
    assert "\trs625443\t" in gt["raw_vcf_line"] and gt["alleles_called"] == "C / T"
    # rs4364901: a FAIL C>CAT 1/1 row must not shadow the PASS C>T 0/1 call
    gt2 = ncbi_service.build_your_genotype("1", 100186848, target_spdi=("NC_000001.11", 100186847, "C", "T"))
    assert gt2["carries_ncbi_alt"] is True and gt2["filter"] == "PASS"
    # 1-bp deletion anchored one base earlier: the rows at the deleted base are read
    gt3 = ncbi_service.build_your_genotype("2", 47414420, target_spdi=(None, 47414419, "T", ""), target_vcf=(47414419, "GT", "G"))
    assert gt3["raw_vcf_line"].startswith("2\t47414420\t") and gt3["is_reference_block"] is False


def test_fail_calls_are_low_quality_not_carried():
    gt = ncbi_service.build_your_genotype("17", 46337318, target_spdi=("NC_000017.11", 46337317, "T", "C"))
    assert gt["filter"] == "FAIL" and gt["genotype"] == "1/1"
    assert gt["carries_ncbi_alt"] is None and gt["zygosity_for_ncbi_alt"] == "low_quality"


def test_haploid_calls_are_hemizygous():
    gt = ncbi_service.build_your_genotype("X", 150596742, target_spdi=("NC_000023.11", 150596741, "T", "C"))
    assert gt["genotype"] == "1" and gt["carries_ncbi_alt"] is True
    assert gt["zygosity_for_ncbi_alt"] == "hemizygous"


# -- ref-allele-records-shown-as-carried

def test_reference_allele_records_rank_after_allele_changes():
    eq_item = _fake_item("226007", "G", ["RCV2"], ref="G", cls="Conflicting classifications of pathogenicity")
    real = _fake_item("642", "A", ["RCV1"], cls="Pathogenic")
    eq_gt = {"carries_ncbi_alt": None, "zygosity_for_ncbi_alt": "reference_allele"}
    real_gt = {"carries_ncbi_alt": False, "zygosity_for_ncbi_alt": "absent"}
    ranked = sorted([(eq_item, eq_gt), (real, real_gt)], key=ncbi_service._rank_with_genotype)
    assert ranked[0][0] is real


# -- classification-normalization-misleading / plp-downgraded-to-lp / risk-filter-misses-risk-allele-terms

def test_classification_normalisation():
    norm = indexer.normalize_classification
    assert norm("Pathogenic/Likely pathogenic") == "Pathogenic/Likely Pathogenic"
    assert norm("Pathogenic/Likely pathogenic/Pathogenic, low penetrance/Established risk allele; risk factor") == "Pathogenic/Likely Pathogenic / Risk Factor"
    assert classification_rank(norm("Pathogenic/Likely pathogenic")) == 0
    assert norm("Benign/Likely benign") == "Benign/Likely Benign" and classification_rank("Benign/Likely Benign") == 7
    assert norm("Conflicting classifications of pathogenicity; other; risk factor") == "Conflicting Interpretations / Risk Factor"
    assert norm("Likely pathogenic") == "Likely Pathogenic" and norm("drug response") == "Drug Response"
    assert normalize_clinvar_classification("") == "Not Classified"
    somatic = _copy.deepcopy(cv("18390"))
    somatic["germline_classification"] = dict(somatic["germline_classification"], description="")
    somatic["oncogenicity_classification"] = {"description": "Likely oncogenic"}
    rec = parse_clinvar_summary(somatic)
    assert rec["classification_normalized"] == "Not Classified"
    assert rec["somatic_classification"] == "Oncogenicity: Likely oncogenic"


def test_risk_allele_terms_in_filter_rank_and_post_filter():
    assert "risk allele" in SIGNIFICANCE_FILTERS["risk"]
    assert classification_rank("Established Risk Allele") == 2 and classification_rank("Likely Risk Allele") == 2
    item = _fake_item("1", "A", ["RCV1"], cls="Established risk allele")
    assert ncbi_service._post_filter([item], "risk") == [item]


# -- cnv-sv-rows-queried-by-anchor-base

def test_structural_rows_skip_anchor_base_lookup(api, fake_ncbi, temp_db):
    row = _db_row("source_file LIKE '%cnv.vcf'")
    b = api.get(f"/api/variants/{row['id']}/ncbi").json()
    assert b["clinvar"] == [] and b["dbsnp"] is None and b["your_genotype"] is None
    assert b["errors"] and "structural" in b["errors"][0]
    assert all(c["match"] is None for c in b["catalog_vs_ncbi"])
    assert fake_ncbi.calls == []  # no lookup of the breakpoint base
    s = api.post(f"/api/variants/{row['id']}/ncbi/sync", headers=APP_HEADERS).json()
    assert s["synced"] is False


# -- cors-wildcard-credentials-and-get-sync

def test_cors_and_state_changing_routes(api, fake_ncbi):
    evil = api.options("/api/ncbi/import", headers={"Origin": "https://evil.example", "Access-Control-Request-Method": "POST"})
    assert "access-control-allow-origin" not in evil.headers
    simple = api.get("/api/ncbi/status", headers={"Origin": "https://evil.example"})
    assert "access-control-allow-origin" not in simple.headers
    ok = api.get("/api/ncbi/status", headers={"Origin": "http://localhost:5173"})
    assert ok.headers.get("access-control-allow-origin") == "http://localhost:5173"
    assert "access-control-allow-credentials" not in ok.headers
    assert api.post("/api/ncbi/import", json={"variation_ids": ["18390"]}).status_code == 403
    assert api.get("/api/variants/1/ncbi", params={"sync": "true"}).status_code == 405
    # DNS rebinding: unknown Host header is rejected
    assert api.get("/api/health", headers={"Host": "attacker.example"}).status_code == 400


# -- no-deadline-when-ncbi-slow

def test_circuit_breaker_budget_and_stale_cache(temp_db, monkeypatch):
    c = NcbiClient()
    c.limiter = RateLimiter(1000)
    c.backoff_base = 0.0
    calls = []

    def ok(url, params, timeout):
        calls.append(1)
        return (200, json.dumps({"esearchresult": {"count": "1", "idlist": ["9"]}}), {})

    def down(url, params, timeout):
        calls.append(1)
        raise httpx.ConnectTimeout("blackholed")

    monkeypatch.setattr(c, "_transport", ok)
    c.esearch("clinvar", "cached-term")
    c.cache.ttl_seconds = -1  # entry is now expired
    monkeypatch.setattr(c, "_transport", down)
    calls.clear()
    # network down: the expired entry is served instead of an error
    assert c.esearch("clinvar", "cached-term")["ids"] == ["9"]
    with pytest.raises(NcbiError):
        c.esearch("clinvar", "uncached-1")
    assert c.breaker_open()  # 2 consecutive failed requests
    n = len(calls)
    t0 = time.monotonic()
    with pytest.raises(NcbiError):
        c.esearch("clinvar", "uncached-2")
    assert len(calls) == n and time.monotonic() - t0 < 0.5  # fails fast without network
    # the status probe bypasses the breaker and closes it once NCBI answers again
    monkeypatch.setattr(c, "_transport", lambda u, p, t: (200, json.dumps({"esearchresult": {"count": "1", "idlist": []}}), {}))
    assert c.ping() == (True, None) and not c.breaker_open()
    # budgets: attempts cap and an exhausted deadline (which is not a connection failure)
    monkeypatch.setattr(c, "_transport", down)
    calls.clear()
    with c.budget(30, attempts=1):
        with pytest.raises(NcbiError):
            c.esearch("clinvar", "uncached-3")
    assert len(calls) == 1
    calls.clear()
    with c.budget(0):
        with pytest.raises(NcbiError):
            c.esearch("clinvar", "uncached-4")
    assert calls == []


# -- rate-limiter-send-time-429s

def test_rate_limiter_margin_and_ratelimit_header(temp_db, monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    c = NcbiClient()
    assert c.limiter.interval >= 0.4  # 3 req/s without a key, counted on arrival
    c.limiter = RateLimiter(1000)
    monkeypatch.setattr(c, "_transport", lambda u, p, t: (200, json.dumps({"esearchresult": {"count": "0", "idlist": []}}), {"X-RateLimit-Remaining": "0"}))
    c.esearch("clinvar", "quota")
    assert c.limiter._next_slot - time.monotonic() > 0.5  # every thread is held back ~1 s


# -- sync-uses-stale-cache

def test_sync_and_import_bypass_response_cache(api, fake_ncbi, temp_db):
    row = _db_row("variant_id = 'rs671'")
    api.get(f"/api/variants/{row['id']}/ncbi")  # warms the cache
    n = len(fake_ncbi.calls)
    api.get(f"/api/variants/{row['id']}/ncbi")
    assert len(fake_ncbi.calls) == n  # read path is served from cache
    api.post(f"/api/variants/{row['id']}/ncbi/sync", headers=APP_HEADERS)
    fresh = [u for u, _p in fake_ncbi.calls[n:]]
    assert any(u.endswith("esearch.fcgi") for u in fresh) and any(u.endswith("esummary.fcgi") for u in fresh)
    n2 = len(fake_ncbi.calls)
    api.post("/api/ncbi/import", headers=APP_HEADERS, json={"variation_ids": ["18390"]})
    assert any(u.endswith("esummary.fcgi") for u, _p in fake_ncbi.calls[n2:])


def test_negative_esearch_results_expire_quickly(temp_db, monkeypatch):
    c = NcbiClient()
    c.limiter = RateLimiter(1000)
    calls = []
    monkeypatch.setattr(c, "_transport", lambda u, p, t: calls.append(1) or (200, json.dumps({"esearchresult": {"count": "0", "idlist": []}}), {}))
    c.esearch("clinvar", "nothing-yet")
    c.esearch("clinvar", "nothing-yet")
    assert len(calls) == 1
    monkeypatch.setattr(ncbi_client, "NEGATIVE_CACHE_TTL_SECONDS", -1)
    c.esearch("clinvar", "nothing-yet")
    assert len(calls) == 2


# -- risk-freq-wrong-indel-allele

def test_risk_frequency_indel_names_and_explicit_zero():
    multi_del = {
        "spdi": "NC_000001.11:100:AAAAAAAAA:AAAAAAAA,NC_000001.11:100:AAAAAAAAA:",
        "alts": ["AAAAAAAA", "-"],
        "frequencies": [{"study": "TOPMED", "allele": "-", "frequency": 4e-06, "count": 1}],
    }
    # '-' could be either deletion: do not report another allele's frequency
    assert risk_allele_frequency_from_dbsnp(multi_del, ("NC_000001.11", 100, "AAAAAAAAA", "AAAAAAAA"), "A") is None
    single_del = dict(multi_del, spdi="NC_000001.11:100:AAAAAAAAA:AAAAAAAA", alts=["AAAAAAAA"])
    assert risk_allele_frequency_from_dbsnp(single_del, ("NC_000001.11", 100, "AAAAAAAAA", "AAAAAAAA"), "A")["frequency"] == 4e-06
    # the preferred study's explicit zero for the exact allele is a real 0
    rs671 = parse_dbsnp_summary(snp("671"))
    zero = dict(rs671, frequencies=[{"study": "ALFA", "allele": "A", "frequency": 0.0, "count": 0}, {"study": "TOPMED", "allele": "A", "frequency": 0.02, "count": 5}])
    assert risk_allele_frequency_from_dbsnp(zero, ("NC_000012.12", 111803961, "G", "A"), "A") == {"allele": "A", "frequency": 0.0, "study": "ALFA"}
