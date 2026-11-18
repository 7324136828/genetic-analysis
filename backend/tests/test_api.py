"""Catalog / VCF explorer API against the synthetic analysis (tests/synthetic.py)."""

import pytest

pytestmark = pytest.mark.usefixtures("job_ws")


def test_health_endpoint(api, bare_api):
    data = api.get("/api/health").json()
    assert data["status"] == "ok" and data["analysis_selected"] is True
    assert data["sample_id"] == "Synthetic sample"
    assert data["reference_genome"] == "GRCh38"
    assert data["indexed_variants"] >= 4
    assert len(data["vcf_files"]) == 3 and all(f["exists"] for f in data["vcf_files"])
    # no analysis selected: still healthy (used by the launcher), without genome data
    bare = bare_api.get("/api/health").json()
    assert bare["status"] == "ok" and bare["analysis_selected"] is False and "vcf_files" not in bare


def test_explorer_requires_a_completed_analysis(bare_api, api):
    assert bare_api.get("/api/summary").status_code == 409
    assert bare_api.get("/api/summary", headers={"X-GenomeScope-Job": "not-a-uuid"}).status_code == 400
    assert bare_api.get("/api/summary", headers={"X-GenomeScope-Job": "00000000-0000-4000-8000-000000000000"}).status_code == 404
    assert api.get("/api/summary").status_code == 200


def test_summary_endpoint(api):
    data = api.get("/api/summary").json()
    assert data["sample_id"] == "Synthetic sample"
    assert data["total_indexed_loci"] >= 4
    assert data["carried_variants"] >= 3
    assert data["structural_cnv_carried"] >= 1
    assert len(data["by_classification"]) > 0 and len(data["by_status"]) > 0 and len(data["by_variant_type"]) > 0
    assert any(g["gene"] == "ALDH2" for g in data["top_carried_genes"])


def test_search_by_all_eleven_fields(api):
    r1 = api.get("/api/variants", params={"rcv_id": "RCV000020060"}).json()
    assert r1["total"] >= 1
    item = r1["items"][0]
    for key in ["rcv_id", "chr", "position", "variant_id", "gene", "risk_frequency_text", "variant_type",
                "your_data", "risk_version", "your_status", "condition", "classification", "evidence"]:
        assert key in item

    r2 = api.get("/api/variants", params={"chr": "12", "position": "111803962"}).json()
    assert r2["total"] >= 1 and r2["items"][0]["gene"] == "ALDH2"

    r3 = api.get("/api/variants", params={"variant_id": "rs671"}).json()
    assert r3["total"] >= 1 and r3["items"][0]["variant_id"] == "rs671"

    assert api.get("/api/variants", params={"gene": "MTHFR"}).json()["total"] >= 1
    assert api.get("/api/variants", params={"risk_frequency": "<0.01", "carried_only": True}).json()["total"] >= 1
    assert api.get("/api/variants", params={"variant_type": "CNV Loss"}).json()["total"] >= 1
    assert api.get("/api/variants", params={"your_status": "Carrier (Heterozygous)"}).json()["total"] >= 1
    assert api.get("/api/variants", params={"condition": "Alcohol flushing"}).json()["total"] >= 1
    assert api.get("/api/variants", params={"classification": "Drug Response", "carried_only": True}).json()["total"] >= 1
    assert api.get("/api/variants", params={"evidence": "4 Stars"}).json()["total"] >= 1


def test_variant_detail_and_direct_vcf_read(api):
    rec_id = api.get("/api/variants", params={"variant_id": "rs671"}).json()["items"][0]["id"]
    detail = api.get(f"/api/variants/{rec_id}").json()
    assert detail["variant_id"] == "rs671" and detail["your_status"] == "Detected (Homozygous)"
    assert len(detail["context_vcf_reads"]) >= 1

    r_vcf = api.get("/api/vcf/read", params={"chr": "12", "position": "111803900-111804050"}).json()
    assert r_vcf["count"] >= 1
    assert any(r["position"] == 111803962 for r in r_vcf["reads"])
    assert api.get("/api/vcf/read", params={"chr": "12", "position": "abc"}).status_code == 400


def test_source_files_use_canonical_names(api):
    files = api.get("/api/filters").json()["source_files"]
    assert set(files) <= {"snp-indel.genome.vcf", "sv.vcf", "cnv.vcf"}


def test_api_responses_are_not_cached(api):
    r = api.get("/api/summary")
    assert r.headers.get("cache-control") == "no-store"
