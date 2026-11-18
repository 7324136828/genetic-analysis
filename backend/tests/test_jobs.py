"""
Analysis jobs: upload validation and staging, the worker pipeline, history, ZIP download,
discard / delete and temp-folder hygiene. Everything runs offline (GENOMESCOPE_OFFLINE=1
or stubbed annotation sources) on the synthetic genome.
"""

from __future__ import annotations

import gzip
import io
import json
import time
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app import worker
from app.main import create_app
from app.services import job_store, staging
from app.services.pipeline import manager
from app.utils import temp_manager
from genomescope import indexer
from tests import synthetic
from tests.conftest import APP_HEADERS

GVCF = ("\n".join(synthetic.gvcf_lines()) + "\n").encode()
SV = ("\n".join(synthetic.sv_lines()) + "\n").encode()
CNV = ("\n".join(synthetic.cnv_lines()) + "\n").encode()

MYVARIANT_HIT = {
    "clinvar": {
        "chrom": "12", "hg38": {"start": 111803962, "end": 111803962}, "rsid": "rs671", "gene": {"symbol": "ALDH2"},
        "ref": "G", "alt": "A",
        "rcv": [{"accession": "RCV000020060", "clinical_significance": "drug response",
                 "review_status": "reviewed by expert panel", "conditions": {"name": "Alcohol flushing response (synthetic)"}}],
    },
    "gnomad_genome": {"af": {"af": 0.0061}},
}


def _upload(client, files, data=None, headers=APP_HEADERS):
    return client.post("/api/jobs", files=files, data=data or {"label": "Test run", "clinvar_mode": "off"}, headers=headers)


@pytest.fixture
def queued_only(monkeypatch):
    """Jobs are staged but never started (no worker subprocess)."""
    submitted = []
    monkeypatch.setattr(manager, "submit", lambda job_id: submitted.append(job_id))
    return submitted


def _cleanup(job_id):
    job = job_store.get(job_id)
    if job:
        temp_manager.purge_job_dir(job["temp_dir"])
        job_store.delete(job_id)


# ---------------------------------------------------------------------------
# Upload validation & staging
# ---------------------------------------------------------------------------

def test_upload_requires_app_header_and_gvcf(bare_api, queued_only):
    assert _upload(bare_api, {"snp_indel": ("a.vcf", GVCF)}, headers={}).status_code == 403
    assert _upload(bare_api, {"sv": ("sv.vcf", SV)}).status_code in (400, 422)
    assert queued_only == []


@pytest.mark.parametrize(
    "content, message",
    [
        (b"hello world\n", "not a VCF"),
        (b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", "no sample"),
        (b"##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n", "GRCh37"),
        (b"##fileformat=VCFv4.2\n##reference=file:///ref/hs37d5.fa\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n", "GRCh37"),
        (b"", "empty"),
    ],
)
def test_invalid_uploads_are_rejected_without_leftovers(bare_api, queued_only, content, message):
    before = {p.name for p in temp_manager.jobs_root().iterdir() if p.is_dir()}
    r = _upload(bare_api, {"snp_indel": ("bad.vcf", content)})
    assert r.status_code == 400 and message.lower() in r.json()["detail"].lower()
    after = {p.name for p in temp_manager.jobs_root().iterdir() if p.is_dir()}
    assert after == before  # the half-created job folder was purged
    assert queued_only == []


def test_gzip_upload_is_staged_under_canonical_names(bare_api, queued_only):
    gz = gzip.compress(GVCF)
    r = _upload(bare_api, {
        "snp_indel": ("KIT123-sample-30x.snp-indel.genome.vcf.gz", gz),
        "sv": ("KIT123-sample.sv.vcf", SV),
        "cnv": ("KIT123-sample.cnv.vcf", CNV),
    })
    assert r.status_code == 202, r.text
    job = r.json()
    try:
        assert job["status"] == "queued" and job["label"] == "Test run" and queued_only == [job["id"]]
        assert {f["slot"] for f in job["files"]} == {"snp_indel", "sv", "cnv"}
        snp = next(f for f in job["files"] if f["slot"] == "snp_indel")
        assert snp["compressed_upload"] is True and snp["size_bytes"] == len(GVCF) and snp["reference"] == "GRCh38"
        assert "temp_dir" not in job and "zip_path" not in job  # internal paths are never exposed
        root = Path(job_store.get(job["id"])["temp_dir"])
        assert root.parent == temp_manager.jobs_root()
        names = sorted(p.name for p in (root / "inputs").iterdir())
        assert names == ["cnv.vcf", "manifest.json", "snp-indel.genome.vcf", "sv.vcf"]
        assert (root / "inputs" / "snp-indel.genome.vcf").read_bytes() == GVCF
        manifest = json.loads((root / "inputs" / "manifest.json").read_text())
        assert manifest["reference"] == "GRCh38"
        # the user's file name is kept only in the manifest / history, never as a data file name
        assert not any("KIT123" in p.name for p in root.rglob("*"))
    finally:
        _cleanup(job["id"])


def test_label_defaults_and_is_sanitised(bare_api, queued_only):
    r = _upload(bare_api, {"snp_indel": ("a.vcf", GVCF)}, data={"label": "  \x00\x07  ", "clinvar_mode": "off"})
    job = r.json()
    try:
        assert job["label"] == f"Analysis {job['id'][:8]}"
    finally:
        _cleanup(job["id"])
    assert _upload(bare_api, {"snp_indel": ("a.vcf", GVCF)}, data={"clinvar_mode": "bogus"}).status_code == 422


# ---------------------------------------------------------------------------
# Worker pipeline (in-process, annotation sources stubbed)
# ---------------------------------------------------------------------------

def test_worker_pipeline_builds_catalog_and_zip(bare_api, queued_only, monkeypatch, capsys):
    monkeypatch.setattr(indexer, "query_myvariant_batch", lambda q, size=1000: [MYVARIANT_HIT])
    monkeypatch.setattr(indexer, "annotate_genomic_interval_genes", lambda c, s, e: ("SYNGENE1", "synthetic gene"))
    job = _upload(bare_api, {"snp_indel": ("g.vcf", GVCF), "sv": ("s.vcf", SV), "cnv": ("c.vcf", CNV)}).json()
    try:
        root = Path(job_store.get(job["id"])["temp_dir"])
        result = worker.run(root, "Test run", "off", {"id": job["id"], "label": "Test run"})
        out = capsys.readouterr().out
        stages = [json.loads(line[8:]) for line in out.splitlines() if line.startswith("@@STAGE ")]
        assert [(s["key"], s["status"]) for s in stages] == [
            ("index", "running"), ("index", "done"), ("catalog", "running"), ("catalog", "done"),
            ("clinvar", "running"), ("clinvar", "skipped"), ("package", "running"), ("package", "done"),
        ]
        assert [s["progress"] for s in stages] == sorted(s["progress"] for s in stages)
        assert result["catalog"]["total_records"] >= 3 and result["clinvar"]["status"] == "skipped"
        zpath = Path(result["package"]["zip_path"])
        assert zpath.parent == root / "archive"
        with zipfile.ZipFile(zpath) as zf:
            names = set(zf.namelist())
            assert {"variants_catalog.csv", "carried_variants.csv", "analysis_summary.json"} <= names
            summary = json.loads(zf.read("analysis_summary.json"))
            assert summary["summary"]["sample_id"] == "Test run"
            csv_text = zf.read("variants_catalog.csv").decode()
            assert "rs671" in csv_text and "ALDH2" in csv_text
    finally:
        _cleanup(job["id"])


# ---------------------------------------------------------------------------
# Full lifecycle through the real job manager + worker subprocess (offline)
# ---------------------------------------------------------------------------

def _wait(client, job_id, states, timeout=120):
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["status"] in states:
            return job
        time.sleep(0.3)
    raise AssertionError(f"job did not reach {states}: {job}")


def test_full_lifecycle_with_worker_subprocess(monkeypatch):
    monkeypatch.setenv("GENOMESCOPE_OFFLINE", "1")
    with TestClient(create_app(start_manager=True)) as client:
        r = _upload(client, {"snp_indel": ("g.vcf.gz", gzip.compress(GVCF)), "cnv": ("c.vcf", CNV)})
        assert r.status_code == 202, r.text
        job_id = r.json()["id"]
        job = _wait(client, job_id, ("completed", "failed"))
        assert job["status"] == "completed", job.get("error_message")
        assert job["progress"] == 100 and job["can_open"] and job["can_download"]
        assert [s["status"] for s in job["stages"]] == ["done", "done", "done", "skipped", "done"]
        assert job["elapsed_seconds"] is not None

        # history lists it; logs are available
        assert any(j["id"] == job_id for j in client.get("/api/jobs").json())
        logs = client.get(f"/api/jobs/{job_id}/logs").json()["lines"]
        assert any("Offline mode" in line for line in logs)

        # explorer works on the finished analysis
        hdr = {"X-GenomeScope-Job": job_id}
        summary = client.get("/api/summary", headers=hdr).json()
        assert summary["total_indexed_loci"] >= 1  # CNV rows (offline: no MyVariant rows)
        assert client.get("/api/vcf/read", params={"chr": "12", "position": "111803962"}, headers=hdr).json()["count"] >= 1

        # download
        z = client.get(f"/api/jobs/{job_id}/download-zip")
        assert z.status_code == 200 and z.headers["content-type"] == "application/zip"
        assert "Test_run_genomescope_results.zip" in z.headers["content-disposition"]
        with zipfile.ZipFile(io.BytesIO(z.content)) as zf:
            assert "pipeline.log" in zf.namelist() and "variants_catalog.csv" in zf.namelist()

        # a completed analysis cannot be "discarded" (only deleted)
        assert client.post(f"/api/jobs/{job_id}/discard", headers=APP_HEADERS).status_code == 409

        # delete removes the folder and the ledger row
        root = Path(job_store.get(job_id)["temp_dir"])
        assert client.delete(f"/api/jobs/{job_id}").status_code == 403  # CSRF header required
        assert client.delete(f"/api/jobs/{job_id}", headers=APP_HEADERS).status_code == 204
        assert not root.exists() and client.get(f"/api/jobs/{job_id}").status_code == 404


def test_discard_queued_job_purges_temp_folder(bare_api, queued_only):
    job = _upload(bare_api, {"snp_indel": ("g.vcf", GVCF)}).json()
    root = Path(job_store.get(job["id"])["temp_dir"])
    assert root.exists()
    r = bare_api.post(f"/api/jobs/{job['id']}/discard", headers=APP_HEADERS)
    assert r.status_code == 200 and r.json()["status"] == "discarded"
    assert not root.exists()
    # the explorer refuses a discarded analysis; download is refused too
    assert bare_api.get("/api/summary", headers={"X-GenomeScope-Job": job["id"]}).status_code == 409
    assert bare_api.get(f"/api/jobs/{job['id']}/download-zip").status_code == 400
    _cleanup(job["id"])


def test_download_rebuilds_zip_after_catalog_changes(api, job_ws):
    job_id = job_ws.root.name
    first = api.get(f"/api/jobs/{job_id}/download-zip")
    assert first.status_code == 200
    zip_path = job_ws.archive_dir / "genomescope_results.zip"
    stamp = zip_path.stat().st_mtime
    time.sleep(1.1)
    conn = indexer.get_db_connection()
    conn.execute("UPDATE variants SET condition = 'edited' WHERE variant_id = 'rs671'")
    conn.commit()
    conn.close()
    api.get(f"/api/jobs/{job_id}/download-zip")
    assert zip_path.stat().st_mtime > stamp


# ---------------------------------------------------------------------------
# Temp-folder hygiene
# ---------------------------------------------------------------------------

def test_purge_refuses_paths_outside_jobs_root(tmp_path):
    outside = tmp_path / "keep-me"
    outside.mkdir()
    with pytest.raises(ValueError):
        temp_manager.purge_job_dir(outside)
    assert outside.exists()
    with pytest.raises(ValueError):
        temp_manager.job_dir("../../etc")


def test_sweep_removes_expired_jobs_and_old_orphans(monkeypatch):
    jid = temp_manager.new_job_id()
    ws = temp_manager.create_job_workspace(jid, "old")
    job_store.create(jid, "old", str(ws.root), [], 0, {}, [])
    job_store.update(jid, status=job_store.COMPLETED, completed_at="2000-01-01T00:00:00Z")
    orphan = temp_manager.create_job_workspace(temp_manager.new_job_id(), "orphan").root
    fresh_orphan = temp_manager.create_job_workspace(temp_manager.new_job_id(), "uploading").root
    import os
    old = time.time() - 7 * 3600
    os.utime(orphan, (old, old))
    removed = manager.sweep_expired()
    assert removed >= 2
    assert not ws.root.exists() and job_store.get(jid) is None and not orphan.exists()
    assert fresh_orphan.exists()  # could be an upload in progress
    temp_manager.purge_job_dir(fresh_orphan)


def test_staging_helpers():
    assert staging.safe_display_name("C:\\Users\\x\\My Genome.vcf") == "My Genome.vcf"
    assert staging.safe_display_name("../../etc/passwd") == "passwd"
    assert staging.safe_display_name(None) == "upload.vcf"
