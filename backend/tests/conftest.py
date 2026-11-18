"""
Shared fixtures. Tests never touch real user data: the jobs root and the ClinVar cache are
redirected to throw-away temp folders before the application is imported, and every test
analysis is built from the synthetic genome in tests/synthetic.py.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import uuid
from pathlib import Path

_TEST_ROOT = Path(tempfile.mkdtemp(prefix="genomescope-tests-"))
os.environ["GENOMESCOPE_JOBS_DIR"] = str(_TEST_ROOT / "jobs")
os.environ["GENOMESCOPE_CLINVAR_DIR"] = str(_TEST_ROOT / "clinvar-cache")
os.environ.pop("NCBI_API_KEY", None)
os.environ.pop("NCBI_EMAIL", None)

import pytest  # noqa: E402
from fastapi.testclient import TestClient  # noqa: E402

from app.main import create_app  # noqa: E402
from app.services import job_store  # noqa: E402
from app.utils.temp_manager import jobs_root, purge_job_dir  # noqa: E402
from genomescope import workspace  # noqa: E402
from genomescope.workspace import Workspace  # noqa: E402
from tests import synthetic  # noqa: E402

JOB_HEADER = "X-GenomeScope-Job"
APP_HEADERS = {"X-Requested-With": "GenomeScope"}


def pytest_sessionfinish(session, exitstatus):
    shutil.rmtree(_TEST_ROOT, ignore_errors=True)


@pytest.fixture(scope="session")
def template_workspace(tmp_path_factory) -> Path:
    """A fully built synthetic analysis (inputs + catalog + sparse index), built once."""
    root = tmp_path_factory.mktemp("template") / str(uuid.uuid4())
    ws = Workspace(root, "Synthetic sample")
    synthetic.write_inputs(ws)
    with workspace.activate(ws):
        synthetic.build_catalog(ws)
    return root


@pytest.fixture
def job_ws(template_workspace) -> Workspace:
    """A fresh copy of the synthetic analysis, registered as a completed job and active."""
    job_id = str(uuid.uuid4())
    root = jobs_root() / job_id
    shutil.copytree(template_workspace, root)
    (root / "work" / "vcf_sparse_index.json").unlink(missing_ok=True)  # index paths are per copy
    ws = Workspace(root, "Synthetic sample")
    job_store.create(job_id=job_id, label=ws.sample_label, temp_dir=str(root), files=[], file_size=0,
                     options={"clinvar_mode": "off"}, stages=[])
    job_store.update(job_id, status=job_store.COMPLETED, progress=100)
    with workspace.activate(ws):
        yield ws
    purge_job_dir(root)
    job_store.delete(job_id)


@pytest.fixture
def app():
    return create_app(start_manager=False)


@pytest.fixture
def api(app, job_ws) -> TestClient:
    """Client bound to the synthetic analysis (sends the analysis header on every request)."""
    return TestClient(app, headers={JOB_HEADER: job_ws.root.name})


@pytest.fixture
def bare_api(app) -> TestClient:
    """Client without a selected analysis."""
    return TestClient(app)
