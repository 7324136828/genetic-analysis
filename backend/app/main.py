"""
GenomeScope API entrypoint.

    uvicorn app.main:app --host 127.0.0.1 --port 8000      (run from backend/)
"""

from __future__ import annotations

from contextlib import asynccontextmanager
from typing import Any, Dict, Optional

from fastapi import Depends, FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.trustedhost import TrustedHostMiddleware
from fastapi.responses import JSONResponse

from app.api import explorer, jobs
from app.api.deps import JOB_HEADER, optional_analysis
from app.config import APP_VERSION, get_settings
from app.services.pipeline import manager
from app.utils.temp_manager import jobs_root
from genomescope import indexer, ncbi_service
from genomescope.vcf_engine import engine_for
from genomescope.workspace import NoActiveWorkspace, Workspace


def create_app(start_manager: bool = True) -> FastAPI:
    settings = get_settings()

    @asynccontextmanager
    async def lifespan(_app: FastAPI):
        jobs_root()
        if start_manager:
            manager.start()
        try:
            yield
        finally:
            if start_manager:
                manager.shutdown()

    app = FastAPI(
        title="GenomeScope API",
        description=(
            "Whole-genome VCF analysis: upload a gVCF (plus optional SV / CNV VCFs), build an annotated "
            "variant catalog in an isolated temp folder, verify it against a local ClinVar release, "
            "explore it, and download the results."
        ),
        version=APP_VERSION,
        lifespan=lifespan,
    )

    app.add_middleware(TrustedHostMiddleware, allowed_hosts=settings.allowed_hosts)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cors_origins,
        allow_credentials=False,
        allow_methods=["GET", "POST", "DELETE"],
        allow_headers=["Content-Type", "X-Requested-With", JOB_HEADER],
    )

    @app.middleware("http")
    async def no_store(request: Request, call_next):
        response = await call_next(request)
        # Genome data must not linger in shared / browser caches.
        if request.url.path.startswith("/api/"):
            response.headers.setdefault("Cache-Control", "no-store")
            response.headers.setdefault("X-Content-Type-Options", "nosniff")
        return response

    @app.exception_handler(NoActiveWorkspace)
    async def _no_workspace(_request: Request, exc: NoActiveWorkspace):
        return JSONResponse(status_code=409, content={"detail": str(exc)})

    @app.get("/api/health", tags=["meta"])
    def health(ws: Optional[Workspace] = Depends(optional_analysis)) -> Dict[str, Any]:
        body: Dict[str, Any] = {
            "status": "ok",
            "version": APP_VERSION,
            "analysis_selected": ws is not None,
            "ncbi": ncbi_service.health_info(),
        }
        if ws is not None:
            conn = indexer.get_db_connection()
            try:
                body["indexed_variants"] = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
            finally:
                conn.close()
            body["sample_id"] = ws.sample_label
            body["reference_genome"] = ws.reference
            body["vcf_files"] = engine_for(ws).get_vcf_files_metadata()
        return body

    app.include_router(jobs.router)
    app.include_router(jobs.meta_router)
    app.include_router(explorer.status_router)
    app.include_router(explorer.router)
    return app


app = create_app()
