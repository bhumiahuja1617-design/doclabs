from contextlib import asynccontextmanager
from pathlib import Path
from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException
from backend.app.config import settings
from backend.app.api.routes_jobs import router as jobs_router
from backend.app.api.routes_examples import router as examples_router
from backend.app.api.routes_benchmark import router as benchmark_router
from backend.app.api.routes_proteins import router as proteins_router
from backend.app.api.routes_ligands import router as ligands_router
from backend.app.db.session import init_db

@asynccontextmanager
async def lifespan(app: FastAPI):
    init_db()
    yield

app = FastAPI(
    title=settings.PROJECT_NAME,
    description="DockLab: AutoDock Vina Docking, Structure QC, and Multi-Seed Reliability Engine",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/openapi.json",
    lifespan=lifespan,
)

app.include_router(jobs_router)
app.include_router(examples_router)
app.include_router(benchmark_router)
app.include_router(proteins_router)
app.include_router(ligands_router)

# CORS Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """Ensure no unhandled internal exception leaks raw stack trace to user."""
    return JSONResponse(
        status_code=500,
        content={
            "error": {
                "code": "INTERNAL_SERVER_ERROR",
                "message": "An unexpected error occurred. Please contact the administrator.",
                "details": {"exception_type": exc.__class__.__name__} if settings.DEBUG else {},
            }
        },
    )

@app.exception_handler(RequestValidationError)
async def request_validation_handler(request: Request, exc: RequestValidationError):
    """Return validation failures in the same structured format as API errors."""
    return JSONResponse(
        status_code=422,
        content={"error": {
            "code": "INVALID_REQUEST",
            "message": "The request does not match the expected input.",
            "details": [{"field": ".".join(str(part) for part in issue["loc"]),
                         "message": issue["msg"], "type": issue["type"]}
                        for issue in exc.errors()],
        }},
    )

@app.get("/health", tags=["System"])
async def health_check():
    """Liveness probe."""
    return {"status": "ok", "app": settings.PROJECT_NAME, "version": "1.0.0"}

@app.get("/api/v1/health", tags=["System"])
async def api_health_check():
    """API prefix liveness probe."""
    return {"status": "ok", "app": settings.PROJECT_NAME, "version": "1.0.0"}

# The standalone Render preview serves the built SPA from this same origin.
# In Docker Compose, the separate Nginx frontend handles static files and proxies /api.
frontend_dist = Path(settings.BASE_DIR) / "frontend" / "dist"
if (frontend_dist / "index.html").is_file():
    class SpaStaticFiles(StaticFiles):
        async def get_response(self, path: str, scope):
            try:
                return await super().get_response(path, scope)
            except StarletteHTTPException as exc:
                reserved = ("api/", "docs", "redoc", "openapi.json", "health")
                if exc.status_code != 404 or path.startswith(reserved):
                    raise
                return await super().get_response("index.html", scope)

    app.mount("/", SpaStaticFiles(directory=frontend_dist, html=True), name="frontend")

