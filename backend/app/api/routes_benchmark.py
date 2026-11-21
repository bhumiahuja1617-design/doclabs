"""Launch and inspect reproducible two-target benchmark runs."""
import json
import subprocess
import sys
import threading
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, HTTPException, Response
from pydantic import BaseModel, Field

from backend.app.config import settings

router = APIRouter(prefix="/api/v1/benchmark", tags=["Benchmark"])
ROOT = settings.BASE_DIR
OUT = ROOT / "benchmark" / "results"
STATE = OUT / "run_status.json"
_lock = threading.Lock()


class BenchmarkRequest(BaseModel):
    target: Literal["all", "1hsg", "egfr"] = "all"
    exhaustiveness: int = Field(default=4, ge=1, le=32)
    seed: int = 42


def _execute(target: str, exhaustiveness: int, seed: int) -> None:
    try:
        OUT.mkdir(parents=True, exist_ok=True)
        STATE.write_text(json.dumps({"status": "running", "target": target, "exhaustiveness": exhaustiveness, "seed": seed}), encoding="utf-8")
        cmd = [sys.executable, str(ROOT / "benchmark" / "run_benchmark.py"), "--target", target,
               "--exhaustiveness", str(exhaustiveness), "--seed", str(seed)]
        completed = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
        STATE.write_text(json.dumps({"status": "completed" if completed.returncode == 0 else "failed",
                                     "target": target, "exhaustiveness": exhaustiveness, "seed": seed,
                                     "return_code": completed.returncode,
                                     "log": (completed.stdout + "\n" + completed.stderr)[-12000:]}), encoding="utf-8")
    except Exception as exc:
        STATE.write_text(json.dumps({"status": "failed", "target": target, "error": str(exc)[:1000]}), encoding="utf-8")
    finally:
        _lock.release()


@router.post("")
async def start_benchmark(payload: BenchmarkRequest):
    if not _lock.acquire(blocking=False):
        raise HTTPException(409, detail={"error": {"code": "BENCHMARK_RUNNING", "message": "A benchmark is already running.", "details": {}}})
    OUT.mkdir(parents=True, exist_ok=True)
    threading.Thread(target=_execute, args=(payload.target, payload.exhaustiveness, payload.seed), daemon=True, name="docklab-benchmark").start()
    return {"status": "running", "target": payload.target, "exhaustiveness": payload.exhaustiveness, "seed": payload.seed}


@router.get("")
async def benchmark_status():
    if STATE.exists():
        try:
            status = json.loads(STATE.read_text(encoding="utf-8"))
        except Exception:
            status = {"status": "unknown"}
    else:
        status = {"status": "idle"}
    result = OUT / "benchmark.json"
    if result.exists() and status.get("status") != "running":
        status["results"] = json.loads(result.read_text(encoding="utf-8"))
    report = ROOT / "benchmark" / "RESULTS.md"
    if report.exists() and status.get("status") != "running":
        status["report"] = report.read_text(encoding="utf-8")
    return status


@router.get("/report")
async def benchmark_report():
    report = ROOT / "benchmark" / "RESULTS.md"
    if not report.is_file():
        raise HTTPException(404, detail={"error": {"code": "BENCHMARK_NOT_RUN", "message": "Run the benchmark before opening its generated report.", "details": {}}})
    return Response(report.read_text(encoding="utf-8"), media_type="text/markdown; charset=utf-8")
