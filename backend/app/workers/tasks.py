"""Per-molecule Celery tasks, job aggregation, and seven-day cleanup."""
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from celery import shared_task
from sqlalchemy import update

from backend.app.config import settings
from backend.app.db.models import Job
from backend.app.db.session import SessionLocal
from backend.app.services.analysis_service import AnalysisService
from backend.app.services.docking_service import DockingService
from backend.app.services.ligand_service import LigandService
from backend.app.services.reliability_service import ReliabilityService


@shared_task(name="docklab.dock_molecule", autoretry_for=(), time_limit=660, soft_time_limit=630)
def dock_molecule(job_id: str, index: int) -> dict:
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if job is None or job.status == "cancelled":
            return {"status": "cancelled", "index": index}
        if job.status == "queued":
            job.status = "running"
            session.commit()
        ligand = job.ligands[index]
        settings_snapshot = dict(job.settings)
        receptor = Path(job.protein_path)
        box = dict(job.box)
    root = Path(settings.DATA_DIR) / job_id
    result_dir = root / "results"
    result_file = result_dir / f"ligand-{index:03d}.json"
    try:
        prep = root / "prepared" / f"{index:03d}.pdbqt"
        pose = result_dir / f"ligand-{index:03d}.pdbqt"
        sdf = result_dir / f"ligand-{index:03d}.sdf"
        prep_result = LigandService().prepare_ligand_from_smiles(ligand["smiles"], prep, seed=int(settings_snapshot.get("seed", 42)))
        docked = DockingService().dock(receptor, prep, box["center"], box["size"], pose, sdf,
            exhaustiveness=int(settings_snapshot.get("exhaustiveness", 8)), seed=int(settings_snapshot.get("seed", 42)),
            timeout_seconds=settings.DEFAULT_VINA_TIMEOUT_SECONDS)
        props = prep_result["properties"]
        heavy_atoms = max(1, int(props["heavy_atom_count"]))
        result = {"name": ligand["name"], "label": ligand["label"], "score": docked["best_score"],
                  "ligand_efficiency": round(float(docked["best_score"]) / heavy_atoms, 3), "properties": props,
                  "prepared_pdbqt": str(prep), "pose": str(pose), "sdf": str(sdf), "status": "docked", "index": index}
    except Exception as exc:
        result = {"name": ligand["name"], "label": ligand["label"], "status": "failed", "error": str(exc)[:500], "index": index}
    result_dir.mkdir(parents=True, exist_ok=True)
    result_file.write_text(json.dumps(result), encoding="utf-8")
    with SessionLocal() as session:
        session.execute(update(Job).where(Job.id == job_id).values(progress_done=Job.progress_done + 1))
        session.commit()
    return result


@shared_task(name="docklab.aggregate_job")
def aggregate_job(results: list[dict], job_id: str) -> dict:
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if job is None:
            return {"error": "job not found"}
        if job.status == "cancelled":
            return {"job_id": job_id, "status": "cancelled"}
        all_rows = sorted(results, key=lambda item: item.get("index", 0))
        docked = [row for row in all_rows if row.get("status") == "docked"]
        docked.sort(key=lambda item: item["score"])
        for rank, item in enumerate(docked, 1):
            item["rank"] = rank
        failed = [row for row in all_rows if row.get("status") != "docked"]
        metrics = None
        labeled = [row for row in docked if row.get("label") in {"active", "decoy"}]
        if sum(row["label"] == "active" for row in labeled) >= 3 and sum(row["label"] == "decoy" for row in labeled) >= 3:
            metrics = AnalysisService().screening_metrics([int(row["label"] == "active") for row in labeled], [row["score"] for row in labeled])
        reliability_summary = None
        if job.settings.get("reliability") and docked:
            job.status = "reliability"
            session.commit()
            top_k = min(int(job.settings.get("top_k", 10)), 20, len(docked))
            seed_count = min(int(job.settings.get("n_seeds", 5)), 10)
            first_seed = int(job.settings.get("seed", 42))
            seeds = list(range(first_seed, first_seed + seed_count))
            service = ReliabilityService(DockingService())
            seed_scores = {}
            for row in docked[:top_k]:
                try:
                    row["reliability"] = service.run_ligand(
                        Path(job.protein_path), Path(row["prepared_pdbqt"]), job.box["center"], job.box["size"],
                        Path(settings.DATA_DIR) / job_id / "reliability" / row["name"], seeds,
                        exhaustiveness=int(job.settings.get("exhaustiveness", 8)), initial_score=float(row["score"]),
                        initial_pose_sdf=Path(row["sdf"]),
                    )
                    seed_scores[row["name"]] = row["reliability"]["scores"]
                except Exception as exc:
                    row["reliability"] = {"label": "unavailable", "error": str(exc)[:300]}
            if seed_scores:
                rank_summary = ReliabilityService.rank_stability(seed_scores)
                for row in docked[:top_k]:
                    if row.get("reliability", {}).get("label") in {"high", "moderate", "low"}:
                        row["reliability"]["rank_range"] = rank_summary["rank_ranges"].get(row["name"])
                        row["reliability"]["mean_spearman"] = rank_summary["mean_spearman"]
                reliability_summary = {"top_k": top_k, "seed_count": seed_count, "mean_spearman": rank_summary["mean_spearman"]}
        for row in docked + failed:
            row.pop("index", None)
        job.results = docked + failed
        job.settings = {**job.settings, "analysis": metrics, "reliability_summary": reliability_summary}
        job.status = "completed" if docked else "failed"
        session.commit()
        return {"job_id": job_id, "status": job.status, "completed": len(docked), "failed": len(failed)}


@shared_task(name="docklab.cleanup_expired")
def cleanup_expired() -> int:
    """Delete expired job folders and rows. Called daily by Celery Beat."""
    now = datetime.now(timezone.utc)
    with SessionLocal() as session:
        expired = session.query(Job).filter(Job.expires_at < now).all()
        ids = [item.id for item in expired]
        for item in expired:
            folder = Path(settings.DATA_DIR) / item.id
            if folder.exists() and folder.resolve().parent == Path(settings.DATA_DIR).resolve():
                import shutil
                shutil.rmtree(folder)
            session.delete(item)
        session.commit()
        return len(ids)
