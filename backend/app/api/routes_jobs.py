"""REST endpoints for draft jobs, input preparation, QC, execution, and exports."""
import csv
import io
import json
import os
import re
import threading
import urllib.request
import uuid
import zipfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

from fastapi import APIRouter, HTTPException, Response, Request
from celery import chord, group

from backend.app.config import settings
from backend.app.db.models import Job
from backend.app.db.session import SessionLocal
from backend.app.db.schemas import BoxModel, LigandsRequest, SettingsRequest, RunRequest
from backend.app.services.protein_service import ProteinService
from backend.app.services.ligand_service import LigandService
from backend.app.services.box_service import BoxService
from backend.app.services.qc_service import QCService
from backend.app.services.report_service import ReportService
from backend.app.workers.tasks import dock_molecule, aggregate_job
from backend.app.api.routes_proteins import _json_request

router = APIRouter(prefix="/api/v1/jobs", tags=["Jobs"])
EXAMPLES = settings.BASE_DIR / "examples"


def _run_local_job(job_id: str) -> None:
    """Run jobs in a background thread for the no-Docker Windows development setup."""
    rows = []
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        count = len(job.ligands) if job else 0
    for index in range(count):
        rows.append(dock_molecule.run(job_id, index))
    aggregate_job.run(rows, job_id)


def _job_data(job: Job) -> dict:
    return {"id": job.id, "created_at": job.created_at.isoformat(), "status": job.status,
            "box": job.box, "settings": job.settings, "progress": {"done": job.progress_done, "total": job.progress_total},
            "error_message": job.error_message, "expires_at": job.expires_at.isoformat() if job.expires_at else None}


def _target_metadata(value: dict | None, pdb_id: str, source: str, structure_type: str) -> dict:
    """Keep a small, report-safe target record with each prepared receptor."""
    value = value if isinstance(value, dict) else {}
    allowed = ("name", "title", "organism", "uniprot_id", "uniprot_ids", "method",
               "resolution", "resolution_angstroms", "function_summary", "disease_relevance",
               "source_url", "source", "confidence", "chain_ids", "pdb_id")
    target = {key: value[key] for key in allowed if value.get(key) not in (None, "", [], {})}
    target.update({"pdb_id": pdb_id, "source": target.get("source") or source,
                   "structure_type": structure_type})
    target.setdefault("name", target.get("title") or pdb_id)
    return target


def _available_job_results(job: Job) -> list[dict]:
    """Read completed per-molecule artifacts while a background job is still running."""
    if job.results:
        return list(job.results)
    root = (Path(settings.DATA_DIR) / job.id / "results").resolve()
    if Path(settings.DATA_DIR).resolve() not in root.parents or not root.is_dir():
        return []
    rows = []
    for path in root.glob("ligand-*.json"):
        try:
            index = int(path.stem.rsplit("-", 1)[1])
            if index < 0 or index >= len(job.ligands or []):
                continue
            row = json.loads(path.read_text(encoding="utf-8"))
            row.pop("index", None)
            rows.append(row)
        except (OSError, ValueError, json.JSONDecodeError):
            continue
    docked = sorted((row for row in rows if row.get("status") == "docked"), key=lambda row: row.get("score", 0))
    for rank, row in enumerate(docked, 1):
        row["rank"] = rank
    failed = [row for row in rows if row.get("status") != "docked"]
    return docked + failed


def _available_screening_metrics(rows: list[dict]) -> dict | None:
    labeled = [row for row in rows if row.get("status") == "docked" and row.get("label") in {"active", "decoy"}]
    if sum(row["label"] == "active" for row in labeled) < 3 or sum(row["label"] == "decoy" for row in labeled) < 3:
        return None
    from backend.app.services.analysis_service import AnalysisService
    return AnalysisService().screening_metrics([int(row["label"] == "active") for row in labeled],
                                               [row["score"] for row in labeled])


@router.post("")
async def create_job():
    job = Job(id=str(uuid.uuid4()), expires_at=datetime.now(timezone.utc) + timedelta(days=settings.JOB_RETENTION_DAYS))
    with SessionLocal() as session:
        session.add(job)
        session.commit()
    return _job_data(job)


@router.post("/{job_id}/protein")
async def set_protein(job_id: str, payload: dict):
    pdb_id = str(payload.get("pdb_id", "")).upper().strip()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise HTTPException(422, detail={"error": {"code": "INVALID_PDB_ID", "message": "Provide a four-character PDB ID.", "details": {}}})
    with SessionLocal() as session:
        if not session.get(Job, job_id):
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
    root = Path(settings.DATA_DIR) / job_id / "protein"
    root.mkdir(parents=True, exist_ok=True)
    try:
        svc = ProteinService()
        raw = svc.fetch_or_load_pdb(pdb_id, root)
        inspection = svc.inspect_structure(raw)
        cleaned = root / f"{pdb_id}_cleaned.pdb"
        pdbqt = root / f"{pdb_id}_cleaned.pdbqt"
        requested_chains = payload.get("chain_ids") or ([payload["chain_id"]] if payload.get("chain_id") else None)
        if requested_chains:
            unknown = sorted(set(requested_chains) - set(inspection["chains"]))
            if unknown:
                raise ValueError(f"Requested chain(s) not found: {', '.join(unknown)}")
        prep = svc.clean_and_prepare_receptor(raw, cleaned, pdbqt, keep_chains=requested_chains)
    except Exception as exc:
        raise HTTPException(422, detail={"error": {"code": "PROTEIN_PREPARATION_FAILED", "message": str(exc)[:300], "details": {}}}) from exc
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        job.protein_path = str(pdbqt)
        target_payload = payload.get("target")
        if not isinstance(target_payload, dict):
            try:
                entry = _json_request(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
                info = entry.get("rcsb_entry_info", {})
                methods = [item.get("method") for item in entry.get("exptl", []) if item.get("method")]
                target_payload = {"name": entry.get("struct", {}).get("title", pdb_id),
                                  "title": entry.get("struct", {}).get("title", ""),
                                  "method": ", ".join(methods),
                                  "resolution": (info.get("resolution_combined") or [None])[0],
                                  "source_url": f"https://www.rcsb.org/structure/{pdb_id}"}
            except Exception:
                target_payload = {"name": pdb_id, "source_url": f"https://www.rcsb.org/structure/{pdb_id}"}
        target = _target_metadata(target_payload, pdb_id, "RCSB PDB", "experimental")
        target["chain_ids"] = requested_chains or inspection["chains"]
        job.settings = {**job.settings, "protein": {**target, "target": target, "cleaned_pdb": str(cleaned), "raw_pdb": str(raw), "chains": requested_chains or inspection["chains"], "available_chains": inspection["chains"], "ligands_detected": inspection["ligands_detected"], "warnings": prep["warnings"]}}
        job.status = "created"
        session.commit()
    return {**target, "chains": inspection["chains"], "ligands_detected": inspection["ligands_detected"], "warnings": prep["warnings"]}


@router.post("/{job_id}/protein/alphafold")
async def set_alphafold_protein(job_id: str, payload: dict):
    """Retrieve the latest AFDB coordinate URL from its metadata API and prepare the model."""
    accession = str(payload.get("accession", "")).strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{6,10}", accession):
        raise HTTPException(422, detail={"error": {"code": "INVALID_UNIPROT_ACCESSION", "message": "Provide a UniProt accession from the selected AlphaFold result.", "details": {}}})
    with SessionLocal() as session:
        if not session.get(Job, job_id):
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
    try:
        models = _json_request(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
        models = models if isinstance(models, list) else [models]
        model = next((x for x in models if x.get("pdbUrl") or x.get("cifUrl")), None)
        if not model:
            raise ValueError("AlphaFold DB has no coordinate model for this accession.")
        coordinate_url = model.get("pdbUrl") or model.get("cifUrl")
        request = urllib.request.Request(coordinate_url, headers={"User-Agent": "DockLab/1.0"})
        with urllib.request.urlopen(request, timeout=30) as response:
            content = response.read(100 * 1024 * 1024 + 1)
        if len(content) > 100 * 1024 * 1024:
            raise ValueError("Predicted structure exceeds the 100 MB model size limit.")
        root = Path(settings.DATA_DIR) / job_id / "protein"
        root.mkdir(parents=True, exist_ok=True)
        raw = root / f"{accession}_alphafold.pdb"
        if model.get("pdbUrl"):
            raw.write_bytes(content)
        else:
            from Bio.PDB import MMCIFParser, PDBIO
            cif_path = root / f"{accession}_alphafold.cif"
            cif_path.write_bytes(content)
            structure = MMCIFParser(QUIET=True).get_structure(accession, str(cif_path))
            writer = PDBIO(); writer.set_structure(structure); writer.save(str(raw))
        svc = ProteinService()
        inspection = svc.inspect_structure(raw)
        cleaned = root / f"{accession}_cleaned.pdb"
        pdbqt = root / f"{accession}_cleaned.pdbqt"
        prep = svc.clean_and_prepare_receptor(raw, cleaned, pdbqt)
    except Exception as exc:
        raise HTTPException(422, detail={"error": {"code": "ALPHAFOLD_PREPARATION_FAILED", "message": str(exc)[:300], "details": {}}}) from exc
    target_info = {"pdb_id": model.get("modelEntityId") or model.get("entryId") or f"AF-{accession}-F1",
                   "uniprot_id": accession, "source": "AlphaFold DB", "structure_type": "predicted",
                   "source_url": f"https://alphafold.ebi.ac.uk/entry/{model.get('modelEntityId') or model.get('entryId') or accession}",
                   "confidence": model.get("globalMetricValue"), "organism": model.get("organismScientificName", "")}
    target_info = _target_metadata({**(payload.get("target") or {}), **target_info}, target_info["pdb_id"], "AlphaFold DB", "predicted")
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        job.protein_path = str(pdbqt)
        job.settings = {**job.settings, "protein": {**target_info, "target": target_info, "raw_pdb": str(raw), "cleaned_pdb": str(cleaned), "chains": inspection["chains"], "ligands_detected": inspection["ligands_detected"], "warnings": prep["warnings"]}}
        job.status = "created"
        session.commit()
    return {**target_info, "chains": inspection["chains"], "warnings": prep["warnings"]}


@router.post("/{job_id}/protein/upload")
async def upload_protein(job_id: str, request: Request):
    """Accept a raw PDB file body; avoids multipart parser dependency for desktop installs."""
    filename = request.headers.get("x-filename", "upload.pdb")
    if Path(filename).suffix.lower() not in {".pdb", ".ent"}:
        raise HTTPException(415, detail={"error": {"code": "UNSUPPORTED_PROTEIN_FORMAT", "message": "Upload a PDB or ENT structure file.", "details": {}}})
    content = await request.body()
    if not content or len(content) > settings.MAX_UPLOAD_SIZE_MB * 1024 * 1024:
        raise HTTPException(413, detail={"error": {"code": "INVALID_PROTEIN_FILE", "message": f"Protein file must be non-empty and no larger than {settings.MAX_UPLOAD_SIZE_MB} MB.", "details": {}}})
    with SessionLocal() as session:
        if not session.get(Job, job_id):
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
    root = Path(settings.DATA_DIR) / job_id / "protein"
    root.mkdir(parents=True, exist_ok=True)
    raw = root / Path(filename).name
    raw.write_bytes(content)
    try:
        svc = ProteinService()
        inspection = svc.inspect_structure(raw)
        cleaned = root / "uploaded_cleaned.pdb"
        pdbqt = root / "uploaded_cleaned.pdbqt"
        prep = svc.clean_and_prepare_receptor(raw, cleaned, pdbqt)
    except Exception as exc:
        raw.unlink(missing_ok=True)
        raise HTTPException(422, detail={"error": {"code": "PROTEIN_PREPARATION_FAILED", "message": str(exc)[:300], "details": {}}}) from exc
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        job.protein_path = str(pdbqt)
        target = _target_metadata({"name": Path(filename).stem, "source_url": ""}, Path(filename).stem, "User upload", "uploaded")
        job.settings = {**job.settings, "protein": {**target, "target": target, "raw_pdb": str(raw), "cleaned_pdb": str(cleaned), "chains": inspection["chains"], "ligands_detected": inspection["ligands_detected"], "warnings": prep["warnings"]}}
        job.status = "created"
        session.commit()
    return {**target, "chains": inspection["chains"], "ligands_detected": inspection["ligands_detected"], "warnings": prep["warnings"]}


@router.post("/{job_id}/ligands/upload")
async def upload_ligands(job_id: str, request: Request):
    """Ingest CSV (name,smiles,label) or line-oriented SMILES files."""
    filename = request.headers.get("x-filename", "ligands.csv")
    suffix = Path(filename).suffix.lower()
    if suffix not in {".csv", ".smi", ".smiles", ".txt"}:
        raise HTTPException(415, detail={"error": {"code": "UNSUPPORTED_LIGAND_FORMAT", "message": "Upload CSV or SMILES text (.smi/.smiles/.txt).", "details": {}}})
    content = await request.body()
    if not content or len(content) > settings.MAX_UPLOAD_SIZE_MB * 1024 * 1024:
        raise HTTPException(413, detail={"error": {"code": "INVALID_LIGAND_FILE", "message": f"Ligand file must be non-empty and no larger than {settings.MAX_UPLOAD_SIZE_MB} MB.", "details": {}}})
    raw = content.decode("utf-8-sig", errors="replace")
    try:
        if suffix == ".csv":
            reader = csv.DictReader(io.StringIO(raw))
            if not reader.fieldnames or not ({"smiles", "SMILES"} & set(reader.fieldnames)):
                raise ValueError("CSV needs a smiles column and may include name and label columns.")
            rows = [{"name": row.get("name") or row.get("Name") or f"ligand_{i+1}", "smiles": row.get("smiles") or row.get("SMILES") or "", "label": row.get("label") or "unknown"} for i, row in enumerate(reader)]
        else:
            rows = []
            for i, line in enumerate(raw.splitlines()):
                line = line.strip()
                if not line or line.startswith("#"): continue
                bits = line.split()
                rows.append({"smiles": bits[0], "name": bits[1] if len(bits) > 1 else f"ligand_{i+1}", "label": "unknown"})
        payload = LigandsRequest(ligands=rows)
    except Exception as exc:
        raise HTTPException(422, detail={"error": {"code": "INVALID_LIGAND_FILE", "message": str(exc)[:300], "details": {}}}) from exc
    return await set_ligands(job_id, payload)


@router.get("/{job_id}/structure")
async def get_structure(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job: raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        path = job.settings.get("protein", {}).get("cleaned_pdb")
        if not path or not Path(path).is_file():
            raise HTTPException(404, detail={"error": {"code": "STRUCTURE_NOT_FOUND", "message": "Prepared PDB coordinates are unavailable.", "details": {}}})
        return Response(Path(path).read_text(errors="replace"), media_type="chemical/x-pdb")


@router.post("/{job_id}/ligands")
async def set_ligands(job_id: str, payload: LigandsRequest):
    svc = LigandService()
    output = []
    for item in payload.ligands:
        try:
            mol = svc.parse_and_standardize_smiles(item.smiles)
            output.append({**item.model_dump(), "valid": True, "properties": svc.compute_drug_likeness(mol)})
        except Exception as exc:
            output.append({**item.model_dump(), "valid": False, "error": str(exc)[:300]})
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        job.ligands = [item.model_dump() for item in payload.ligands]
        session.commit()
    return {"molecules": output}


@router.put("/{job_id}/box")
async def set_box(job_id: str, payload: BoxModel):
    box = BoxService().validate_manual_box(payload.center, payload.size)
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        job.box = box
        session.commit()
    return box


@router.get("/{job_id}/pockets")
async def get_pockets(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        raw = job.settings.get("protein", {}).get("raw_pdb")
        if not raw or not Path(raw).exists():
            raise HTTPException(409, detail={"error": {"code": "PROTEIN_REQUIRED", "message": "Set a protein before requesting pocket suggestions.", "details": {}}})
        result = BoxService().detect_pockets_fpocket(Path(raw))
        return {**result, "crystal_ligand_box": job.box}


@router.put("/{job_id}/settings")
async def set_settings(job_id: str, payload: SettingsRequest):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        job.settings = {**job.settings, **payload.model_dump()}
        session.commit()
    return {"settings": payload.model_dump()}


@router.post("/{job_id}/qc")
async def run_qc(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        if not job.protein_path or not job.box or not job.ligands:
            raise HTTPException(409, detail={"error": {"code": "INCOMPLETE_JOB", "message": "Set a protein, ligands, and search box before QC.", "details": {}}})
        raw_protein = job.settings.get("protein", {}).get("cleaned_pdb")
        report = QCService().run(Path(raw_protein), job.ligands, job.box)
        job.qc_report = report
        job.status = "ready" if not report["requires_acknowledgement"] else "qc"
        session.commit()
    return report


@router.get("/{job_id}/qc")
async def get_qc(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        return job.qc_report or {"summary": "QC has not been run.", "checks": []}


@router.post("/{job_id}/run")
async def run_job(job_id: str, payload: RunRequest):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        if not job.qc_report:
            raise HTTPException(409, detail={"error": {"code": "QC_REQUIRED", "message": "Run pre-docking QC first.", "details": {}}})
        if job.qc_report.get("requires_acknowledgement") and not payload.acknowledge_qc_failures:
            raise HTTPException(409, detail={"error": {"code": "QC_ACK_REQUIRED", "message": "Acknowledge the QC failures before docking.", "details": {}}})
        if not job.protein_path or not job.box:
            raise HTTPException(409, detail={"error": {"code": "INCOMPLETE_JOB", "message": "A prepared protein and search box are required.", "details": {}}})
        job.status = "queued"
        job.settings = {**job.settings, "qc_failures_acknowledged": bool(payload.acknowledge_qc_failures)}
        job.progress_total = len(job.ligands)
        job.progress_done = 0
        session.commit()
    if os.getenv("DOCKLAB_USE_CELERY", "0").lower() not in {"1", "true", "yes"}:
        threading.Thread(target=_run_local_job, args=(job_id,), daemon=True, name=f"docklab-{job_id[:8]}").start()
        return {"job_id": job_id, "status": "queued", "execution": "local-background-thread"}
    task = chord([dock_molecule.s(job_id, i) for i in range(len(job.ligands))], aggregate_job.s(job_id)).apply_async()
    return {"job_id": job_id, "status": "queued", "task_id": task.id, "execution": "celery"}


@router.get("/{job_id}")
async def get_job(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        return _job_data(job)


@router.post("/{job_id}/cancel")
async def cancel_job(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        job.status = "cancelled"
        session.commit()
    return {"id": job_id, "status": "cancelled"}


@router.get("/{job_id}/results")
async def get_results(job_id: str):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        rows = _available_job_results(job)
        analysis = job.settings.get("analysis") or _available_screening_metrics(rows)
        return {"results": rows, "analysis": analysis,
                "is_partial": job.status not in {"completed", "failed", "cancelled"},
                "progress": {"done": job.progress_done, "total": job.progress_total},
                "status": job.status}


@router.get("/{job_id}/molecules/{molecule_id}/poses/{pose_number}")
async def get_pose(job_id: str, molecule_id: str, pose_number: int = 1):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        matches = [r for r in job.results if str(r.get("name")) == molecule_id]
        if not matches or not matches[0].get("sdf"):
            raise HTTPException(404, detail={"error": {"code": "POSE_NOT_FOUND", "message": "Pose was not found.", "details": {}}})
        path = Path(matches[0]["sdf"]).resolve()
        if Path(settings.DATA_DIR).resolve() not in path.parents:
            raise HTTPException(403, detail={"error": {"code": "INVALID_PATH", "message": "Pose path is outside job storage.", "details": {}}})
        return Response(path.read_text(errors="replace"), media_type="chemical/x-mdl-sdfile")


@router.get("/{job_id}/export")
async def export_job(job_id: str, type: str = "csv"):
    with SessionLocal() as session:
        job = session.get(Job, job_id)
        if not job:
            raise HTTPException(404, detail={"error": {"code": "JOB_NOT_FOUND", "message": "Job was not found.", "details": {}}})
        rows = _available_job_results(job)
        if type == "html":
            analysis = job.settings.get("analysis") or _available_screening_metrics(rows)
            original_results = job.results
            try:
                job.results = rows
                if analysis and not job.settings.get("analysis"):
                    job.settings = {**job.settings, "analysis": analysis}
                report = ReportService.render_html(job)
            finally:
                job.results = original_results
            return Response(report, media_type="text/html; charset=utf-8", headers={"Content-Disposition": "inline; filename=docklab-report.html"})
        if type == "pdf":
            try:
                original_results = job.results
                try:
                    job.results = rows
                    analysis = job.settings.get("analysis") or _available_screening_metrics(rows)
                    if analysis and not job.settings.get("analysis"):
                        job.settings = {**job.settings, "analysis": analysis}
                    pdf = ReportService.render_pdf(job)
                finally:
                    job.results = original_results
            except Exception as exc:
                raise HTTPException(501, detail={"error": {"code": "PDF_RENDER_FAILED", "message": "PDF rendering is unavailable in this environment.", "details": {"reason": str(exc)[:200]}}}) from exc
            return Response(pdf, media_type="application/pdf", headers={"Content-Disposition": "attachment; filename=docklab-report.pdf"})
        if type == "settings":
            return Response(json.dumps({"settings": job.settings, "box": job.box}, indent=2), media_type="application/json", headers={"Content-Disposition": "attachment; filename=settings.json"})
        if type == "sdf":
            poses = [Path(row["sdf"]) for row in rows if row.get("sdf") and Path(row["sdf"]).exists()]
            if not poses:
                return Response("", media_type="chemical/x-mdl-sdfile", headers={"Content-Disposition": "attachment; filename=poses.sdf"})
            return Response("\n".join(path.read_text(errors="replace") for path in poses), media_type="chemical/x-mdl-sdfile", headers={"Content-Disposition": "attachment; filename=poses.sdf"})
        if type == "csv":
            buf = io.StringIO()
            fields = sorted({k for row in rows for k in row})
            writer = csv.DictWriter(buf, fieldnames=fields or ["name", "status"])
            writer.writeheader()
            for row in rows:
                writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})
            return Response(buf.getvalue(), media_type="text/csv", headers={"Content-Disposition": "attachment; filename=results.csv"})
        if type == "zip":
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("settings.json", json.dumps({"settings": job.settings, "box": job.box}, indent=2))
                table = io.StringIO()
                fields = sorted({k for row in rows for k in row}) or ["name", "status"]
                writer = csv.DictWriter(table, fieldnames=fields); writer.writeheader()
                for row in rows: writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})
                archive.writestr("results.csv", table.getvalue())
                archive.writestr("qc.json", json.dumps(job.qc_report, indent=2))
                for row in rows:
                    if row.get("sdf") and Path(row["sdf"]).exists(): archive.write(row["sdf"], f"poses/{Path(row['sdf']).name}")
            return Response(buf.getvalue(), media_type="application/zip", headers={"Content-Disposition": "attachment; filename=docklab-export.zip"})
        raise HTTPException(422, detail={"error": {"code": "UNSUPPORTED_EXPORT", "message": "Supported export types: csv, sdf, pdf, settings, zip.", "details": {}}})
