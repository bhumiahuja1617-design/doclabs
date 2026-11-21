"""Static example target discovery and example data retrieval."""
import csv
import json
import shutil
import uuid
from datetime import datetime, timedelta, timezone
from fastapi import APIRouter, HTTPException
from backend.app.config import settings
from backend.app.db.models import Job
from backend.app.db.session import SessionLocal
from backend.app.services.box_service import BoxService

router = APIRouter(prefix="/api/v1/examples", tags=["Examples"])


@router.get("")
async def list_examples():
    root = settings.BASE_DIR / "examples"
    result = []
    for target in root.iterdir():
        card = target / "target.json"
        if target.is_dir() and card.exists():
            result.append({"name": target.name, "target": json.loads(card.read_text(encoding="utf-8"))})
    return {"examples": result}


@router.get("/{name}")
async def get_example(name: str):
    root = (settings.BASE_DIR / "examples" / name).resolve()
    if settings.BASE_DIR.resolve() not in root.parents or not root.is_dir():
        raise HTTPException(404, detail={"error": {"code": "EXAMPLE_NOT_FOUND", "message": "Example was not found.", "details": {}}})
    card = root / "target.json"
    if not card.exists():
        raise HTTPException(404, detail={"error": {"code": "EXAMPLE_NOT_FOUND", "message": "Example was not found.", "details": {}}})
    ligands = []
    for filename in ("actives.csv", "decoys.csv"):
        path = root / filename
        if path.exists():
            with path.open(encoding="utf-8", newline="") as stream:
                ligands.extend(list(csv.DictReader(stream)))
    return {"name": name, "target": json.loads(card.read_text(encoding="utf-8")), "ligands": ligands}


@router.post("/{name}/load")
async def load_example(name: str):
    """Create a runnable draft using the bundled prepared protein and example ligands."""
    root = (settings.BASE_DIR / "examples" / name).resolve()
    if settings.BASE_DIR.resolve() not in root.parents or not (root / "target.json").exists():
        raise HTTPException(404, detail={"error": {"code": "EXAMPLE_NOT_FOUND", "message": "Example was not found.", "details": {}}})
    target = json.loads((root / "target.json").read_text(encoding="utf-8"))
    data = settings.DATA_DIR / str(uuid.uuid4())
    protein_dir = data / "protein"
    protein_dir.mkdir(parents=True, exist_ok=True)
    pdb_name = target.get("pdb_id", name).lower()
    candidates = [root / f"{pdb_name}_cleaned.pdbqt", root / f"{name}_cleaned.pdbqt"]
    prepared = next((path for path in candidates if path.exists()), None)
    if prepared is None:
        raise HTTPException(422, detail={"error": {"code": "EXAMPLE_NOT_PREPARED", "message": "This example has no prepared receptor file.", "details": {}}})
    receptor = protein_dir / prepared.name
    shutil.copy2(prepared, receptor)
    cleaned_candidates = [root / f"{pdb_name}_cleaned.pdb", root / f"{name}_cleaned.pdb"]
    cleaned = next((path for path in cleaned_candidates if path.exists()), None)
    if cleaned:
        shutil.copy2(cleaned, protein_dir / cleaned.name)
    ligand_file = []
    for filename in ("actives.csv", "decoys.csv"):
        path = root / filename
        if path.exists():
            with path.open(encoding="utf-8", newline="") as stream:
                ligand_file.extend(list(csv.DictReader(stream)))
    crystal = root / "crystal_ligand.pdb"
    if not crystal.exists():
        raise HTTPException(422, detail={"error": {"code": "EXAMPLE_MISSING_LIGAND", "message": "Crystal ligand coordinates are missing.", "details": {}}})
    box = BoxService().box_from_crystal_ligand(crystal)
    job = Job(id=data.name, protein_path=str(receptor), box=box, ligands=ligand_file,
        settings={"protein": {"pdb_id": target.get("pdb_id", name), "cleaned_pdb": str(protein_dir / cleaned.name) if cleaned else None, "target": target}, "seed": 42, "exhaustiveness": 8},
        expires_at=datetime.now(timezone.utc) + timedelta(days=settings.JOB_RETENTION_DAYS))
    with SessionLocal() as session:
        session.add(job)
        session.commit()
    return {"job": {"id": job.id, "status": job.status, "box": job.box}, "target": target, "ligand_count": len(ligand_file)}
