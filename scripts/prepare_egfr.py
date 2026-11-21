"""Prepare local EGFR 1M17 receptor and its AQ4 co-crystal ligand."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.app.services.protein_service import ProteinService
from backend.app.services.box_service import BoxService

OUT = ROOT / "examples" / "egfr"
RAW = OUT / "1m17.pdb"
protein = ProteinService()
cleaned, receptor = OUT / "1m17_cleaned.pdb", OUT / "1m17_cleaned.pdbqt"
ligand = OUT / "crystal_ligand.pdb"
protein.extract_crystal_ligand(RAW, "AQ4", ligand)
prep = protein.clean_and_prepare_receptor(RAW, cleaned, receptor)
obabel = protein._get_obabel_cmd()
crystal_pdbqt = OUT / "crystal_ligand.pdbqt"
crystal_sdf = OUT / "crystal_ligand.sdf"
subprocess.run([obabel, "-ipdb", str(ligand), "-opdbqt", "-O", str(crystal_pdbqt), "-h"], check=True, capture_output=True, text=True)
subprocess.run([obabel, "-ipdbqt", str(crystal_pdbqt), "-osdf", "-O", str(crystal_sdf)], check=True, capture_output=True, text=True)
box = BoxService().box_from_crystal_ligand(ligand)
(OUT / "box.json").write_text(json.dumps(box, indent=2), encoding="utf-8")
print(json.dumps({"receptor": str(receptor), "ligand_resname": "AQ4", "box": box, "warnings": prep["warnings"]}, indent=2))
