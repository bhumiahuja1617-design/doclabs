#!/usr/bin/env python3
"""
scripts/finish_batch.py
Docks the remaining ligands (indinavir, lopinavir, atazanavir) that failed or
were not completed in the first run, then assembles a full results.csv.

Uses exhaustiveness=4 (Fast preset) with 180s timeout — sufficient for real
scores while keeping runtime practical. All existing poses are reused.
"""
import sys
import csv
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from backend.app.services.ligand_service import LigandService
from backend.app.services.docking_service import DockingService

SEED = 42
EXHAUSTIVENESS = 4   # Fast preset — still real Vina scores
TIMEOUT = 180        # seconds

OUT_DIR = ROOT / "results"
PREP_DIR = OUT_DIR / "prep"
POSES_DIR = OUT_DIR / "poses"
POSES_DIR.mkdir(parents=True, exist_ok=True)

RECEPTOR_PDBQT = PREP_DIR / "1hsg_cleaned.pdbqt"
BOX = {
    "center": {"x": 13.073, "y": 22.467, "z": 5.557},
    "size":   {"x": 24.854, "y": 20.214, "z": 18.635}
}

# Already-docked results from task-406 log
EXISTING = {
    "darunavir":  {"score_kcal_mol": -8.88,  "ligand_efficiency": -0.23, "mw": 547.7, "logp": 3.4,  "hbd": 3, "hba": 8, "rot_bonds": 10, "lipinski_violations": 1, "pains_alerts": 0, "brenk_alerts": 0, "label": "active", "status": "docked"},
    "saquinavir": {"score_kcal_mol": -10.83, "ligand_efficiency": -0.22, "mw": 670.8, "logp": 3.8,  "hbd": 4, "hba": 7, "rot_bonds": 11, "lipinski_violations": 2, "pains_alerts": 0, "brenk_alerts": 0, "label": "active", "status": "docked"},
    "ritonavir":  {"score_kcal_mol": -9.79,  "ligand_efficiency": -0.19, "mw": 720.9, "logp": 5.0,  "hbd": 3, "hba": 8, "rot_bonds": 14, "lipinski_violations": 2, "pains_alerts": 0, "brenk_alerts": 0, "label": "active", "status": "docked"},
    "nelfinavir": {"score_kcal_mol": -10.45, "ligand_efficiency": -0.26, "mw": 567.8, "logp": 4.4,  "hbd": 4, "hba": 6, "rot_bonds": 8,  "lipinski_violations": 1, "pains_alerts": 0, "brenk_alerts": 0, "label": "active", "status": "docked"},
}

# Ligands still needing docking
PENDING = [
    {"name": "indinavir",  "label": "active", "smiles": "CC(C)(C)NC(=O)C1CN(CC1)Cc1ccc2ccccc2n1.OC(=O)C(O)(CC(=O)O)CC(=O)O"},
    {"name": "lopinavir",  "label": "active", "smiles": "CC1=NC=C(C(=O)NC(CC(C)C)C(=O)NC2CC3CCCC(C3)C2)C=C1.CC(C)CC1NC(=O)C(Cc2ccccc2)NC(=O)c2ccc(N)cc2"},
    {"name": "atazanavir", "label": "active", "smiles": "COC(=O)NC(C(=O)NC(CC1=CC=CC=C1)C(=O)NC(C(C)(C)C)C(=O)NC2=CC=C(C=C2)C3=CC=NC=C3)CC4=CC=CC=C4"},
]

# Read actives CSV for SMILES (more reliable than hard-coded)
ACTIVES_CSV = ROOT / "examples" / "1hsg" / "actives.csv"
smiles_map = {}
with open(ACTIVES_CSV) as f:
    for row in csv.DictReader(f):
        smiles_map[row["name"]] = row["smiles"]

ligand_svc = LigandService()
docking_svc = DockingService()

new_results = []

for lig in PENDING:
    name = lig["name"]
    # Try to get SMILES from actives.csv first
    smiles = smiles_map.get(name, lig.get("smiles", ""))
    label = lig["label"]
    cand_pdbqt = PREP_DIR / f"{name}.pdbqt"
    cand_dock_pdbqt = POSES_DIR / f"{name}_docked.pdbqt"
    cand_sdf = POSES_DIR / f"{name}_pose.sdf"

    print(f"\n--- Docking {name} (exhaustiveness={EXHAUSTIVENESS}, timeout={TIMEOUT}s) ---")

    # Prepare if pdbqt not already there
    if not cand_pdbqt.exists():
        print(f"  Preparing ligand from SMILES...")
        try:
            prep = ligand_svc.prepare_ligand_from_smiles(smiles=smiles, output_pdbqt=cand_pdbqt, seed=SEED)
            props = prep["properties"]
        except Exception as e:
            print(f"  PREP FAILED: {e}")
            new_results.append({
                "name": name, "label": label,
                "score_kcal_mol": 0.0, "ligand_efficiency": 0.0,
                "mw": 0, "logp": 0, "hbd": 0, "hba": 0, "rot_bonds": 0,
                "lipinski_violations": 0, "pains_alerts": 0, "brenk_alerts": 0,
                "status": f"failed_prep: {e}", "pose_path": "", "smiles": smiles
            })
            continue
    else:
        print(f"  Using existing {cand_pdbqt.name}")
        # We still need props, re-run prep to get them
        try:
            prep = ligand_svc.prepare_ligand_from_smiles(smiles=smiles, output_pdbqt=cand_pdbqt, seed=SEED)
            props = prep["properties"]
        except Exception:
            props = {"heavy_atom_count": 0, "mw": 0, "logp": 0, "hbd": 0, "hba": 0,
                     "rotatable_bonds": 0, "lipinski_violations": 0, "pains_alerts": [], "brenk_alerts": []}

    # Dock
    try:
        dock_res = docking_svc.dock(
            receptor_pdbqt=RECEPTOR_PDBQT,
            ligand_pdbqt=cand_pdbqt,
            box_center=BOX["center"],
            box_size=BOX["size"],
            output_pdbqt=cand_dock_pdbqt,
            output_sdf=cand_sdf,
            exhaustiveness=EXHAUSTIVENESS,
            seed=SEED,
            timeout_seconds=TIMEOUT
        )
        score = dock_res["best_score"]
        heavy_atoms = props.get("heavy_atom_count", 1) or 1
        le = round(score / heavy_atoms, 3)
        print(f"  Score: {score:.2f} kcal/mol | LE: {le:.3f}")
        new_results.append({
            "name": name, "label": label,
            "score_kcal_mol": score, "ligand_efficiency": le,
            "mw": props.get("mw", 0), "logp": props.get("logp", 0),
            "hbd": props.get("hbd", 0), "hba": props.get("hba", 0),
            "rot_bonds": props.get("rotatable_bonds", 0),
            "lipinski_violations": props.get("lipinski_violations", 0),
            "pains_alerts": len(props.get("pains_alerts", [])),
            "brenk_alerts": len(props.get("brenk_alerts", [])),
            "status": "docked", "pose_path": str(cand_dock_pdbqt.relative_to(OUT_DIR)),
            "smiles": smiles
        })
    except Exception as e:
        print(f"  DOCK FAILED: {e}")
        new_results.append({
            "name": name, "label": label,
            "score_kcal_mol": 0.0, "ligand_efficiency": 0.0,
            "mw": props.get("mw", 0), "logp": props.get("logp", 0),
            "hbd": props.get("hbd", 0), "hba": props.get("hba", 0),
            "rot_bonds": props.get("rotatable_bonds", 0),
            "lipinski_violations": props.get("lipinski_violations", 0),
            "pains_alerts": len(props.get("pains_alerts", [])),
            "brenk_alerts": len(props.get("brenk_alerts", [])),
            "status": f"failed: {e}", "pose_path": "", "smiles": smiles
        })

# Merge with existing results
all_records = []
for name, data in EXISTING.items():
    cand_dock_pdbqt = POSES_DIR / f"{name}_docked.pdbqt"
    all_records.append({
        "name": name, "label": data["label"],
        "score_kcal_mol": data["score_kcal_mol"],
        "ligand_efficiency": data["ligand_efficiency"],
        "mw": data["mw"], "logp": data["logp"],
        "hbd": data["hbd"], "hba": data["hba"],
        "rot_bonds": data["rot_bonds"],
        "lipinski_violations": data["lipinski_violations"],
        "pains_alerts": data["pains_alerts"],
        "brenk_alerts": data["brenk_alerts"],
        "status": "docked",
        "pose_path": str(cand_dock_pdbqt.relative_to(OUT_DIR)),
        "smiles": smiles_map.get(name, "")
    })

all_records.extend(new_results)

# Sort docked ones by score
docked = sorted([r for r in all_records if r["status"] == "docked"],
                key=lambda x: x["score_kcal_mol"])
failed = [r for r in all_records if r["status"] != "docked"]

for idx, r in enumerate(docked, 1):
    r["rank"] = idx
for r in failed:
    r["rank"] = 999

final = docked + failed

# Write results.csv
FIELDNAMES = ["rank","name","label","score_kcal_mol","ligand_efficiency",
              "mw","logp","hbd","hba","rot_bonds","lipinski_violations",
              "pains_alerts","brenk_alerts","status","pose_path","smiles"]

out_csv = OUT_DIR / "results.csv"
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=FIELDNAMES)
    w.writeheader()
    w.writerows(final)

print("\n" + "=" * 70)
print("FINAL DOCKING RESULTS (Ranked by Vina Score)")
print("=" * 70)
print(f" {'Rank':<4} | {'Name':<12} | {'Score':>8} | {'LE':>6} | {'MW':>7} | Ro5")
print("-" * 60)
for r in docked:
    print(f" {r['rank']:<4} | {r['name']:<12} | {r['score_kcal_mol']:>8.2f} | {r['ligand_efficiency']:>6.3f} | {r['mw']:>7.1f} | {r['lipinski_violations']}")
if failed:
    print("\nFailed ligands:")
    for r in failed:
        print(f"  {r['name']}: {r['status']}")
print("=" * 70)
print(f"Results written to: {out_csv}")

