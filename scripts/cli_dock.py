#!/usr/bin/env python3
"""
scripts/cli_dock.py
Command-Line Science Pipeline for DockLab.

Executes protein prep, crystal ligand extraction, box generation, redocking validation,
and batch small-molecule docking without web dependencies.
"""

import sys
import os
import argparse
import csv
import json
import time
from pathlib import Path

# Add project root to path
ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

from backend.app.services.protein_service import ProteinService
from backend.app.services.ligand_service import LigandService
from backend.app.services.box_service import BoxService
from backend.app.services.docking_service import DockingService
from backend.app.services.validation_service import ValidationService
from backend.app.services.qc_service import QCService
from backend.app.services.reliability_service import ReliabilityService
from backend.app.services.analysis_service import AnalysisService

def main():
    started = time.perf_counter()
    parser = argparse.ArgumentParser(description="DockLab Science Core CLI Docking Pipeline")
    parser.add_argument("--pdb", type=str, required=True, help="PDB ID (e.g. 1HSG) or path to PDB file")
    parser.add_argument("--ligands", type=str, required=True, help="Path to ligands CSV (name, smiles, label)")
    parser.add_argument("--out", type=str, default="results", help="Output directory for poses and results CSV")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina search exhaustiveness (default: 8)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)")
    parser.add_argument("--ligand-resname", type=str, default="MK1", help="Residue name of co-crystallized ligand")
    parser.add_argument("--qc", action="store_true", help="Run pre-docking QC and print its checklist")
    parser.add_argument("--reliability", action="store_true", help="Repeat the redocking run with multiple seeds")
    parser.add_argument("--seeds", type=int, default=5, help="Number of independent runs for reliability (2-10)")
    parser.add_argument("--topk", type=int, default=10, help="Number of top-ranked ligands for reliability (max 20)")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "prep"
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 75)
    print("  DockLab Science Core CLI Docking Engine")
    print(f"  Target: {args.pdb} | Ligands: {args.ligands} | Out: {args.out}")
    print(f"  Exhaustiveness: {args.exhaustiveness} | Seed: {args.seed}")
    print("=" * 75)

    # 1. Initialize Services
    protein_svc = ProteinService()
    ligand_svc = LigandService()
    box_svc = BoxService()
    docking_svc = DockingService()
    validation_svc = ValidationService()
    qc = None
    rel = None

    if not 2 <= args.seeds <= 10 or not 1 <= args.topk <= 20:
        parser.error("--seeds must be 2-10 and --topk must be 1-20")

    # 2. Fetch and Inspect Protein
    print("\n[Step 1/5] Ingesting and Preparing Receptor...")
    raw_pdb = protein_svc.fetch_or_load_pdb(args.pdb, cache_dir)
    inspection = protein_svc.inspect_structure(raw_pdb)
    print(f"  Receptor Chains: {inspection['chains']} | Atoms: {inspection['total_atoms']}")
    print(f"  Heteroatoms detected: {[l['resname'] for l in inspection['ligands_detected']]}")

    # 3. Extract Crystal Ligand directly from experimental coordinates
    crystal_pdb = cache_dir / "crystal_ligand.pdb"
    crystal_pdbqt = cache_dir / "crystal_ligand.pdbqt"
    crystal_sdf = cache_dir / "crystal_ligand.sdf"

    print(f"\n[Step 2/5] Extracting Co-Crystallized Ligand '{args.ligand_resname}' from Crystal Coordinates...")
    try:
        protein_svc.extract_crystal_ligand(raw_pdb, args.ligand_resname, crystal_pdb)
        ligand_svc.prepare_crystal_ligand_pdbqt(crystal_pdb, crystal_pdbqt)
        # Convert to SDF for ProLIF
        subprocess_obabel = protein_svc._get_obabel_cmd()
        import subprocess
        subprocess.run([subprocess_obabel, "-ipdbqt", str(crystal_pdbqt), "-osdf", "-O", str(crystal_sdf)],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("  Crystal ligand extracted and parameterized to PDBQT.")
    except Exception as exc:
        print(f"  WARNING: Could not extract crystal ligand: {exc}")
        crystal_pdbqt = None

    # 4. Clean Receptor with PDBFixer
    cleaned_pdb = cache_dir / f"{raw_pdb.stem}_cleaned.pdb"
    receptor_pdbqt = cache_dir / f"{raw_pdb.stem}_cleaned.pdbqt"
    prep_info = protein_svc.clean_and_prepare_receptor(raw_pdb, cleaned_pdb, receptor_pdbqt)
    print(f"  Receptor cleaned & protonated at pH 7.4 -> {receptor_pdbqt.name}")
    if prep_info["warnings"]:
        for w in prep_info["warnings"]:
            print(f"  [QC Note] {w}")

    # 5. Define Search Bounding Box
    print("\n[Step 3/5] Computing Active Site Grid Box...")
    if crystal_pdb and crystal_pdb.exists():
        box = box_svc.box_from_crystal_ligand(crystal_pdb)
    else:
        # Fallback to active site catalytic Asp25 residues
        box = box_svc.box_from_residues(cleaned_pdb, [("A", 25), ("B", 25)])
    
    print(f"  Box Center: X={box['center']['x']}, Y={box['center']['y']}, Z={box['center']['z']}")
    print(f"  Box Size:   X={box['size']['x']}, Y={box['size']['y']}, Z={box['size']['z']}")

    ligand_rows = []
    with open(args.ligands, "r", encoding="utf-8", newline="") as f:
        ligand_rows = list(csv.DictReader(f))
    if args.qc:
        qc = QCService().run(cleaned_pdb, ligand_rows, box, crystal_pdb if crystal_pdb.exists() else None)
        print("\n[Pre-docking QC]")
        for check in qc["checks"]:
            print(f"  {check['severity'].upper():4s} {check['title']}: {check['detail']}")
        print("  " + qc["summary"])

    # 6. Redocking Validation
    redock_result = None
    if crystal_pdbqt and crystal_pdbqt.exists():
        print("\n[Step 4/5] Executing Redocking Validation of Crystal Ligand...")
        redock_pdbqt = out_dir / "redock_pose.pdbqt"
        redock_sdf = out_dir / "redock_pose.sdf"
        
        redock_dock = docking_svc.dock(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=crystal_pdbqt,
            box_center=box["center"],
            box_size=box["size"],
            output_pdbqt=redock_pdbqt,
            output_sdf=redock_sdf,
            exhaustiveness=args.exhaustiveness,
            seed=args.seed
        )
        
        rmsd_info = validation_svc.compute_symmetric_rmsd(crystal_sdf, redock_sdf)
        
        # Key residue recovery
        recovery_info = None
        try:
            recovery_info = validation_svc.compute_key_residue_recovery(cleaned_pdb, crystal_sdf, redock_sdf)
            recovery_pct = recovery_info["recovery_percentage"]
            common_count = recovery_info["common_residues_count"]
            total_count = recovery_info["crystal_residues_count"]
        except Exception as e:
            print(f"  Key-residue recovery unavailable: {e}")

        redock_result = {
            "score": redock_dock["best_score"],
            "rmsd": rmsd_info["rmsd_angstroms"],
            "status": rmsd_info["status"],
            "recovery_pct": recovery_pct if recovery_info else None,
            "recovery_str": f"{common_count}/{total_count}" if recovery_info else None
        }

        print(f"  >> Redocked Best Score: {redock_result['score']:.2f} kcal/mol")
        print(f"  >> Heavy-Atom RMSD:     {redock_result['rmsd']:.3f} Angstroms [{redock_result['status']}]")
        if redock_result["recovery_pct"] is not None:
            print(f"  >> Key-Residue Recovery:{redock_result['recovery_pct']:.1f}% ({redock_result['recovery_str']} residues)")
        print(f"  >> Validation Status:   {'PASSED (<2.0 A)' if rmsd_info['passed'] else 'FAILED'}")

    # 7. Batch Small-Molecule Docking
    print("\n[Step 5/5] Docking Candidate Small Molecules...")
    print(f"  Loaded {len(ligand_rows)} candidate molecules from {args.ligands}")

    results_records = []
    poses_dir = out_dir / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)

    for idx, lig in enumerate(ligand_rows, 1):
        name = lig.get("name", f"mol_{idx}")
        smiles = lig.get("smiles", "").strip()
        label = lig.get("label", "unknown")
        print(f"  [{idx}/{len(ligand_rows)}] Preparing & docking: {name:<20}", end="", flush=True)

        try:
            cand_pdbqt = cache_dir / f"{name}.pdbqt"
            cand_sdf = poses_dir / f"{name}_pose.sdf"
            cand_dock_pdbqt = poses_dir / f"{name}_docked.pdbqt"

            prep_res = ligand_svc.prepare_ligand_from_smiles(
                smiles=smiles,
                output_pdbqt=cand_pdbqt,
                seed=args.seed
            )
            props = prep_res["properties"]

            dock_res = docking_svc.dock(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=cand_pdbqt,
                box_center=box["center"],
                box_size=box["size"],
                output_pdbqt=cand_dock_pdbqt,
                output_sdf=cand_sdf,
                exhaustiveness=args.exhaustiveness,
                seed=args.seed
            )

            score = dock_res["best_score"]
            heavy_atoms = props["heavy_atom_count"]
            le = round(score / heavy_atoms, 3) if heavy_atoms > 0 else 0.0

            results_records.append({
                "rank": 0,  # Will assign after sort
                "name": name,
                "label": label,
                "score_kcal_mol": score,
                "ligand_efficiency": le,
                "mw": props["mw"],
                "logp": props["logp"],
                "hbd": props["hbd"],
                "hba": props["hba"],
                "rot_bonds": props["rotatable_bonds"],
                "lipinski_violations": props["lipinski_violations"],
                "pains_alerts": len(props["pains_alerts"]),
                "brenk_alerts": len(props["brenk_alerts"]),
                "pose_path": str(cand_dock_pdbqt.relative_to(out_dir)),
                "prepared_path": str(cand_pdbqt),
                "pose_sdf_path": str(cand_sdf),
                "smiles": smiles,
                "status": "docked"
            })
            print(f" -> Score: {score:6.2f} kcal/mol | LE: {le:5.2f}")
        except Exception as exc:
            print(f" -> FAILED: {exc}")
            results_records.append({
                "rank": 999,
                "name": name,
                "label": label,
                "score_kcal_mol": None,
                "ligand_efficiency": None,
                "mw": None,
                "logp": None,
                "hbd": None,
                "hba": None,
                "rot_bonds": None,
                "lipinski_violations": None,
                "pains_alerts": None,
                "brenk_alerts": None,
                "pose_path": "",
                "smiles": smiles,
                "status": f"failed: {exc}"
            })

    # Sort results by score (more negative is better)
    docked_only = [r for r in results_records if r["status"] == "docked"]
    failed_only = [r for r in results_records if r["status"] != "docked"]
    docked_only.sort(key=lambda x: x["score_kcal_mol"])

    for rank_idx, record in enumerate(docked_only, 1):
        record["rank"] = rank_idx

    final_results = docked_only + failed_only

    # Optional score ranking metrics use measured scores and labels from the CSV.
    labeled = [r for r in docked_only if r.get("label") in {"active", "decoy"}]
    active_count = sum(r["label"] == "active" for r in labeled)
    decoy_count = sum(r["label"] == "decoy" for r in labeled)
    metrics = None
    if active_count >= 3 and decoy_count >= 3:
        metrics = AnalysisService().screening_metrics(
            [int(r["label"] == "active") for r in labeled],
            [r["score_kcal_mol"] for r in labeled],
        )
        print(f"\n  Screening metrics: ROC-AUC={metrics['roc_auc']:.3f}; EF5={metrics['ef5']:.2f}")

    reliability_report = None
    if args.reliability and docked_only:
        rel_dir = out_dir / "reliability"
        seed_numbers = list(range(args.seed, args.seed + args.seeds))
        seed_scores = {}
        for record in docked_only[:args.topk]:
            rel = ReliabilityService(docking_svc).run_ligand(
                receptor_pdbqt, Path(record["prepared_path"]), box["center"], box["size"], rel_dir / record["name"],
                seeds=seed_numbers, exhaustiveness=args.exhaustiveness,
                initial_score=float(record["score_kcal_mol"]), initial_pose_sdf=Path(record["pose_sdf_path"]),
            )
            seed_scores[record["name"]] = rel["scores"]
            record["reliability"] = rel
        if seed_scores:
            ranks = ReliabilityService.rank_stability(seed_scores)
            for record in docked_only[:args.topk]:
                if record.get("reliability"):
                    record["reliability"]["rank_range"] = ranks["rank_ranges"].get(record["name"])
                    record["reliability"]["mean_spearman"] = ranks["mean_spearman"]
            reliability_report = {"top_k": min(args.topk, len(seed_scores)), "n_seeds": args.seeds,
                "mean_spearman": ranks["mean_spearman"],
                "ligands": {r["name"]: r["reliability"] for r in docked_only[:args.topk] if r.get("reliability")}}
            print(f"  Top-ligand reliability: {len(seed_scores)} ligands × {args.seeds} seeds | mean rank Spearman={ranks['mean_spearman']:.3f}")

    # Write results CSV
    results_csv = out_dir / "results.csv"
    fieldnames = [
        "rank", "name", "label", "score_kcal_mol", "ligand_efficiency",
        "mw", "logp", "hbd", "hba", "rot_bonds", "lipinski_violations",
        "pains_alerts", "brenk_alerts", "status", "pose_path", "smiles"
        , "reliability_label", "reliability_score_sd", "pose_consistency", "rank_range"
    ]
    with open(results_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for record in final_results:
            rel = record.get("reliability") or {}
            record.update({"reliability_label": rel.get("label"), "reliability_score_sd": rel.get("score_sd"),
                "pose_consistency": rel.get("pose_consistency"), "rank_range": json.dumps(rel.get("rank_range")) if rel.get("rank_range") is not None else None})
            writer.writerow({k: record.get(k) for k in fieldnames})

    # Print Final Summary Table
    print("\n" + "=" * 80)
    print("  FINAL DOCKING RESULTS (Ranked by Vina Score)")
    print("=" * 80)
    print(f" {'Rank':<4} | {'Name':<16} | {'Label':<7} | {'Score':<8} | {'LE':<6} | {'MW':<7} | {'LogP':<6} | {'Ro5 Viol':<8}")
    print("-" * 80)
    for r in docked_only[:15]:
        print(f" {r['rank']:<4} | {r['name']:<16} | {r['label']:<7} | {r['score_kcal_mol']:6.2f}   | {r['ligand_efficiency']:5.2f}  | {r['mw']:6.1f}  | {r['logp']:5.2f}  | {r['lipinski_violations']:<8}")

    print("=" * 80)
    if redock_result:
        recovery_text = f"{redock_result['recovery_pct']:.1f}%" if redock_result["recovery_pct"] is not None else "unavailable"
        print(f"  REDOCKING VALIDATION: RMSD = {redock_result['rmsd']:.3f} A | Score = {redock_result['score']:.2f} kcal/mol | Key Residues = {recovery_text}")
    print(f"  Results saved to: {results_csv.resolve()}")
    report = {"target": args.pdb, "ligand_count": len(ligand_rows), "docked_count": len(docked_only), "failed_count": len(failed_only),
              "redocking": redock_result, "screening": metrics if active_count >= 3 and decoy_count >= 3 else None,
              "qc": qc, "reliability": reliability_report, "runtime_seconds": round(time.perf_counter() - started, 2),
              "seed": args.seed, "exhaustiveness": args.exhaustiveness}
    metrics_path = out_dir / "metrics.json"
    metrics_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  Metrics saved to: {metrics_path.resolve()}")
    print("=" * 80)

if __name__ == "__main__":
    main()

