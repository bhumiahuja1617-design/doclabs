#!/usr/bin/env python3
"""Run the HIV-1 protease inhibitor/mutation comparison from local structures."""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.app.services.protein_service import ProteinService
from backend.app.services.ligand_service import LigandService
from backend.app.services.box_service import BoxService
from backend.app.services.docking_service import DockingService
from backend.app.services.validation_service import ValidationService
from backend.app.services.reliability_service import ReliabilityService

EXAMPLE = ROOT / "examples" / "1hsg"
OUT = ROOT / "docs" / "case_study"
VARIANTS = [
    {"name": "wild_type", "pdb": EXAMPLE / "1hsg.pdb", "resname": "MK1", "chain": "A", "chains": ["A", "B"]},
    {"name": "I50V_A71V", "pdb": EXAMPLE / "mutants" / "3oxv.pdb", "resname": "478", "chain": "A", "chains": ["A", "B"]},
    {"name": "I84V_L90M", "pdb": EXAMPLE / "mutants" / "2i4v.pdb", "resname": "DJR", "chain": "B", "chains": ["A", "B"]},
    {"name": "I50V", "pdb": EXAMPLE / "mutants" / "2r43.pdb", "resname": "G3G", "chain": "B", "chains": ["A", "B"]},
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--skip-reliability", action="store_true", help="Run only the single-seed comparison.")
    parser.add_argument("--limit-ligands", type=int, default=None, help="Optional short diagnostic run using the first N inhibitors.")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    ligand_service, protein_service = LigandService(), ProteinService()
    box_service, docking, validation = BoxService(), DockingService(), ValidationService()
    reliability = ReliabilityService(docking)

    with (EXAMPLE / "actives.csv").open(encoding="utf-8", newline="") as f:
        ligands = list(csv.DictReader(f))
    if args.limit_ligands is not None:
        ligands = ligands[:max(1, args.limit_ligands)]
    ligand_preps = {}
    for ligand in ligands:
        prep_dir = OUT / "prepared_ligands"
        pdbqt, sdf = prep_dir / f"{ligand['name']}.pdbqt", prep_dir / f"{ligand['name']}.sdf"
        ligand_service.prepare_ligand_from_smiles(ligand["smiles"], pdbqt, sdf, seed=args.seed)
        ligand_preps[ligand["name"]] = (pdbqt, sdf)

    prepared = {}
    for variant in VARIANTS:
        variant_dir = OUT / "runs" / variant["name"]
        variant_dir.mkdir(parents=True, exist_ok=True)
        crystal = variant_dir / "crystal_ligand.pdb"
        protein_service.extract_crystal_ligand(variant["pdb"], variant["resname"], crystal, chain_id=variant["chain"])
        receptor_pdb = variant_dir / "receptor.pdb"
        receptor_pdbqt = variant_dir / "receptor.pdbqt"
        prep = protein_service.clean_and_prepare_receptor(variant["pdb"], receptor_pdb, receptor_pdbqt, keep_chains=variant["chains"])
        box = box_service.box_from_crystal_ligand(crystal)
        prepared[variant["name"]] = {"dir": variant_dir, "receptor_pdb": receptor_pdb, "receptor_pdbqt": receptor_pdbqt, "box": box, "warnings": prep["warnings"]}

    records = []
    by_pair = {}
    for variant in VARIANTS:
        p = prepared[variant["name"]]
        for ligand in ligands:
            name = ligand["name"]
            pdbqt, source_sdf = ligand_preps[name]
            run_dir = p["dir"] / name
            run_dir.mkdir(parents=True, exist_ok=True)
            pose, pose_sdf = run_dir / "pose.pdbqt", run_dir / "pose.sdf"
            try:
                result = docking.dock(p["receptor_pdbqt"], pdbqt, p["box"]["center"], p["box"]["size"], pose, pose_sdf,
                    exhaustiveness=args.exhaustiveness, seed=args.seed)
                score = float(result["best_score"])
                interactions = validation.interaction_summary(p["receptor_pdb"], pose_sdf)
                clashes = validation.steric_clash_count(p["receptor_pdb"], pose_sdf)
                reliability_result = None
                if not args.skip_reliability:
                    reliability_result = reliability.run_ligand(p["receptor_pdbqt"], pdbqt, p["box"]["center"], p["box"]["size"],
                        run_dir / "reliability", list(range(args.seed, args.seed + args.seeds)), args.exhaustiveness)
                row = {"variant": variant["name"], "ligand": name, "label": ligand.get("label", "active"), "score_kcal_mol": score,
                    "contacts": interactions["residues"], "hbond_interactions": interactions["hydrogen_bond_interactions"],
                    "steric_clashes": clashes, "reliability_label": reliability_result["label"] if reliability_result else "not_run",
                    "reliability_score_sd": reliability_result["score_sd"] if reliability_result else None,
                    "pose_consistency": reliability_result["pose_consistency"] if reliability_result else None,
                    "preparation_warnings": "; ".join(p["warnings"])}
                by_pair[(variant["name"], name)] = row
                records.append(row)
                print(f"{variant['name']:14s} {name:14s} {score:7.2f} kcal/mol | {row['reliability_label']}", flush=True)
            except Exception as exc:
                row = {"variant": variant["name"], "ligand": name, "label": ligand.get("label", "active"), "status": "failed", "error": str(exc)[:500]}
                records.append(row)
                by_pair[(variant["name"], name)] = row
                print(f"{variant['name']} {name}: FAILED: {exc}", flush=True)

    output_fields = sorted({key for row in records for key in row})
    with (OUT / "results.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=output_fields); writer.writeheader()
        for row in records:
            row = dict(row)
            if isinstance(row.get("contacts"), list): row["contacts"] = ";".join(row["contacts"])
            writer.writerow(row)

    contact_rows = []
    for ligand in ligands:
        wt = by_pair.get(("wild_type", ligand["name"]), {})
        wt_contacts = set(wt.get("contacts", []))
        for variant in VARIANTS[1:]:
            row = by_pair.get((variant["name"], ligand["name"]), {})
            contacts = set(row.get("contacts", []))
            contact_rows.append({"variant": variant["name"], "ligand": ligand["name"],
                "delta_score_vs_wt": (row["score_kcal_mol"] - wt["score_kcal_mol"]) if row.get("score_kcal_mol") is not None and wt.get("score_kcal_mol") is not None else None,
                "contacts_lost_vs_wt": ";".join(sorted(wt_contacts - contacts)), "contacts_gained_vs_wt": ";".join(sorted(contacts - wt_contacts)),
                "wt_contact_retention": len(wt_contacts & contacts) / len(wt_contacts) if wt_contacts else None,
                "hbond_delta_vs_wt": row.get("hbond_interactions", 0) - wt.get("hbond_interactions", 0) if row.get("hbond_interactions") is not None else None,
                "key_residue_recovery": "not_applicable: this is screening, not redocking of each ligand's crystal pose"})
    with (OUT / "contact_changes.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=contact_rows[0].keys()); writer.writeheader(); writer.writerows(contact_rows)

    plot_rows = [row for row in records if row.get("score_kcal_mol") is not None]
    if plot_rows:
        variants = [v["name"] for v in VARIANTS]
        names = [ligand["name"] for ligand in ligands]
        vals = [row["score_kcal_mol"] for row in plot_rows]
        lo, hi = min(vals) - .5, max(vals) + .5
        width, height = max(760, 120 + 65 * len(names)), 390
        xstep = (width - 145) / max(1, len(names) - 1)
        y = lambda value: 48 + (hi - value) / max(.01, hi - lo) * 250
        colors = ["#161616", "#a63a28", "#68645c", "#9a6b00"]
        svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Measured Vina score by HIV protease structure">', '<rect width="100%" height="100%" fill="#f3f0e8"/>', '<text x="24" y="25" font-family="Georgia" font-size="17">Vina score by receptor variant (kcal/mol)</text>']
        for i, variant in enumerate(variants):
            points = []
            for j, name in enumerate(names):
                row = by_pair.get((variant, name), {})
                if row.get("score_kcal_mol") is not None:
                    points.append((80 + j * xstep, y(row["score_kcal_mol"]), name))
            if points:
                svg_points = " ".join(f"{x:.1f},{yy:.1f}" for x, yy, _ in points)
                svg.append(f'<polyline fill="none" stroke="{colors[i]}" stroke-width="1.5" points="{svg_points}"/>')
                for x, yy, name in points:
                    svg.append(f'<circle cx="{x:.1f}" cy="{yy:.1f}" r="3" fill="{colors[i]}"/>')
            svg.append(f'<text x="{width-130}" y="{50+i*18}" font-family="monospace" font-size="10" fill="{colors[i]}">{variant}</text>')
        for i, name in enumerate(names):
            svg.append(f'<text transform="translate({80+i*xstep:.1f},320) rotate(-40)" font-family="monospace" font-size="9">{name}</text>')
        svg.append('</svg>')
        (OUT / "scores_by_variant.svg").write_text("\n".join(svg), encoding="utf-8")

    # Fold-change phenotypic values must be curated before a correlation is emitted.
    references = []
    with (OUT / "reference_data.csv").open(encoding="utf-8", newline="") as f:
        references = [r for r in csv.DictReader(f) if r.get("fold_change")]
    summary = {"structures": len(VARIANTS), "ligands": len(ligands), "rows": len(records), "diagnostic_ligand_limit": args.limit_ligands,
               "successful_docks": sum(r.get("status") != "failed" for r in records),
               "failed_docks": sum(r.get("status") == "failed" for r in records),
               "reference_fold_change_rows": len(references), "spearman_delta_score_vs_log_fold_change": None,
               "reference_correlation_status": "TODO: no numeric fold-change data were present in reference_data.csv" if not references else "pending numeric analysis"}
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["failed_docks"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
