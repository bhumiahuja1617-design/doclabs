"""Build a reproducible EGFR actives/weak-activity comparison from cached ChEMBL data."""
import csv
import json
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "egfr"
SOURCES = [OUT / name for name in ("chembl_ic50.json", "chembl_weak.json", "chembl_weak_2.json")]


def properties(smiles: str) -> tuple[float, float, int, int, int]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    return (Descriptors.MolWt(mol), Descriptors.MolLogP(mol), rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol), rdMolDescriptors.CalcNumRotatableBonds(mol))


def rows():
    data = [activity for source in SOURCES for activity in json.loads(source.read_text(encoding="utf-8"))["activities"]]
    unique: dict[str, dict] = {}
    for item in data:
        if (item.get("target_chembl_id") != "CHEMBL203" or item.get("standard_type") != "IC50"
                or item.get("standard_units") != "nM" or item.get("standard_relation") != "="
                or item.get("assay_type") != "B" or item.get("data_validity_comment")
                or not item.get("canonical_smiles") or not item.get("standard_value")):
            continue
        try:
            value = float(item["standard_value"])
            p = properties(item["canonical_smiles"])
        except (ValueError, TypeError):
            continue
        key = item["molecule_chembl_id"]
        if key not in unique or value < float(unique[key]["standard_value"]):
            unique[key] = {**item, "_value": value, "_props": p}
    actives = sorted((r for r in unique.values() if r["_value"] <= 1000), key=lambda r: (r["_value"], r["molecule_chembl_id"]))
    weak = [r for r in unique.values() if r["_value"] >= 10000]
    chosen_active, smiles_seen = [], set()
    for item in actives:
        mol = Chem.MolFromSmiles(item["canonical_smiles"])
        key = Chem.MolToSmiles(mol, isomericSmiles=True)
        if key not in smiles_seen:
            chosen_active.append(item)
            smiles_seen.add(key)
        if len(chosen_active) == 20:
            break
    if len(chosen_active) < 10:
        raise RuntimeError(f"Only {len(chosen_active)} unique active compounds met the threshold.")

    # Match weakly active measured molecules to the nearest active by normalized
    # RDKit descriptors (MW, logP, HBD, HBA, rotatable bonds); do not invent SMILES.
    active_props = [r["_props"] for r in chosen_active]
    scales = [max(1.0, max(p[i] for p in active_props) - min(p[i] for p in active_props)) for i in range(5)]
    scored_weak, weak_smiles = [], set()
    for item in weak:
        mol = Chem.MolFromSmiles(item["canonical_smiles"])
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        if canonical in smiles_seen or canonical in weak_smiles:
            continue
        nearest = min(chosen_active, key=lambda active: sum(((item["_props"][i] - active["_props"][i]) / scales[i]) ** 2 for i in range(5)))
        distance = sum(((item["_props"][i] - nearest["_props"][i]) / scales[i]) ** 2 for i in range(5))
        item["_nearest_active"] = nearest
        item["_match_distance"] = distance
        scored_weak.append(item)
        weak_smiles.add(canonical)
    chosen_weak = sorted(scored_weak, key=lambda r: (r["_match_distance"], r["_value"], r["molecule_chembl_id"]))[:100]
    return chosen_active, chosen_weak


def write(path: Path, items: list[dict], label: str):
    fields = ["name", "chembl_id", "smiles", "label", "standard_type", "standard_relation", "standard_value_nM", "assay_chembl_id", "document_chembl_id", "document_year", "source_url", "assay_url", "nearest_active_chembl_id"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, r in enumerate(items, 1):
            writer.writerow({"name": f"{label}_{i:02d}", "chembl_id": r["molecule_chembl_id"], "smiles": r["canonical_smiles"],
                "label": label, "standard_type": "IC50", "standard_relation": "=", "standard_value_nM": r["standard_value"],
                "assay_chembl_id": r["assay_chembl_id"], "document_chembl_id": r["document_chembl_id"], "document_year": r["document_year"],
                "source_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{r['molecule_chembl_id']}/",
                "assay_url": f"https://www.ebi.ac.uk/chembl/assay_report_card/{r['assay_chembl_id']}/",
                "nearest_active_chembl_id": r.get("_nearest_active", {}).get("molecule_chembl_id", "")})


if __name__ == "__main__":
    active, decoy = rows()
    write(OUT / "actives.csv", active, "active")
    write(OUT / "decoys.csv", decoy, "decoy")
    print(f"Wrote {len(active)} unique EGFR IC50 actives (<=1,000 nM) and {len(decoy)} measured weak compounds (>=10,000 nM).")
    print("Decoys are descriptor-matched (RDKit MW, logP, HBD, HBA, rotatable bonds); all labels retain assay/source URLs.")
