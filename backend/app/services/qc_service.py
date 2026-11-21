"""Pre-docking input checks. Findings are diagnostics, never silent filters."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


@dataclass
class QCCheck:
    id: str
    title: str
    severity: str
    detail: str
    how_detected: str
    suggestion: str


class QCService:
    """Inspect receptor, ligand, and box files before docking."""

    def run(
        self,
        receptor_pdb: Path,
        ligands: list[dict[str, Any]],
        box: dict[str, Any],
        crystal_ligand_pdb: Path | None = None,
    ) -> dict[str, Any]:
        checks: list[QCCheck] = []
        structure = PDBParser(QUIET=True).get_structure("qc", str(receptor_pdb))
        atoms = [a for a in structure.get_atoms() if a.element not in {"H", "D"}]
        if not atoms:
            checks.append(QCCheck("protein_atoms", "Protein atoms", "fail", "No heavy atoms were parsed.", "Parsed ATOM records with Biopython.", "Supply a valid protein structure."))
            coords = np.empty((0, 3))
        else:
            coords = np.asarray([a.coord for a in atoms], dtype=float)
            checks.append(QCCheck("protein_atoms", "Protein atoms", "pass", f"Parsed {len(atoms)} heavy atoms.", "Parsed coordinates and element fields with Biopython.", ""))

        residues = [(c.id, r) for c in structure[0] for r in c]
        pdb_lines = receptor_pdb.read_text(errors="ignore").splitlines()
        gaps = sum(1 for line in pdb_lines if line.startswith("REMARK 465") and len(line) > 15 and line[15:18].strip())
        checks.append(QCCheck("missing_residues", "Missing residues", "warn" if gaps else "pass", f"PDB header lists {gaps} REMARK 465 residue records." if gaps else "No REMARK 465 missing-residue records were found.", "Counted REMARK 465 records in the supplied PDB header; this does not infer absent sequence positions.", "Review missing loops near the binding site." if gaps else ""))
        missing_atom_rows = [line for line in pdb_lines if line.startswith("REMARK 470") and len(line) > 15 and line[15:18].strip()]
        checks.append(QCCheck("missing_atoms", "Missing heavy atoms", "warn" if missing_atom_rows else "pass", f"PDB header lists {len(missing_atom_rows)} REMARK 470 atom records." if missing_atom_rows else "No REMARK 470 missing-atom records were found.", "Counted REMARK 470 records; proximity to the search box is not inferred here.", "Review missing atoms in residues near the pocket." if missing_atom_rows else ""))
        chains = sorted({c.id for c in structure[0]})
        biological_assembly = any(line.startswith("REMARK 350") for line in pdb_lines)
        checks.append(QCCheck("chains_assembly", "Chains and biological assembly", "warn" if len(chains) > 1 and not biological_assembly else "pass", f"Parsed chains: {', '.join(chains) or 'none'}; assembly remarks {'present' if biological_assembly else 'absent'}.", "Read Biopython chain IDs and scanned REMARK 350; this does not predict removed-chain contacts.", "Confirm the selected chains represent the biological assembly." if len(chains) > 1 or not biological_assembly else ""))

        waters = [r for _, r in residues if r.resname in {"HOH", "WAT"}]
        lig_atoms = []
        if crystal_ligand_pdb and crystal_ligand_pdb.exists():
            lp = PDBParser(QUIET=True).get_structure("lig", str(crystal_ligand_pdb))
            lig_atoms = [a for a in lp.get_atoms() if a.element not in {"H", "D"}]
        bridging = 0
        if lig_atoms and atoms:
            protein_xyz = np.asarray([a.coord for a in atoms])
            for water in waters:
                wat_xyz = np.asarray([a.coord for a in water if a.element not in {"H", "D"}])
                if not len(wat_xyz):
                    continue
                d_lig = np.linalg.norm(np.asarray([a.coord for a in lig_atoms])[:, None, :] - wat_xyz[None, :, :], axis=2)
                d_prot = np.linalg.norm(protein_xyz[:, None, :] - wat_xyz[None, :, :], axis=2)
                if float(d_lig.min()) <= 3.5 and int((d_prot.min(axis=1) <= 3.5).sum()) >= 2:
                    bridging += 1
        checks.append(QCCheck("waters", "Potential bridging waters", "warn" if bridging else "pass", f"{bridging} waters meet the geometric bridging rule." if lig_atoms else "Crystal ligand unavailable; bridging-water geometry was not evaluated.", "Checked water-to-ligand and water-to-protein heavy-atom distances at 3.5 Å.", "Inspect these waters before removing them." if bridging else ""))

        het = sorted({r.resname for _, r in residues if r.id[0].startswith("H_") and r.resname not in {"HOH", "WAT", "SO4", "GOL", "EDO", "PEG", "ACT", "CL", "NA"}})
        checks.append(QCCheck("heterogens", "Cofactors and other heterogens", "warn" if het else "pass", "Nearby heterogen residue names: " + ", ".join(het) if het else "No non-additive heterogens were detected.", "Inspected parsed HETATM residues; proximity is not assessed in this version.", "Confirm whether any listed cofactor belongs in the receptor." if het else ""))

        center, size = box.get("center", {}), box.get("size", {})
        dims = [float(size.get(a, 0)) for a in "xyz"]
        valid = all(np.isfinite(dims)) and all(0 < d <= 30 for d in dims) and all(np.isfinite(float(center.get(a, 0))) for a in "xyz")
        checks.append(QCCheck("box_bounds", "Search box dimensions", "fail" if not valid else ("warn" if max(dims) > 25 else "pass"), f"Box size: {dims[0]:g} × {dims[1]:g} × {dims[2]:g} Å.", "Validated finite center and positive dimensions no larger than 30 Å per side.", "Set a finite center and dimensions from 0 to 30 Å." if not valid else ("Consider whether a box above 25 Å is adequately sampled." if max(dims) > 25 else "")))
        if valid and len(coords):
            low = np.array([center[a] for a in "xyz"], dtype=float) - np.array(dims) / 2
            high = np.array([center[a] for a in "xyz"], dtype=float) + np.array(dims) / 2
            inside = np.all((coords >= low) & (coords <= high), axis=1).sum()
            severity = "pass" if inside else "fail"
            checks.append(QCCheck("box_protein_overlap", "Protein overlaps box", severity, f"{inside} receptor heavy atoms fall inside the box.", "Counted receptor heavy-atom coordinates inside the axis-aligned box.", "Recenter or enlarge the box to include the binding site." if not inside else ""))
        if valid and lig_atoms:
            c = np.asarray([center[a] for a in "xyz"], dtype=float)
            half = np.asarray(dims) / 2
            contained = sum(bool(np.all(np.abs(a.coord - c) <= half)) for a in lig_atoms)
            checks.append(QCCheck("box_contains_crystal", "Crystal ligand inside box", "pass" if contained == len(lig_atoms) else "warn", f"{contained} of {len(lig_atoms)} crystal-ligand heavy atoms are inside the box.", "Checked crystal-ligand coordinates against box faces.", "Expand or recenter the box to contain the complete reference ligand." if contained != len(lig_atoms) else ""))

        if ligands:
            unassigned: list[str] = []
            invalid: list[str] = []
            sanity: list[str] = []
            for item in ligands:
                mol = Chem.MolFromSmiles(str(item.get("smiles", "")))
                name = str(item.get("name", "unnamed"))
                if mol is None:
                    invalid.append(name)
                    continue
                if any(label == "?" for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True)):
                    unassigned.append(name)
                if Descriptors.MolWt(mol) > 700 or rdMolDescriptors.CalcNumRotatableBonds(mol) > 20 or any(atom.GetSymbol() not in {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"} for atom in mol.GetAtoms()):
                    sanity.append(name)
            sev = "fail" if invalid else ("warn" if unassigned else "pass")
            detail = f"Invalid SMILES: {', '.join(invalid)}." if invalid else (f"Unassigned stereocenters in: {', '.join(unassigned)}." if unassigned else f"All {len(ligands)} ligand SMILES parsed and assigned stereochemistry.")
            checks.append(QCCheck("ligand_inputs", "Ligand validity and stereochemistry", sev, detail, "Parsed SMILES with RDKit and checked FindMolChiralCenters(includeUnassigned=True).", "Correct invalid SMILES and specify stereochemistry where known." if invalid or unassigned else ""))
            checks.append(QCCheck("ligand_sanity", "Ligand size, flexibility, and elements", "warn" if sanity else "pass", f"Above MW/flexibility limits or contains a nonstandard element: {', '.join(sanity)}." if sanity else "No ligand exceeds 700 Da or 20 rotatable bonds; no unhandled elements found.", "Calculated RDKit molecular weight, rotatable bonds, and atom symbols.", "Review these molecules; drug-likeness metrics are informational, not filters." if sanity else ""))
        else:
            checks.append(QCCheck("ligand_inputs", "Ligand inputs", "fail", "No ligands were supplied.", "Checked ligand input list length.", "Add at least one ligand."))

        if valid and len(coords):
            c = np.asarray([center[a] for a in "xyz"], dtype=float)
            near = [res for _, res in residues if res.resname in {"HIS", "ASP", "GLU", "LYS"} and any(np.linalg.norm(a.coord - c) <= 8 for a in res)]
            names = sorted({r.resname for r in near})
            checks.append(QCCheck("protonation_note", "Pocket protonation note", "warn" if near else "pass", f"{len(near)} potentially ionizable residues lie within 8 Å of the box center ({', '.join(names)})." if near else "No HIS, ASP, GLU, or LYS residues were found within 8 Å of the box center.", "Counted selected ionizable residue types using Biopython atom-to-center distances; pKa values are not calculated.", "Review histidine tautomers and catalytic carboxylates if they affect the binding site." if near else ""))

        counts = {s: sum(c.severity == s for c in checks) for s in ("pass", "warn", "fail")}
        return {"checks": [asdict(c) for c in checks], "counts": counts, "summary": f"QC: {counts['pass']} passed, {counts['warn']} warnings, {counts['fail']} failures", "requires_acknowledgement": counts["fail"] > 0}
