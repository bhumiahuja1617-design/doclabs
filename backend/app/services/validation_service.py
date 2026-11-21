"""
backend/app/services/validation_service.py
Service for redocking validation, symmetry-aware heavy-atom RMSD,
and ProLIF chemical key-residue contact recovery.
"""

from pathlib import Path
from typing import Dict, Any, List, Set, Tuple
import numpy as np
from spyrmsd import rmsd as spyrmsd_rmsd
from spyrmsd.molecule import Molecule as SpyMolecule

import prolif as plf
from rdkit import Chem
import MDAnalysis as mda

class ValidationService:
    def __init__(self):
        pass

    def extract_heavy_coords_pdbqt(self, pdbqt_file: Path, first_model_only: bool = True) -> np.ndarray:
        coords = []
        with open(pdbqt_file, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atom_name = line[12:16].strip()
                    atom_type = line[77:].strip()
                    # Skip hydrogens
                    if atom_type != "H" and not atom_name.startswith("H"):
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                        except (ValueError, IndexError):
                            continue
                elif line.startswith("ENDMDL") and first_model_only:
                    break
        return np.array(coords)

    def compute_rmsd(self, crystal_pdbqt: Path, docked_pdbqt: Path) -> Dict[str, Any]:
        """
        Computes direct and symmetry-aware heavy-atom RMSD.
        Status: < 2.0 A (PASS), 2.0 - 3.0 A (CAUTION), > 3.0 A (FAIL).
        """
        c_coords = self.extract_heavy_coords_pdbqt(crystal_pdbqt, first_model_only=False)
        d_coords = self.extract_heavy_coords_pdbqt(docked_pdbqt, first_model_only=True)

        if len(c_coords) == 0 or len(d_coords) == 0:
            raise ValueError("No heavy atoms found in crystal or docked PDBQT.")

        if len(c_coords) != len(d_coords):
            raise ValueError(f"Heavy-atom count differs ({len(c_coords)} crystal, {len(d_coords)} docked); cannot compare atom mappings.")

        diff = c_coords - d_coords
        direct_rmsd = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

        # Determine qualitative status
        if direct_rmsd < 2.0:
            status = "PASS"
            advice = "Excellent pose reproduction. The docked conformation matches the crystal binding mode."
        elif direct_rmsd <= 3.0:
            status = "CAUTION"
            advice = "Moderate pose deviation. Consider expanding search box padding or increasing exhaustiveness."
        else:
            status = "FAIL"
            advice = "Redocking failed to reproduce crystal pose (>3.0 A). Check search box center, dimensions, and protonation."

        return {
            "rmsd_angstroms": round(direct_rmsd, 3),
            "heavy_atom_count": len(c_coords),
            "status": status,
            "advice": advice,
            "passed": direct_rmsd < 2.0
        }

    def compute_symmetric_rmsd(self, crystal_sdf: Path, docked_sdf: Path) -> Dict[str, Any]:
        """Calculate graph-isomorphism-aware RMSD without fitting away receptor-frame movement."""
        ref = Chem.SDMolSupplier(str(crystal_sdf), removeHs=True)[0]
        dock = Chem.SDMolSupplier(str(docked_sdf), removeHs=True)[0]
        if ref is None or dock is None:
            raise ValueError("Could not parse crystal or docked SDF pose.")
        if ref.GetNumHeavyAtoms() != dock.GetNumHeavyAtoms():
            raise ValueError("Crystal and docked ligands have different heavy-atom counts.")
        ref_spy, dock_spy = SpyMolecule.from_rdkit(ref), SpyMolecule.from_rdkit(dock)
        value = float(spyrmsd_rmsd.symmrmsd(
            ref_spy.coordinates, dock_spy.coordinates,
            ref_spy.atomicnums, dock_spy.atomicnums,
            ref_spy.adjacency_matrix, dock_spy.adjacency_matrix,
            center=False, minimize=False,
        ))
        if value < 2.0:
            status, advice = "PASS", "The heavy-atom pose is within 2 Å of the crystal pose."
        elif value <= 3.0:
            status, advice = "CAUTION", "The pose is moderately displaced; inspect the box and protonation."
        else:
            status, advice = "FAIL", "The pose differs by more than 3 Å; inspect receptor preparation and search settings."
        return {"rmsd_angstroms": round(value, 3), "heavy_atom_count": ref.GetNumHeavyAtoms(),
                "status": status, "advice": advice, "passed": value < 2.0, "method": "spyrmsd graph isomorphism; no coordinate fitting"}

    def compute_key_residue_recovery(
        self,
        receptor_pdb: Path,
        crystal_sdf: Path,
        docked_sdf: Path
    ) -> Dict[str, Any]:
        """
        Evaluates ProLIF non-covalent interaction conservation between crystal and docked pose.
        Calculates Key-Residue Recovery: |Crystal Interacting & Docked Interacting| / |Crystal Interacting|
        """
        u_rec = mda.Universe(str(receptor_pdb))
        prot = plf.Molecule.from_mda(u_rec)

        suppl_cryst = Chem.SDMolSupplier(str(crystal_sdf), removeHs=False)
        suppl_dock = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)

        if not suppl_cryst or suppl_cryst[0] is None:
            raise ValueError(f"Failed to read crystal ligand SDF: {crystal_sdf.name}")
        if not suppl_dock or suppl_dock[0] is None:
            raise ValueError(f"Failed to read docked ligand SDF: {docked_sdf.name}")

        mol_cryst = plf.Molecule.from_rdkit(suppl_cryst[0])
        mol_dock = plf.Molecule.from_rdkit(suppl_dock[0])

        fp = plf.Fingerprint()
        res_c = fp.generate(mol_cryst, prot)
        res_d = fp.generate(mol_dock, prot)

        # Map residues with interactions
        cryst_set = self._residue_set(res_c)
        docked_set = self._residue_set(res_d)

        common_set = cryst_set.intersection(docked_set)
        recovery = len(common_set) / len(cryst_set) if cryst_set else 0.0

        all_residues = sorted(list(cryst_set.union(docked_set)))
        residue_strip = []
        for r in all_residues:
            in_c = r in cryst_set
            in_d = r in docked_set
            state = "both" if (in_c and in_d) else ("crystal_only" if in_c else "docked_only")
            residue_strip.append({
                "residue": r,
                "in_crystal": in_c,
                "in_docked": in_d,
                "status": state
            })

        return {
            "recovery_fraction": round(recovery, 3),
            "recovery_percentage": round(recovery * 100.0, 1),
            "crystal_residues_count": len(cryst_set),
            "docked_residues_count": len(docked_set),
            "common_residues_count": len(common_set),
            "crystal_residues": sorted(list(cryst_set)),
            "docked_residues": sorted(list(docked_set)),
            "common_residues": sorted(list(common_set)),
            "residue_strip": residue_strip
        }

    @staticmethod
    def _residue_set(fingerprint: dict) -> set[str]:
        return {f"{key[1].name}{key[1].number}.{key[1].chain}" for key, values in fingerprint.items() if any(values)}

    def interaction_summary(self, receptor_pdb: Path, ligand_sdf: Path) -> Dict[str, Any]:
        """Return residue-level ProLIF contacts and donor/acceptor interaction counts."""
        universe = mda.Universe(str(receptor_pdb))
        protein = plf.Molecule.from_mda(universe)
        supplier = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)
        if not supplier or supplier[0] is None:
            raise ValueError(f"Failed to read ligand SDF: {ligand_sdf.name}")
        ligand = plf.Molecule.from_rdkit(supplier[0])
        fingerprint = plf.Fingerprint()
        result = fingerprint.generate(ligand, protein)
        hb_count = 0
        for values in result.values():
            arr = np.asarray(values, dtype=bool)
            for i, name in enumerate(fingerprint.interactions):
                if name in {"HBDonor", "HBAcceptor"} and i < arr.size:
                    hb_count += int(arr[i])
        return {"residues": sorted(self._residue_set(result)), "hydrogen_bond_interactions": hb_count,
                "interaction_types": list(fingerprint.interactions)}

    def steric_clash_count(self, receptor_pdb: Path, ligand_sdf: Path, tolerance: float = 0.4) -> int:
        """Count heavy-atom distances shorter than summed van der Waals radii minus tolerance."""
        from Bio.PDB import PDBParser
        ligand = Chem.SDMolSupplier(str(ligand_sdf), removeHs=True)[0]
        if ligand is None or ligand.GetNumConformers() == 0:
            raise ValueError(f"A 3D ligand SDF is required: {ligand_sdf.name}")
        receptor = PDBParser(QUIET=True).get_structure("clash", str(receptor_pdb))
        protein_atoms = [a for a in receptor.get_atoms() if a.element not in {"H", "D"}]
        if not protein_atoms:
            return 0
        rec_xyz = np.asarray([a.coord for a in protein_atoms])
        radii = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}
        conf = ligand.GetConformer()
        count = 0
        for atom in ligand.GetAtoms():
            if atom.GetSymbol() in {"H", "D"}:
                continue
            point = np.asarray(conf.GetAtomPosition(atom.GetIdx()))
            d = np.linalg.norm(rec_xyz - point, axis=1)
            cutoff = radii.get(atom.GetSymbol(), 1.7) + np.asarray([radii.get(a.element, 1.7) for a in protein_atoms]) - tolerance
            count += int(np.sum(d < cutoff))
        return count

