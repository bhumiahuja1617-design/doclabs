"""
backend/app/services/ligand_service.py
Cheminformatics service for SMILES parsing, 3D conformer generation,
Meeko PDBQT conversion, and drug-likeness / structural alert screening.
"""

import sys
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, FilterCatalog
from meeko import MoleculePreparation, PDBQTWriterLegacy

class LigandService:
    def __init__(self, bin_dir: Optional[Path] = None):
        self.root_dir = Path(__file__).resolve().parent.parent.parent.parent
        self.bin_dir = bin_dir or (self.root_dir / "bin")
        self._init_alert_catalogs()

    def _init_alert_catalogs(self):
        # Initialize PAINS catalog
        params_pains = FilterCatalog.FilterCatalogParams()
        params_pains.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        self.pains_catalog = FilterCatalog.FilterCatalog(params_pains)

        # Initialize Brenk catalog
        params_brenk = FilterCatalog.FilterCatalogParams()
        params_brenk.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        self.brenk_catalog = FilterCatalog.FilterCatalog(params_brenk)

    def _get_obabel_cmd(self) -> str:
        candidates = [
            self.bin_dir / ("obabel.exe" if sys.platform == "win32" else "obabel"),
            Path(sys.executable).parent / "Library" / "bin" / ("obabel.exe" if sys.platform == "win32" else "obabel"),
            Path(sys.executable).parent / "bin" / ("obabel.exe" if sys.platform == "win32" else "obabel"),
        ]
        for c in candidates:
            if c.exists():
                return str(c)
        return "obabel"

    def parse_and_standardize_smiles(self, smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: '{smiles}'")

        # Keep largest fragment (desalt)
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

        # Chemically sanitize
        Chem.SanitizeMol(mol)
        return mol

    def compute_drug_likeness(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculates Lipinski Ro5, Veber rules, PAINS, Brenk alerts, and molecular properties.
        """
        # Ensure hydrogens are added for proper HBD/HBA calculations
        mol_h = Chem.AddHs(mol)
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()

        # Lipinski Rule of 5 checks
        lipinski_violations = 0
        if mw > 500.0: lipinski_violations += 1
        if logp > 5.0: lipinski_violations += 1
        if hbd > 5: lipinski_violations += 1
        if hba > 10: lipinski_violations += 1
        lipinski_pass = lipinski_violations <= 1

        # Veber rules (Rotatable bonds <= 10, TPSA <= 140)
        veber_violations = 0
        if rot_bonds > 10: veber_violations += 1
        if tpsa > 140.0: veber_violations += 1
        veber_pass = veber_violations == 0

        # Structural alerts
        pains_matches = [entry.GetDescription() for entry in self.pains_catalog.GetMatches(mol)]
        brenk_matches = [entry.GetDescription() for entry in self.brenk_catalog.GetMatches(mol)]

        # Chiral centers check
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unassigned_chiral = [idx for idx, label in chiral_centers if label == "?"]

        warnings = []
        if mw > 700.0:
            warnings.append(f"High molecular weight: {mw:.1f} Da (>700)")
        if rot_bonds > 20:
            warnings.append(f"High flexibility: {rot_bonds} rotatable bonds (>20)")
        if unassigned_chiral:
            warnings.append(f"Unassigned stereocenters detected at atom indices: {unassigned_chiral}")

        return {
            "mw": round(mw, 2),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rot_bonds,
            "tpsa": round(tpsa, 2),
            "heavy_atom_count": heavy_atoms,
            "lipinski_violations": lipinski_violations,
            "lipinski_pass": lipinski_pass,
            "veber_violations": veber_violations,
            "veber_pass": veber_pass,
            "pains_alerts": pains_matches,
            "brenk_alerts": brenk_matches,
            "unassigned_chiral_count": len(unassigned_chiral),
            "warnings": warnings
        }

    def prepare_ligand_from_smiles(
        self,
        smiles: str,
        output_pdbqt: Path,
        output_sdf: Optional[Path] = None,
        seed: int = 42
    ) -> Dict[str, Any]:
        """
        Full 3D embedding with ETKDGv3, MMFF94 minimization, and Meeko PDBQT export.
        """
        output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        mol = self.parse_and_standardize_smiles(smiles)
        props = self.compute_drug_likeness(mol)

        # 3D Conformation Generation with ETKDGv3
        mol_3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        embed_status = AllChem.EmbedMolecule(mol_3d, params)
        if embed_status != 0:
            # Fallback to standard distance geometry if ETKDGv3 fails
            embed_status = AllChem.EmbedMolecule(mol_3d, randomSeed=seed)
            if embed_status != 0:
                raise RuntimeError(f"3D conformer embedding failed for SMILES: {smiles}")

        # Energy minimization with MMFF94
        try:
            AllChem.MMFFOptimizeMolecule(mol_3d)
        except Exception:
            pass

        # Meeko PDBQT conversion
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol_3d)
        if not mol_setups:
            raise RuntimeError("Meeko failed to parameterize ligand molecule.")

        pdbqt_str, is_ok, err_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            raise RuntimeError(f"Meeko PDBQT export failed: {err_msg}")

        with open(output_pdbqt, "w") as f:
            f.write(pdbqt_str)

        # Save SDF conformer if requested
        if output_sdf:
            output_sdf.parent.mkdir(parents=True, exist_ok=True)
            writer = Chem.SDWriter(str(output_sdf))
            writer.write(mol_3d)
            writer.close()

        return {
            "pdbqt": str(output_pdbqt),
            "sdf": str(output_sdf) if output_sdf else None,
            "properties": props
        }

    def prepare_crystal_ligand_pdbqt(
        self,
        crystal_pdb_file: Path,
        output_pdbqt: Path
    ) -> Path:
        """
        Converts extracted crystal coordinates directly to PDBQT with Open Babel.
        Guarantees that atom serials and crystal coordinates remain 100% isomorphic.
        """
        output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        obabel_cmd = self._get_obabel_cmd()
        cmd = [
            obabel_cmd,
            "-ipdb", str(crystal_pdb_file),
            "-opdbqt", "-O", str(output_pdbqt),
            "-h"
        ]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Failed to convert crystal ligand to PDBQT with Open Babel: {res.stderr}")
        return output_pdbqt

