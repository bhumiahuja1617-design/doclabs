"""
backend/app/services/protein_service.py
Pure Python service for protein ingestion, cleaning, ligand extraction, and PDBQT conversion.
"""

import sys
import os
import shutil
import urllib.request
import subprocess
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

from Bio.PDB import PDBParser, PDBIO, Select
import pdbfixer
from openmm.app import PDBFile

# Standard crystallization additives and solvents to ignore
COMMON_ADDITIVES = {
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "CL", "NA",
    "MG", "ZN", "CA", "K", "DMS", "TRS", "BME", "MPD", "IPA", "FMT"
}

class ProteinService:
    def __init__(self, bin_dir: Optional[Path] = None):
        self.root_dir = Path(__file__).resolve().parent.parent.parent.parent
        self.bin_dir = bin_dir or (self.root_dir / "bin")

    def _get_obabel_cmd(self) -> str:
        # Check local bin, conda Library/bin, or system PATH
        candidates = [
            self.bin_dir / ("obabel.exe" if sys.platform == "win32" else "obabel"),
            Path(sys.executable).parent / "Library" / "bin" / ("obabel.exe" if sys.platform == "win32" else "obabel"),
            Path(sys.executable).parent / "bin" / ("obabel.exe" if sys.platform == "win32" else "obabel"),
        ]
        for c in candidates:
            if c.exists():
                return str(c)
        return "obabel"

    def fetch_or_load_pdb(self, pdb_id_or_path: str, cache_dir: Path) -> Path:
        cache_dir.mkdir(parents=True, exist_ok=True)
        raw_path = Path(pdb_id_or_path)
        if raw_path.exists() and raw_path.is_file():
            dest = cache_dir / raw_path.name
            if dest != raw_path:
                shutil.copy(raw_path, dest)
            return dest
        
        # Treat as PDB ID
        pdb_id = pdb_id_or_path.strip().upper()
        dest_file = cache_dir / f"{pdb_id}.pdb"
        # A previous interrupted download must not poison every later fetch.
        # Only reuse cached files that contain actual coordinate records.
        if dest_file.exists() and self._has_pdb_coordinates(dest_file):
            return dest_file
        dest_file.unlink(missing_ok=True)

        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (DockLab/1.0)"})
        temporary = cache_dir / f".{pdb_id}.download"
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                content = resp.read()
            temporary.write_bytes(content)
            if not self._has_pdb_coordinates(temporary):
                raise ValueError(f"RCSB returned no PDB coordinates for {pdb_id}. Check the PDB ID or try Search public proteins.")
            os.replace(temporary, dest_file)
        except Exception as exc:
            temporary.unlink(missing_ok=True)
            if isinstance(exc, ValueError):
                raise
            raise RuntimeError(f"Could not download PDB {pdb_id} from RCSB ({exc.__class__.__name__}). Check your internet connection or use Search public proteins.") from exc
        return dest_file

    @staticmethod
    def _has_pdb_coordinates(path: Path) -> bool:
        try:
            has_coordinates = False
            has_end_record = False
            with path.open("r", encoding="ascii", errors="ignore") as pdb:
                for line in pdb:
                    has_coordinates = has_coordinates or line.startswith(("ATOM  ", "HETATM"))
                    has_end_record = has_end_record or line.strip() in {"END", "ENDMDL"}
            return has_coordinates and has_end_record
        except OSError:
            return False

    def inspect_structure(self, pdb_file: Path) -> Dict[str, Any]:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", str(pdb_file))
        
        chains = []
        hetatms = []
        waters_count = 0
        total_atoms = 0

        for model in structure:
            for chain in model:
                chains.append(chain.id)
                for res in chain:
                    res_name = res.resname.strip()
                    atom_count = len(res)
                    total_atoms += atom_count
                    if res.id[0] == "W" or res_name in {"HOH", "WAT"}:
                        waters_count += 1
                    elif res.id[0].startswith("H_"):
                        if res_name not in COMMON_ADDITIVES:
                            hetatms.append({
                                "resname": res_name,
                                "chain": chain.id,
                                "res_num": res.id[1],
                                "atom_count": atom_count
                            })

        if total_atoms == 0:
            raise ValueError(f"Protein file {pdb_file.name} contains zero atoms.")

        return {
            "chains": sorted(list(set(chains))),
            "ligands_detected": hetatms,
            "waters_count": waters_count,
            "total_atoms": total_atoms
        }

    def extract_crystal_ligand(
        self,
        raw_pdb_file: Path,
        ligand_resname: str,
        output_pdb: Path,
        chain_id: Optional[str] = None,
    ) -> Path:
        """
        Extract the co-crystallized ligand directly from original experimental crystal coordinates.
        Preserves original coordinate space and CONECT records without converting through SMILES.
        """
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        hetatm_lines = []
        lig_serials = set()

        with open(raw_pdb_file, "r") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand_resname and (chain_id is None or line[21:22].strip() == chain_id):
                    hetatm_lines.append(line)
                    try:
                        lig_serials.add(int(line[6:11].strip()))
                    except ValueError:
                        pass

        if not hetatm_lines:
            raise ValueError(f"No ligand with residue name '{ligand_resname}' found in {raw_pdb_file.name}.")

        conect_lines = []
        with open(raw_pdb_file, "r") as f:
            for line in f:
                if line.startswith("CONECT"):
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            src_serial = int(parts[1])
                            if src_serial in lig_serials:
                                conect_lines.append(line)
                        except ValueError:
                            pass

        with open(output_pdb, "w") as f:
            f.writelines(hetatm_lines)
            f.writelines(conect_lines)

        return output_pdb

    def clean_and_prepare_receptor(
        self,
        raw_pdb_file: Path,
        cleaned_pdb_file: Path,
        pdbqt_file: Path,
        keep_chains: Optional[List[str]] = None,
        pH: float = 7.4
    ) -> Dict[str, Any]:
        """
        Cleans the receptor with PDBFixer: strips heterogens/waters, fixes missing heavy atoms,
        adds hydrogens at pH 7.4, and converts to PDBQT via Open Babel with Meeko fallback.
        """
        cleaned_pdb_file.parent.mkdir(parents=True, exist_ok=True)
        pdbqt_file.parent.mkdir(parents=True, exist_ok=True)

        warnings = []
        fixer_input = raw_pdb_file
        if keep_chains:
            selected_path = cleaned_pdb_file.with_name(cleaned_pdb_file.stem + "_selected_chains.pdb")
            structure = PDBParser(QUIET=True).get_structure("selected", str(raw_pdb_file))
            chain_ids = set(keep_chains)

            class ChainSelect(Select):
                def accept_chain(self, chain):
                    return chain.id in chain_ids

            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(str(selected_path), ChainSelect())
            fixer_input = selected_path
        fixer = pdbfixer.PDBFixer(filename=str(fixer_input))
        
        # Remove heterogens and crystallization waters
        fixer.removeHeterogens(keepWater=False)

        fixer.findMissingResidues()
        missing_res = getattr(fixer, "missingResidues", {})
        if missing_res:
            warnings.append(f"Missing residues detected: {len(missing_res)} gap regions.")

        fixer.findMissingAtoms()
        missing_atoms = getattr(fixer, "missingAtoms", {})
        if missing_atoms:
            warnings.append(f"Missing heavy atoms detected in {len(missing_atoms)} residues. Restored via PDBFixer.")
            fixer.addMissingAtoms()

        fixer.addMissingHydrogens(pH)

        with open(cleaned_pdb_file, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        # Convert to PDBQT using Open Babel with -xr (rigid receptor)
        obabel_cmd = self._get_obabel_cmd()
        cmd = [
            obabel_cmd,
            "-ipdb", str(cleaned_pdb_file),
            "-opdbqt", "-O", str(pdbqt_file),
            "-xr"
        ]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Failed to convert receptor to PDBQT with Open Babel: {res.stderr}")

        return {
            "cleaned_pdb": str(cleaned_pdb_file),
            "pdbqt": str(pdbqt_file),
            "warnings": warnings,
            "missing_residues_count": len(missing_res),
            "missing_atoms_count": len(missing_atoms)
        }

