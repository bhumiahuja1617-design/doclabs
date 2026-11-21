"""
backend/app/services/box_service.py
Service for calculating 3D docking search grid coordinates and dimensions.
Supports crystal ligand centroid, residue-based pockets, manual input, and fpocket integration.
"""

import sys
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import numpy as np
from Bio.PDB import PDBParser

class BoxService:
    def __init__(self, bin_dir: Optional[Path] = None):
        self.root_dir = Path(__file__).resolve().parent.parent.parent.parent
        self.bin_dir = bin_dir or (self.root_dir / "bin")

    def box_from_crystal_ligand(
        self,
        crystal_pdb_or_pdbqt: Path,
        padding: float = 10.0,
        min_size: float = 15.0,
        max_size: float = 30.0
    ) -> Dict[str, Any]:
        """
        Computes bounding box centered on crystal ligand centroid, padded by 10 Angstroms,
        clamped to [15.0, 30.0] Angstroms per side.
        """
        coords = []
        with open(crystal_pdb_or_pdbqt, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue

        if not coords:
            raise ValueError(f"No atomic coordinates found in {crystal_pdb_or_pdbqt.name}")

        arr = np.array(coords)
        centroid = np.mean(arr, axis=0)
        min_pt = np.min(arr, axis=0)
        max_pt = np.max(arr, axis=0)
        extent = max_pt - min_pt
        
        # Dimensions with padding and clamping
        dims = np.clip(extent + padding, min_size, max_size)

        return {
            "center": {
                "x": round(float(centroid[0]), 3),
                "y": round(float(centroid[1]), 3),
                "z": round(float(centroid[2]), 3),
            },
            "size": {
                "x": round(float(dims[0]), 3),
                "y": round(float(dims[1]), 3),
                "z": round(float(dims[2]), 3),
            },
            "source": "crystal_ligand",
            "extent": [round(float(e), 3) for e in extent],
        }

    def box_from_residues(
        self,
        protein_pdb: Path,
        residue_ids: List[Tuple[str, int]],  # [(chain, resnum)]
        padding: float = 10.0,
        min_size: float = 15.0,
        max_size: float = 30.0
    ) -> Dict[str, Any]:
        """
        Computes bounding box from specified active site residues.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("rec", str(protein_pdb))
        
        target_res = set(residue_ids)
        coords = []

        for model in structure:
            for chain in model:
                for res in chain:
                    key = (chain.id, res.id[1])
                    if key in target_res:
                        for atom in res:
                            coords.append(atom.get_coord())

        if not coords:
            raise ValueError("None of the specified residues were found in the protein structure.")

        arr = np.array(coords)
        centroid = np.mean(arr, axis=0)
        min_pt = np.min(arr, axis=0)
        max_pt = np.max(arr, axis=0)
        extent = max_pt - min_pt
        dims = np.clip(extent + padding, min_size, max_size)

        return {
            "center": {
                "x": round(float(centroid[0]), 3),
                "y": round(float(centroid[1]), 3),
                "z": round(float(centroid[2]), 3),
            },
            "size": {
                "x": round(float(dims[0]), 3),
                "y": round(float(dims[1]), 3),
                "z": round(float(dims[2]), 3),
            },
            "source": "residues",
            "residues": [f"{r[0]}:{r[1]}" for r in residue_ids]
        }

    def validate_manual_box(
        self,
        center: Dict[str, float],
        size: Dict[str, float]
    ) -> Dict[str, Any]:
        """
        Validates manual box boundaries. Maximum 30 Angstroms per side, warns if > 25 Angstroms.
        """
        warnings = []
        for axis in ["x", "y", "z"]:
            dim = size.get(axis, 0.0)
            if dim <= 0.0:
                raise ValueError(f"Box dimension '{axis}' must be greater than 0.")
            if dim > 30.0:
                raise ValueError(f"Box dimension '{axis}' ({dim:.1f} A) exceeds maximum permitted limit of 30.0 A.")
            if dim > 25.0:
                warnings.append(f"Box dimension '{axis}' ({dim:.1f} A) is unusually large (>25 A), which may dilute sampling exhaustiveness.")

        return {
            "center": {k: round(float(v), 3) for k, v in center.items()},
            "size": {k: round(float(v), 3) for k, v in size.items()},
            "source": "manual",
            "warnings": warnings
        }

    def detect_pockets_fpocket(
        self,
        protein_pdb: Path
    ) -> Dict[str, Any]:
        """
        Executes fpocket if available on system. Degrades gracefully on Windows.
        """
        if sys.platform == "win32":
            return {
                "available": False,
                "message": "fpocket is not supported natively on Windows. Please use crystal ligand centroid or residue-defined bounding boxes.",
                "pockets": []
            }

        # On Linux/container:
        try:
            res = subprocess.run(
                ["fpocket", "-f", str(protein_pdb)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=60
            )
            if res.returncode != 0:
                return {"available": True, "error": res.stderr, "pockets": []}
            return {"available": True, "pockets": [], "message": "Pocket detection executed."}
        except FileNotFoundError:
            return {
                "available": False,
                "message": "fpocket binary not found on system PATH.",
                "pockets": []
            }

