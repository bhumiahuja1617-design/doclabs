"""
backend/app/services/docking_service.py
Pure Python service for AutoDock Vina molecular docking simulations.
Handles presets, subprocess execution with timeouts, pose parsing, and SDF conversions.
"""

import sys
import subprocess
import re
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple

class DockingService:
    PRESETS = {
        "fast": {"exhaustiveness": 4, "num_modes": 9, "energy_range": 3.0},
        "balanced": {"exhaustiveness": 8, "num_modes": 9, "energy_range": 3.0},
        "thorough": {"exhaustiveness": 32, "num_modes": 9, "energy_range": 3.0},
    }

    def __init__(self, bin_dir: Optional[Path] = None):
        self.root_dir = Path(__file__).resolve().parent.parent.parent.parent
        self.bin_dir = bin_dir or (self.root_dir / "bin")

    def _get_vina_cmd(self) -> str:
        candidates = [
            self.bin_dir / ("vina.exe" if sys.platform == "win32" else "vina"),
            Path(sys.executable).parent / ("vina.exe" if sys.platform == "win32" else "vina"),
            Path(sys.executable).parent / "Scripts" / ("vina.exe" if sys.platform == "win32" else "vina"),
        ]
        for c in candidates:
            if c.exists():
                return str(c)
        return "vina"

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

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box_center: Dict[str, float],
        box_size: Dict[str, float],
        output_pdbqt: Path,
        output_sdf: Optional[Path] = None,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        seed: int = 42,
        cpu: int = 0,
        timeout_seconds: int = 600
    ) -> Dict[str, Any]:
        """
        Runs AutoDock Vina simulation. Parses affinity scores and exports SDF poses.
        """
        output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        vina_cmd = self._get_vina_cmd()

        cmd = [
            vina_cmd,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(box_center["x"]),
            "--center_y", str(box_center["y"]),
            "--center_z", str(box_center["z"]),
            "--size_x", str(box_size["x"]),
            "--size_y", str(box_size["y"]),
            "--size_z", str(box_size["z"]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range),
            "--seed", str(seed),
            "--cpu", str(cpu),
            "--out", str(output_pdbqt)
        ]

        try:
            res = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout_seconds
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(f"Docking calculation timed out after {timeout_seconds} seconds.")

        if res.returncode != 0:
            err_msg = res.stderr.strip() or res.stdout.strip()
            raise RuntimeError(f"AutoDock Vina failed with code {res.returncode}: {err_msg}")

        # Parse modes table from stdout
        modes = self._parse_vina_output(res.stdout)
        if not modes:
            raise RuntimeError("AutoDock Vina completed but no binding modes were reported.")

        # Convert PDBQT poses to SDF for Mol* and ProLIF
        sdf_path_str = None
        if output_sdf:
            output_sdf.parent.mkdir(parents=True, exist_ok=True)
            obabel_cmd = self._get_obabel_cmd()
            conv_res = subprocess.run([
                obabel_cmd,
                "-ipdbqt", str(output_pdbqt),
                "-osdf", "-O", str(output_sdf)
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if conv_res.returncode == 0:
                sdf_path_str = str(output_sdf)

        return {
            "best_score": modes[0]["score"],
            "modes": modes,
            "pdbqt": str(output_pdbqt),
            "sdf": sdf_path_str,
            "settings": {
                "exhaustiveness": exhaustiveness,
                "seed": seed,
                "box_center": box_center,
                "box_size": box_size
            }
        }

    def _parse_vina_output(self, stdout: str) -> List[Dict[str, Any]]:
        modes = []
        table_pattern = re.compile(r"^\s*(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", re.MULTILINE)
        for match in table_pattern.finditer(stdout):
            mode_num = int(match.group(1))
            score = float(match.group(2))
            rmsd_lb = float(match.group(3))
            rmsd_ub = float(match.group(4))
            modes.append({
                "mode": mode_num,
                "score": score,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub
            })
        return modes

