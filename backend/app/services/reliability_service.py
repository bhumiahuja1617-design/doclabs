"""Multi-seed Vina repeatability analysis."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import shutil

import numpy as np
from scipy.stats import spearmanr
from scipy.stats import rankdata

from .docking_service import DockingService
from .validation_service import ValidationService


class ReliabilityService:
    def __init__(self, docking: DockingService | None = None):
        self.docking = docking or DockingService()

    @staticmethod
    def rank_stability(seed_scores: dict[str, list[float]]) -> dict[str, Any]:
        """Calculate each ligand's rank range and mean pairwise seed-rank Spearman rho."""
        if not seed_scores:
            return {"rank_ranges": {}, "mean_spearman": None, "seed_count": 0}
        names = list(seed_scores)
        lengths = {len(scores) for scores in seed_scores.values()}
        if len(lengths) != 1 or next(iter(lengths)) < 2:
            raise ValueError("Each ligand must have the same scores for at least two seeds.")
        n_seeds = next(iter(lengths))
        ranks = []
        for seed_idx in range(n_seeds):
            ordered = np.asarray([seed_scores[name][seed_idx] for name in names], dtype=float)
            ranks.append(rankdata(ordered, method="average"))
        rank_ranges = {name: [float(min(row[i] for row in ranks)), float(max(row[i] for row in ranks))] for i, name in enumerate(names)}
        correlations = [float(spearmanr(ranks[i], ranks[j]).statistic) for i in range(n_seeds) for j in range(i + 1, n_seeds)]
        return {"rank_ranges": rank_ranges, "mean_spearman": float(np.mean(correlations)) if correlations else 1.0, "seed_count": n_seeds}

    @staticmethod
    def summarize_seed_runs(scores: list[float], poses: list[Path], validation: ValidationService | None = None) -> dict[str, Any]:
        if len(scores) < 2 or len(scores) != len(poses):
            raise ValueError("Provide matching scores and at least two pose files.")
        validation = validation or ValidationService()
        clusters: list[list[int]] = []
        for idx, pose in enumerate(poses):
            placed = False
            for cluster in clusters:
                ref = poses[cluster[0]]
                rmsd = validation.compute_symmetric_rmsd(ref, pose) if pose.suffix.lower() == ".sdf" else validation.compute_rmsd(ref, pose)
                if rmsd["rmsd_angstroms"] <= 2.0:
                    cluster.append(idx)
                    placed = True
                    break
            if not placed:
                clusters.append([idx])
        consistency = max(map(len, clusters)) / len(poses)
        sd = float(np.std(scores, ddof=1))
        label = "high" if consistency >= .8 and sd <= .3 else ("moderate" if consistency >= .5 and sd <= .7 else "low")
        return {"n_seeds": len(scores), "scores": scores, "score_mean": float(np.mean(scores)), "score_sd": sd,
                "pose_consistency": consistency, "largest_cluster": max(map(len, clusters)), "rank_range": None,
                "label": label, "explanation": f"Poses agreed in {max(map(len, clusters))} of {len(poses)} runs; score SD was {sd:.2f} kcal/mol. This measures computational repeatability, not biology."}

    def run_ligand(self, receptor: Path, ligand: Path, center: dict[str, float], size: dict[str, float], out_dir: Path,
                   seeds: list[int], exhaustiveness: int = 8, reference_scores: list[float] | None = None,
                   initial_score: float | None = None, initial_pose_sdf: Path | None = None) -> dict[str, Any]:
        if not 2 <= len(seeds) <= 10:
            raise ValueError("Reliability analysis requires 2 to 10 explicit seeds.")
        out_dir.mkdir(parents=True, exist_ok=True)
        scores, poses = [], []
        for seed in seeds:
            pose = out_dir / f"seed_{seed}.pdbqt"
            sdf = out_dir / f"seed_{seed}.sdf"
            if seed == seeds[0] and initial_score is not None and initial_pose_sdf is not None and initial_pose_sdf.exists():
                shutil.copy2(initial_pose_sdf, sdf)
                scores.append(float(initial_score))
                poses.append(sdf)
                continue
            result = self.docking.dock(receptor, ligand, center, size, pose, sdf, exhaustiveness=exhaustiveness, seed=seed)
            scores.append(float(result["best_score"]))
            poses.append(sdf if result.get("sdf") else pose)
        summary = self.summarize_seed_runs(scores, poses)
        summary["seeds"] = seeds
        if reference_scores is not None and len(reference_scores) == len(scores):
            summary["mean_spearman"] = float(spearmanr(reference_scores, scores).statistic)
        return summary
