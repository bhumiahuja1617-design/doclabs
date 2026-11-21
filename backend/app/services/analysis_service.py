"""Interaction summaries and labeled screening metrics."""
from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve


class AnalysisService:
    @staticmethod
    def enrichment_factor(labels: list[int], scores: list[float], fraction: float) -> float:
        if not 0 < fraction <= 1 or len(labels) != len(scores) or not labels:
            raise ValueError("Labels and scores must have equal nonzero length and fraction in (0, 1].")
        labels_arr = np.asarray(labels, dtype=int)
        if labels_arr.sum() == 0:
            raise ValueError("At least one active is required to calculate enrichment.")
        n = max(1, int(np.ceil(len(labels) * fraction)))
        top = np.argsort(np.asarray(scores))[:n]  # more negative Vina score ranks first
        return float((labels_arr[top].sum() / n) / (labels_arr.sum() / len(labels)))

    def screening_metrics(self, labels: list[int], scores: list[float]) -> dict[str, Any]:
        if len(labels) != len(scores) or len(set(labels)) != 2:
            raise ValueError("ROC-AUC requires equal-length scores with both active and decoy labels.")
        # Lower Vina energies indicate higher rank, so invert for ROC's positive direction.
        y = np.asarray(labels, dtype=int)
        positive_score = -np.asarray(scores, dtype=float)
        fpr, tpr, thresholds = roc_curve(y, positive_score)
        return {"n": len(y), "actives": int(y.sum()), "decoys": int(len(y) - y.sum()),
                "roc_auc": float(roc_auc_score(y, positive_score)), "ef1": self.enrichment_factor(list(y), list(scores), .01),
                "ef5": self.enrichment_factor(list(y), list(scores), .05), "ef10": self.enrichment_factor(list(y), list(scores), .10),
                "roc": [{"fpr": float(a), "tpr": float(b), "threshold": float(c)} for a, b, c in zip(fpr, tpr, thresholds)]}
