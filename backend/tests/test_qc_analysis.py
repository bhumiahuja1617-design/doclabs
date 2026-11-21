from pathlib import Path

import pytest

from backend.app.services.analysis_service import AnalysisService
from backend.app.services.qc_service import QCService
from backend.app.services.reliability_service import ReliabilityService


def test_screening_metrics_match_expected_ordering():
    metrics = AnalysisService().screening_metrics(
        [1, 1, 0, 0], [-9.0, -8.0, -6.0, -5.0]
    )
    assert metrics["roc_auc"] == 1.0
    assert metrics["ef5"] == 2.0


def test_ef_is_invalid_without_actives():
    with pytest.raises(ValueError, match="At least one active"):
        AnalysisService.enrichment_factor([0, 0], [-4.0, -5.0], .05)


def test_qc_flags_invalid_ligand_and_out_of_range_box(tmp_path: Path):
    protein = tmp_path / "protein.pdb"
    protein.write_text(
        "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C\nEND\n"
    )
    report = QCService().run(
        protein,
        [{"name": "bad", "smiles": "not-a-smiles"}],
        {"center": {"x": 10, "y": 10, "z": 10}, "size": {"x": 31, "y": 20, "z": 20}},
    )
    by_id = {item["id"]: item for item in report["checks"]}
    assert by_id["ligand_inputs"]["severity"] == "fail"
    assert by_id["box_bounds"]["severity"] == "fail"
    assert report["requires_acknowledgement"] is True


def test_rank_stability_handles_identical_and_shuffled_rankings():
    stable = ReliabilityService.rank_stability({"a": [-9, -8], "b": [-7, -6], "c": [-5, -4]})
    assert stable["mean_spearman"] == pytest.approx(1.0)
    assert stable["rank_ranges"]["a"] == [1.0, 1.0]
    changed = ReliabilityService.rank_stability({"a": [-9, -5], "b": [-7, -7], "c": [-5, -9]})
    assert changed["mean_spearman"] == pytest.approx(-1.0)
