"""
backend/tests/test_science_core.py
Unit tests for DockLab science core services: box geometry, SMILES parsing,
drug-likeness descriptors, and RMSD calculations.
"""

import pytest
import numpy as np
from pathlib import Path

from backend.app.services.box_service import BoxService
from backend.app.services.ligand_service import LigandService
from backend.app.services.validation_service import ValidationService

def test_box_calculation_math(tmp_path):
    # Create synthetic PDB coordinates
    pdb_path = tmp_path / "test_ligand.pdb"
    with open(pdb_path, "w") as f:
        f.write("HETATM    1  C1  LIG A   1      10.000  20.000  30.000  1.00  0.00           C\n")
        f.write("HETATM    2  C2  LIG A   1      20.000  25.000  35.000  1.00  0.00           C\n")

    box_svc = BoxService()
    box = box_svc.box_from_crystal_ligand(pdb_path, padding=10.0)

    # Centroid: (15.0, 22.5, 32.5)
    assert box["center"]["x"] == 15.0
    assert box["center"]["y"] == 22.5
    assert box["center"]["z"] == 32.5

    # Extent: (10.0, 5.0, 5.0) -> +10 padding = (20.0, 15.0, 15.0) -> clamped to min 15.0
    assert box["size"]["x"] == 20.0
    assert box["size"]["y"] == 15.0
    assert box["size"]["z"] == 15.0

def test_manual_box_validation():
    box_svc = BoxService()
    # Valid box
    res = box_svc.validate_manual_box(
        center={"x": 10.0, "y": 10.0, "z": 10.0},
        size={"x": 20.0, "y": 20.0, "z": 20.0}
    )
    assert res["size"]["x"] == 20.0

    # Overly large box triggers exception (> 30 A)
    with pytest.raises(ValueError, match="exceeds maximum permitted limit"):
        box_svc.validate_manual_box(
            center={"x": 0.0, "y": 0.0, "z": 0.0},
            size={"x": 35.0, "y": 20.0, "z": 20.0}
        )

def test_smiles_parsing_and_drug_likeness():
    lig_svc = LigandService()
    # Aspirin: CC(=O)Oc1ccccc1C(=O)O
    mol = lig_svc.parse_and_standardize_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert mol is not None

    props = lig_svc.compute_drug_likeness(mol)
    assert props["mw"] == pytest.approx(180.16, abs=0.5)
    assert props["lipinski_violations"] == 0
    assert props["lipinski_pass"] is True
    assert props["heavy_atom_count"] == 13

def test_invalid_smiles_rejection():
    lig_svc = LigandService()
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        lig_svc.parse_and_standardize_smiles("INVALID_CHEMISTRY_XYZ123")

def test_rmsd_direct_calculation(tmp_path):
    # Two identical sets of coordinates should produce 0.0 RMSD
    pdbqt_1 = tmp_path / "pose1.pdbqt"
    pdbqt_2 = tmp_path / "pose2.pdbqt"

    content = (
        "ATOM      1  C   UNL     1      10.000  10.000  10.000  1.00  0.00    +0.000 C \n"
        "ATOM      2  N   UNL     1      11.000  11.000  11.000  1.00  0.00    +0.000 N \n"
    )
    with open(pdbqt_1, "w") as f: f.write(content)
    with open(pdbqt_2, "w") as f: f.write(content)

    val_svc = ValidationService()
    res = val_svc.compute_rmsd(pdbqt_1, pdbqt_2)
    assert res["rmsd_angstroms"] == 0.0
    assert res["status"] == "PASS"

