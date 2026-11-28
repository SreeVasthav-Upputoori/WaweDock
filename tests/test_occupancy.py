from types import SimpleNamespace

import pytest

from mds.chem.occupancy import (
    classify_residue,
    split_pdb_frames,
    write_occupancy_csv,
)


def test_occupancy_classification_uses_documented_thresholds():
    assert classify_residue(0.90, 1) == "anchor"
    assert classify_residue(0.60, 2) == "hotspot"
    assert classify_residue(0.40, 1) == "accessible"
    assert classify_residue(0.20, 1) == "transient"


def test_split_multi_model_pdb(tmp_path):
    pdb = tmp_path / "ensemble.pdb"
    pdb.write_text(
        "MODEL        1\nATOM      1  CA  ALA A   1       0.000   0.000   0.000\n"
        "ENDMDL\nMODEL        2\nATOM      1  CA  ALA A   1       1.000   0.000   0.000\n"
        "ENDMDL\n"
    )
    frames = split_pdb_frames(pdb)
    assert len(frames) == 2
    assert "MODEL" not in frames[0]
    assert "ATOM" in frames[1]


def test_split_directory_of_pdb_frames(tmp_path):
    (tmp_path / "frame-2.pdb").write_text("ATOM      1\n")
    (tmp_path / "frame-1.ent").write_text("ATOM      1\n")
    (tmp_path / "ignored.txt").write_text("ATOM      1\n")
    frames = split_pdb_frames(tmp_path)
    assert len(frames) == 2
    assert frames[0].strip() == "ATOM      1"


def test_analyze_occupancy_aggregates_frames(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from mds.chem import occupancy

    frame = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "HETATM    2  C1  LIG B   2       1.000   0.000   0.000  1.00 20.00           C\n"
    )
    pdb = tmp_path / "frames.pdb"
    pdb.write_text(
        "MODEL        1\n" + frame + "ENDMDL\n"
        "MODEL        2\n" + frame + "ENDMDL\n"
    )

    monkeypatch.setattr(
        "mds.chem.interactions.protein_to_rdkit",
        lambda path: object(),
    )
    monkeypatch.setattr(
        "mds.chem.interactions.compute_interactions",
        lambda protein, ligand, diagram=False: SimpleNamespace(
            residues=[
                SimpleNamespace(residue="ASP25:A", interactions=["H-bond (ligand donor)"]),
                SimpleNamespace(residue="TYR82:A", interactions=["Hydrophobic"]),
            ],
            n_interactions=2,
        ),
    )
    monkeypatch.setattr(Chem, "MolFromPDBFile", lambda *args, **kwargs: object())
    monkeypatch.setattr(
        "mds.chem.poses.assign_bond_orders_from_template",
        lambda template, raw: object(),
    )

    result = occupancy.analyze_occupancy(pdb, "LIG", "[C]")
    assert result.n_frames == 2
    assert len(result.residues) == 2
    assert result.residues[0].occupancy == 1.0
    assert result.residues[0].classification == "anchor"

    csv_path = write_occupancy_csv(result, tmp_path / "occupancy.csv")
    assert csv_path.exists()
    assert "Occupancy" in csv_path.read_text(encoding="utf-8")
