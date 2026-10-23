"""Regression checks for interaction classification and receptor identity guards."""

import json
import sys
from types import SimpleNamespace

import pytest
Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend", exc_type=ImportError)

from mds.chem import interactions
from mds.chem.report import PreparationError


def _molecule(smiles, position):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (position + i * 0.1, 0, 0))
    mol.AddConformer(conf)
    return mol


def test_metal_coordination_matches_unfiltered_analysis(monkeypatch):
    import prolif

    protein = _molecule("N", 0)
    for atom in protein.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("LYS")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)
    ligand = _molecule("[Zn+2]", 2)
    filtered = interactions.compute_interactions(protein, ligand, diagram=False)
    assert "Metal coordination" in {
        label for residue in filtered.residues for label in residue.interactions
    }

    original = prolif.Fingerprint.run_from_iterable

    def unfiltered(self, *args, **kwargs):
        return original(self, *args, **kwargs, residues="all")

    monkeypatch.setattr(prolif.Fingerprint, "run_from_iterable", unfiltered)
    assert interactions.compute_interactions(protein, ligand, diagram=False) == filtered


def test_receptor_preserves_insertion_codes_with_unique_internal_ids(tmp_path, monkeypatch):
    converted = SimpleNamespace(write=lambda _: "mock mol2")
    pybel = SimpleNamespace(readfile=lambda *_: iter([converted]))
    monkeypatch.setitem(sys.modules, "openbabel", SimpleNamespace(pybel=pybel))
    mol = _molecule("[C]", 0)
    monkeypatch.setattr(Chem, "MolFromMol2Block", lambda *_, **__: mol)
    pdb = tmp_path / "insertion.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A  25A      0.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    )
    protein = interactions.protein_to_rdkit(pdb)
    info = protein.GetAtomWithIdx(0).GetPDBResidueInfo()
    assert info.GetResidueNumber() == 26
    assert info.GetInsertionCode() == "A"
    assert json.loads(protein.GetProp("_mds_residue_labels")) == {"ALA26.A": "ALA25A.A"}


def test_inserted_residue_labels_in_table_and_diagram():
    protein = _molecule("N", 0)
    for atom in protein.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("LYS")
        info.SetResidueNumber(26)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)
    protein.SetProp("_mds_residue_labels", json.dumps({"LYS26.A": "LYS25A.A"}))
    result = interactions.compute_interactions(protein, _molecule("[Zn+2]", 2))
    assert result.residues[0].residue == "LYS25A:A"
    assert result.network_html and '"LYS25A.A"' in result.network_html
    assert '"LYS26.A"' not in result.network_html


def test_receptor_rejects_empty_structure(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "openbabel", SimpleNamespace(pybel=SimpleNamespace()))
    pdb = tmp_path / "empty.pdb"
    pdb.write_text("END\n")
    with pytest.raises(PreparationError, match="No protein atoms"):
        interactions.protein_to_rdkit(pdb)


def test_receptor_rejects_atom_position_mismatch(tmp_path, monkeypatch):
    converted = SimpleNamespace(write=lambda _: "mock mol2")
    pybel = SimpleNamespace(readfile=lambda *_: iter([converted]))
    monkeypatch.setitem(sys.modules, "openbabel", SimpleNamespace(pybel=pybel))
    mol = _molecule("[C]", 5)
    monkeypatch.setattr(Chem, "MolFromMol2Block", lambda *_, **__: mol)
    pdb = tmp_path / "reordered.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A  25       0.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    )
    with pytest.raises(PreparationError, match="Atom order mismatch"):
        interactions.protein_to_rdkit(pdb)
