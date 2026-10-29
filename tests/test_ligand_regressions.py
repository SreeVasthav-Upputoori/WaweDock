"""Regression coverage for ligand chemistry and failure reporting."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem", exc_type=ImportError)
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

from mds.chem.ligand import (
    LigandEntry, _apply_simple_ionisation, _optimise,
    entry_from_smiles, prepare_ligand, read_ligands,
)
from mds.chem.report import PreparationError, Report


@pytest.mark.parametrize("smiles, expected", [
    ("NCCN", "[NH3+]CC[NH3+]"),
    ("[NH2]CC", "CC[NH3+]"),
    ("CC(=O)[OH]", "CC(=O)[O-]"),
    ("N=C(N)c1ccccc1", "NC(=[NH2+])c1ccccc1"),
    ("CC(=O)N", "CC(N)=O"),
    ("C[NH2+]C", "C[NH2+]C"),
])
def test_ionisation_hydrogens_all_sites_and_idempotence(smiles, expected):
    original = Chem.MolFromSmiles(smiles)
    before = Chem.MolToSmiles(original)
    result = _apply_simple_ionisation(original, Report())
    assert Chem.MolToSmiles(result) == expected
    assert Chem.MolToSmiles(_apply_simple_ionisation(result, Report())) == expected
    assert Chem.MolToSmiles(original) == before
    assert all(a.GetNumRadicalElectrons() == 0 for a in result.GetAtoms())


def test_csv_missing_fields_are_invalid_entries_not_batch_errors(tmp_path):
    path = tmp_path / "ligands.csv"
    path.write_text("id,smiles\nmissing\nblank,\ngood,CCO\n")
    entries = read_ligands(path)
    assert len(entries) == 3
    assert all(e.mol is None and e.error for e in entries[:2])
    assert entries[2].mol is not None


def test_empty_smiles_and_empty_molecule_rejected():
    assert entry_from_smiles("  ").mol is None
    with pytest.raises(PreparationError, match="empty ligand"):
        prepare_ligand(LigandEntry("empty", Chem.Mol()))


@pytest.mark.parametrize("options", [{"ionize": "typo"}, {"n_conformers": 0}, {"n_conformers": 1.5}])
def test_invalid_options_rejected(options):
    with pytest.raises(PreparationError) as error:
        prepare_ligand(entry_from_smiles("CCO"), **options)
    assert error.value.report.errors


def embedded():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert rdDistGeom.EmbedMolecule(mol, randomSeed=42) == 0
    return mol


def test_preserve_coordinates_even_with_default_optimisation():
    pytest.importorskip("meeko", exc_type=ImportError)
    mol = embedded()
    before = mol.GetConformer().GetPositions().copy()
    result = prepare_ligand(LigandEntry("ethanol", mol), keep_input_3d=True)
    np.testing.assert_array_equal(result.mol.GetConformer().GetPositions(), before)
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), before)
    assert "Optimisation skipped" in result.report.render()


def test_ionisation_report_describes_final_state():
    pytest.importorskip("meeko", exc_type=ImportError)
    result = prepare_ligand(entry_from_smiles("CCN"), ionize="simple")
    assert result.properties["FormalCharge"] == 1
    assert "input protonation state was kept" not in result.report.render()


def test_meeko_failure_becomes_reported_preparation_error(monkeypatch):
    meeko = pytest.importorskip("meeko", exc_type=ImportError)
    def fail(*args, **kwargs):
        raise RuntimeError("unsupported chemistry")
    monkeypatch.setattr(meeko.MoleculePreparation, "prepare", fail)
    with pytest.raises(PreparationError) as error:
        prepare_ligand(entry_from_smiles("CCO"))
    assert "unsupported chemistry" in error.value.report.render()


def test_uff_missing_parameters_rejected(monkeypatch):
    monkeypatch.setattr(rdForceFieldHelpers, "MMFFHasAllMoleculeParams", lambda m: False)
    monkeypatch.setattr(rdForceFieldHelpers, "UFFHasAllMoleculeParams", lambda m: False)
    with pytest.raises(PreparationError, match="missing force-field parameters"):
        _optimise(embedded(), Report())


def test_uff_nonconvergence_reported(monkeypatch):
    monkeypatch.setattr(rdForceFieldHelpers, "MMFFHasAllMoleculeParams", lambda m: False)
    monkeypatch.setattr(rdForceFieldHelpers, "UFFHasAllMoleculeParams", lambda m: True)
    monkeypatch.setattr(rdForceFieldHelpers, "UFFOptimizeMoleculeConfs", lambda *a, **kw: [(1, 10.0)])
    report = Report()
    _optimise(embedded(), report)
    assert any("did not converge" in w and "UFF" in w for w in report.warnings)
