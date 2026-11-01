"""Pose Consistency tab: reference choice, comparison, ranking, detail pane and export."""

import csv

import pytest
from PySide6.QtCore import Qt

Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend",
                           exc_type=ImportError)

from mds.chem.pose_consistency import CONSERVED
from mds.gui.controller import Controller
from mds.gui.panels.results_panel import ResultsPanel
from mds.jobs import JobRunner
from mds.project import Project
from conftest import wait_until


def _outcomes():
    """Two ligands: one docked on top of the reference, one displaced 9 A out of the site."""
    from rdkit.Chem import AllChem

    from mds.chem.workflow import LigandOutcome
    from mds.engines.base import Pose

    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1CCNC(=O)C"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0

    far = Chem.Mol(mol)
    conformer = far.GetConformer()
    for i in range(far.GetNumAtoms()):
        position = conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i, (position.x + 9.0, position.y, position.z))

    on_site = LigandOutcome("on site", "ok", poses=[Pose(1, -8.0, 0, 0, "")],
                            mols=[Chem.Mol(mol)], reference_rmsd=[None])
    off_site = LigandOutcome("off site", "ok", poses=[Pose(1, -7.0, 0, 0, "")],
                             mols=[far], reference_rmsd=[None])
    return [on_site, off_site], mol


@pytest.fixture
def panel(qapp, tmp_path, monkeypatch):
    controller = Controller(JobRunner())
    controller.set_project(Project.create(tmp_path / "project"))
    monkeypatch.setattr(controller, "contacts_for", lambda *_args: [], raising=False)
    monkeypatch.setattr(controller, "_require_chemistry", lambda: True)
    outcomes, mol = _outcomes()
    widget = ResultsPanel(controller)
    controller.outcomes = outcomes
    widget.show_outcomes(outcomes)
    widget._reference_mol = mol
    yield widget
    widget.close()


def _with_crystal_reference(panel):
    """Offer the embedded molecule as the co-crystal reference, as a prepared receptor would."""
    panel.c.reference = panel._reference_mol
    panel.c.reference_name = "LIG"
    panel.consistency_panel.refresh_references()
    return panel.consistency_panel


def _compare(tab, timeout=60):
    tab.compare_btn.click()
    assert wait_until(lambda: tab.table.rowCount() > 0, timeout), tab.status.text()
    return tab


def test_without_a_reference_the_panel_says_so_and_cannot_compare(panel):
    tab = panel.consistency_panel
    assert tab.reference.count() == 0
    assert not tab.compare_btn.isEnabled()
    assert "co-crystal ligand" in tab.status.text()


def test_crystal_ligand_is_offered_once_a_receptor_provides_one(panel):
    tab = _with_crystal_reference(panel)
    assert tab.reference.count() == 1
    assert tab.reference.currentData() == "crystal"
    assert "LIG" in tab.reference.currentText()
    assert tab.compare_btn.isEnabled()


def test_any_selected_pose_can_become_the_reference(panel):
    tab = panel.consistency_panel
    tab.set_pose(panel.c.outcomes[0], 0)
    assert tab.use_pose_btn.isEnabled()
    tab.use_pose_btn.click()
    assert tab.reference.currentData() == "pose"
    assert "on site #1" in tab.reference.currentText()
    assert tab.compare_btn.isEnabled()


def test_the_reference_pose_is_not_compared_against_itself(panel):
    tab = panel.consistency_panel
    tab.set_pose(panel.c.outcomes[0], 0)
    tab.use_pose_btn.click()
    _compare(tab)
    assert tab.table.rowCount() == 1
    assert tab.table.item(0, 0).text() == "off site #1"


def test_comparison_ranks_poses_by_consensus(panel):
    tab = _compare(_with_crystal_reference(panel))
    assert tab.table.rowCount() == 2
    assert tab.table.item(0, 0).text() == "on site #1"
    assert tab.table.item(0, 1).text() == CONSERVED
    assert float(tab.table.item(0, 2).text()) > float(tab.table.item(1, 2).text())
    assert tab.table.item(1, 0).text() == "off site #1"
    # Column 3 is the whole-molecule RMSD: zero on site, ~9 A for the displaced pose.
    assert float(tab.table.item(0, 3).text()) == pytest.approx(0.0, abs=0.01)
    assert float(tab.table.item(1, 3).text()) == pytest.approx(9.0, abs=0.1)
    assert tab.export_btn.isEnabled()


def test_the_interaction_axis_reports_why_it_is_unavailable(panel):
    """With no prepared receptor the fingerprint cannot run; it must say so, not score zero."""
    tab = _compare(_with_crystal_reference(panel))
    assert tab.table.item(0, 6).text() == "—"
    notes = tab.table.item(0, 0).toolTip()
    assert "Interaction axis unavailable" in notes
    # The remaining axes still rank the on-site pose as conserved.
    assert tab.table.item(0, 1).text() == CONSERVED


def test_selecting_a_row_shows_that_pose_and_its_detail(panel):
    tab = _compare(_with_crystal_reference(panel))
    shown = []
    panel.pose_selected.connect(lambda outcome, index: shown.append((outcome.name, index)))
    tab.table.selectRow(1)
    assert shown and shown[-1] == ("off site", 0)
    assert "off site #1" in tab.detail.toPlainText()
    assert "Notes" in tab.detail.toPlainText()


def test_new_results_clear_the_comparison(panel):
    tab = _compare(_with_crystal_reference(panel))
    outcomes, _mol = _outcomes()
    panel.show_outcomes(outcomes)
    assert tab.table.rowCount() == 0
    assert not tab.export_btn.isEnabled()
    assert tab.detail.toPlainText() == ""
    assert tab._pose_reference is None


def test_export_writes_every_axis_and_the_contact_lists(panel, monkeypatch, tmp_path):
    tab = _compare(_with_crystal_reference(panel))
    target = tmp_path / "consistency.csv"
    monkeypatch.setattr("mds.gui.panels.consistency_panel.QFileDialog.getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "CSV (*.csv)")))
    tab.export_btn.click()
    with target.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert {"Pose", "Verdict", "Consensus", "RMSD (A)", "MCS RMSD (A)", "Volume overlap (%)",
            "Conserved contacts", "Lost contacts", "New contacts", "Notes"} <= set(rows[0])
    assert {r["Pose"] for r in rows} == {"on site #1", "off site #1"}


def test_switching_axes_off_leaves_their_columns_blank(panel):
    tab = _with_crystal_reference(panel)
    tab.shape.setChecked(False)
    tab.features.setChecked(False)
    _compare(tab)
    row = next(r for r in range(tab.table.rowCount()) if tab.table.item(r, 0).text() == "on site #1")
    assert tab.table.item(row, 10).text() == "—"      # Volume (%)
    assert tab.table.item(row, 11).text() == "—"      # USR
    assert tab.table.item(row, 7).text() == "—"       # Features


def test_results_for_a_changed_pose_list_are_ignored(panel):
    tab = _with_crystal_reference(panel)
    tab._rows = [(panel.c.outcomes[0], 0)]        # pretend only one pose was requested
    panel.c.consistency_ready.emit("LIG", [("a", object()), ("b", object())])
    assert tab.table.rowCount() == 0


def test_a_failed_comparison_says_so_instead_of_staying_silent(panel):
    tab = _with_crystal_reference(panel)
    panel.c.consistency_busy.emit(True)
    panel.c.consistency_busy.emit(False)          # as a failed or cancelled job would
    assert "stopped before any result arrived" in tab.status.text()


def test_a_successful_comparison_never_shows_the_failure_message(panel):
    tab = _compare(_with_crystal_reference(panel))
    assert "stopped before any result arrived" not in tab.status.text()
    assert "compared with LIG" in tab.status.text()
    assert tab.compare_btn.isEnabled()
