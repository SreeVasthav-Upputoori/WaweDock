"""Pharmacophore tab: deriving from the selected pose, scoring poses, viewer display."""

import json

import pytest
from PySide6.QtCore import Qt

Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend",
                           exc_type=ImportError)

from mds.chem.interactions import InteractionContact, InteractionResult
from mds.chem.pharmacophore import derive
from mds.gui.controller import Controller
from mds.gui.panels.results_panel import ResultsPanel
from mds.jobs import JobRunner
from mds.project import Project
from conftest import wait_until


def _phenol_outcomes():
    """Two ligands whose poses are a phenol and the same phenol shifted out of the site."""
    from rdkit.Chem import AllChem

    from mds.chem.workflow import LigandOutcome
    from mds.engines.base import Pose

    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0

    far = Chem.Mol(mol)
    conformer = far.GetConformer()
    for i in range(far.GetNumAtoms()):
        position = conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i, (position.x + 8.0, position.y, position.z))

    on_site = LigandOutcome("on site", "ok", poses=[Pose(1, -8.0, 0, 0, "")],
                            mols=[Chem.Mol(mol)], reference_rmsd=[None])
    off_site = LigandOutcome("off site", "ok", poses=[Pose(1, -7.0, 0, 0, "")],
                             mols=[far], reference_rmsd=[None])
    return [on_site, off_site], mol


def _interaction_result(mol):
    """What ProLIF would report for the phenol: a donor, an acceptor, a stacked ring."""
    oxygen = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O")
    ring = tuple(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    contacts = [
        InteractionContact("ASP25:A", "HBDonor", "H-bond (ligand donor)", (oxygen,)),
        InteractionContact("THR26:A", "HBAcceptor", "H-bond (ligand acceptor)", (oxygen,)),
        InteractionContact("PHE53:A", "PiStacking", "Pi-stacking", ring),
        InteractionContact("ILE50:A", "VdWContact", "Van der Waals contact", (oxygen,)),
    ]
    return InteractionResult(residues=[], network_html=None, n_interactions=len(contacts),
                             contacts=contacts)


@pytest.fixture
def panel(qapp, tmp_path, monkeypatch):
    controller = Controller(JobRunner())
    controller.set_project(Project.create(tmp_path / "project"))
    monkeypatch.setattr(controller, "contacts_for", lambda *_args: [])
    monkeypatch.setattr(controller, "_require_chemistry", lambda: True)
    outcomes, mol = _phenol_outcomes()
    widget = ResultsPanel(controller)
    controller.outcomes = outcomes
    widget.show_outcomes(outcomes)
    widget._mol = mol
    yield widget
    widget.close()


def _derive_through_panel(panel, monkeypatch):
    """Press Derive with the ProLIF step stubbed, and hand back the interaction result."""
    tab = panel.pharmacophore_panel
    requested = []
    monkeypatch.setattr(panel.c, "analyze_interactions",
                        lambda outcome, index: requested.append((outcome, index)))
    tab.derive_btn.click()
    assert requested, "Derive should ask the controller for the ProLIF classification"
    outcome, index = requested[0]
    panel.c.interactions_ready.emit(outcome, index, _interaction_result(panel._mol))
    return tab


def test_derive_builds_features_from_the_selected_pose(panel, monkeypatch):
    tab = _derive_through_panel(panel, monkeypatch)
    assert tab.features.rowCount() == 3            # the van der Waals contact is not a feature
    labels = {tab.features.item(r, 0).text() for r in range(3)}
    assert labels == {"H-bond donor", "H-bond acceptor", "Aromatic ring"}
    assert "on site pose #1" in tab.status.text()
    assert "-8.00 kcal/mol, predicted" in tab._pharmacophore.source
    for row in range(3):
        assert float(tab.features.item(row, 2).text()) > 0
    assert tab.score_btn.isEnabled() and tab.export_btn.isEnabled()
    assert tab.show_in_viewer.isEnabled()


def test_interaction_results_for_other_poses_are_ignored(panel, monkeypatch):
    """A ProLIF result the user requested from the Interactions tab must not become a pharmacophore."""
    tab = panel.pharmacophore_panel
    other = panel.c.outcomes[1]
    panel.c.interactions_ready.emit(other, 0, _interaction_result(panel._mol))
    assert tab._pharmacophore is None and tab.features.rowCount() == 0
    assert not tab.score_btn.isEnabled()


def test_scoring_ranks_poses_by_feature_recovery(panel, monkeypatch):
    tab = _derive_through_panel(panel, monkeypatch)
    tab.score_btn.click()
    assert wait_until(lambda: tab.matches.rowCount() == 2, 20)

    # Sorted by recovery: the pose the pharmacophore came from first, the displaced one last.
    assert tab.matches.item(0, 0).text() == "on site"
    assert tab.matches.item(0, 2).text() == "100%"
    assert tab.matches.item(0, 3).text() == "3 / 3"
    assert tab.matches.item(0, 4).text() == "—"
    assert tab.matches.item(1, 0).text() == "off site"
    assert tab.matches.item(1, 2).text() == "0%"
    assert set(tab.matches.item(1, 4).text().split(", ")) == {
        "H-bond donor", "H-bond acceptor", "Aromatic ring"}
    assert "not a predicted affinity" in tab.status.text()


def test_selecting_a_scored_row_shows_that_pose(panel, monkeypatch):
    tab = _derive_through_panel(panel, monkeypatch)
    tab.score_btn.click()
    assert wait_until(lambda: tab.matches.rowCount() == 2, 20)
    shown = []
    panel.pose_selected.connect(lambda outcome, index: shown.append((outcome.name, index)))
    tab.matches.selectRow(1)          # the displaced pose
    assert shown and shown[-1] == ("off site", 0)


def test_viewer_display_is_opt_in_and_reversible(panel, monkeypatch):
    tab = _derive_through_panel(panel, monkeypatch)
    events = []
    tab.display_requested.connect(lambda pharmacophore, show: events.append(
        (None if pharmacophore is None else len(pharmacophore.features), show)))
    tab.show_in_viewer.setChecked(True)
    assert events[-1] == (3, True)
    tab.show_in_viewer.setChecked(False)
    assert events[-1] == (None, False)


def test_new_results_clear_the_pharmacophore(panel, monkeypatch):
    tab = _derive_through_panel(panel, monkeypatch)
    cleared = []
    tab.display_requested.connect(lambda pharmacophore, show: cleared.append((pharmacophore, show)))
    outcomes, _mol = _phenol_outcomes()
    panel.show_outcomes(outcomes)
    assert tab._pharmacophore is None
    assert tab.features.rowCount() == 0 and tab.matches.rowCount() == 0
    assert not tab.score_btn.isEnabled() and not tab.export_btn.isEnabled()
    assert cleared[-1] == (None, False)   # the viewer is told to clear its spheres


def test_export_writes_features_and_the_caveat(panel, monkeypatch, tmp_path):
    tab = _derive_through_panel(panel, monkeypatch)
    target = tmp_path / "pharm.json"
    monkeypatch.setattr("mds.gui.panels.pharmacophore_panel.QFileDialog.getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "JSON (*.json)")))
    tab.export_btn.click()
    data = json.loads(target.read_text(encoding="utf-8"))
    assert len(data["features"]) == 3
    assert "not a validated activity model" in data["caveat"]
    assert {f["kind"] for f in data["features"]} == {"HBD", "HBA", "Aromatic"}
    assert all({"x", "y", "z", "tolerance", "residues"} <= set(f) for f in data["features"])


def test_viewer_payload_carries_position_and_tolerance(qapp, monkeypatch):
    """The viewer call sends plain feature dicts, so it needs no chemistry types in JS."""
    from mds.render.mol3d import Mol3DViewer

    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    pharmacophore = derive(_interaction_result(mol).contacts, mol)

    calls = []
    monkeypatch.setattr(Mol3DViewer, "_call", lambda self, fn, *args: calls.append((fn, args)))
    viewer = Mol3DViewer.__new__(Mol3DViewer)
    viewer.show_pharmacophore(pharmacophore.features)
    function, (payload, labelled) = calls[-1]
    assert function == "showPharmacophore" and labelled is True
    assert len(payload) == 3 and all(isinstance(entry, dict) for entry in payload)
    assert {"kind", "label", "tolerance", "x", "y", "z"} == set(payload[0])
    viewer.clear_pharmacophore()
    assert calls[-1] == ("clearPharmacophore", ())


def test_pharmacophore_renders_in_the_native_viewer(qapp, caplog):
    """Live WebGL check: a JS typo in showPharmacophore would otherwise fail silently."""
    from rdkit.Chem import AllChem

    from mds.render.mol3d import Mol3DViewer

    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    pharmacophore = derive(_interaction_result(mol).contacts, mol)

    viewer = Mol3DViewer()
    viewer.resize(500, 400)
    viewer.show()
    caplog.set_level("ERROR", logger="mds.render")

    def js(expression):
        result = []
        viewer.page().runJavaScript(expression, result.append)
        assert wait_until(lambda: bool(result), 10)
        return result[0]

    try:
        assert wait_until(lambda: viewer._ready, 20)
        assert js("typeof showPharmacophore") == "function"
        assert js("pharmShapes.length") == 0

        viewer.show_pharmacophore(pharmacophore.features)
        # one sphere + one label per feature
        assert wait_until(lambda: js("pharmShapes.length") == 6, 10)

        viewer.show_pharmacophore(pharmacophore.features, labelled=False)
        assert wait_until(lambda: js("pharmShapes.length") == 3, 10)   # redraw replaces, not appends

        viewer.clear_pharmacophore()
        assert wait_until(lambda: js("pharmShapes.length") == 0, 10)
        script_errors = [r.getMessage() for r in caplog.records
                         if "Viewer script error" in r.getMessage()]
        assert not script_errors, f"viewer reported JS errors: {script_errors}"
    finally:
        viewer.close()
