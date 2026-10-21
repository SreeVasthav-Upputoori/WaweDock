"""End-to-end test of the GUI workflow through the controller (real preparation + real Vina)."""

from __future__ import annotations

from pathlib import Path

import pytest
from conftest import wait_until

from mds.chem.receptor import ReceptorOptions
from mds.chem.workflow import DockingParameters
from mds.engines import VinaEngine
from mds.gui.main_window import MainWindow
from mds.project import Project

DATA = Path(__file__).parent / "data"
pytestmark = pytest.mark.skipif(not VinaEngine().available(), reason="vina.exe not bundled")


@pytest.fixture
def window(qapp, settings, tmp_path):
    w = MainWindow(settings, use_viewer=False)
    w.controller.set_project(Project.create(tmp_path / "proj"))
    yield w
    w.job_runner.wait_for_done(10000)
    w.close()


def test_full_docking_workflow_and_project_reload(window, tmp_path):
    c = window.controller
    errors = []
    c.error.connect(lambda t, m: errors.append((t, m)))

    c.load_structure(DATA / "3PTB.pdb")
    assert c.source_pdb is not None and (c.project.root / "originals" / "3PTB.pdb").exists()
    assert "3PTB" in window.project_explorer.object_names("Proteins")

    c.prepare_receptor(ReceptorOptions())
    assert wait_until(lambda: c.receptor is not None, 60), errors
    assert any(h.name == "BEN" for h in c.hetero)

    assert c.add_smiles("N=C(N)c1ccccc1", "benzamidine")
    c.prepare_ligands(ionize="simple", n_conformers=1, seed=1)
    assert wait_until(lambda: len(c.prepared) == 1, 60), errors

    ben = next(h for h in c.hetero if h.name == "BEN")
    c.box_from_hetero(ben, 8.0)
    assert c.validate_box().errors == []

    c.run_docking(DockingParameters(exhaustiveness=8, num_modes=5, cpu=2, seed=3))
    assert wait_until(lambda: bool(c.outcomes), 120), errors
    outcome = c.outcomes[0]
    assert outcome.status == "ok" and outcome.poses[0].score < -4
    assert window.results_panel.table.rowCount() == len(outcome.poses)
    assert not errors

    # contacts are distance-based candidates with a caveat, never "interactions"
    contacts = c.contacts_for(outcome, 0)
    assert contacts and any("candidate" in x for cc in contacts for x in cc.candidates)

    # ProLIF interaction fingerprint + 2D diagram, run through the controller as a background job
    interaction_results = []
    c.interactions_ready.connect(lambda o, i, r: interaction_results.append(r))
    c.analyze_interactions(outcome, 0)
    assert wait_until(lambda: bool(interaction_results), 60), errors
    result = interaction_results[0]
    assert result.n_interactions > 0 and result.network_html and "unpkg.com" not in result.network_html
    assert not errors

    # export formats
    c.export_csv(tmp_path / "r.csv")
    c.export_sdf(tmp_path / "r.sdf")
    c.export_excel(tmp_path / "r.xlsx")
    assert all((tmp_path / n).stat().st_size > 0 for n in ("r.csv", "r.sdf", "r.xlsx"))

    # save, then reopen the project in a fresh controller state
    c.save_project()
    root = c.project.root
    reopened = Project.open(root)
    c.reset_state()
    c.set_project(reopened)
    c.load_project_state()
    assert c.receptor is not None and c.box is not None
    assert [e.name for e in c.entries] == ["benzamidine"]
    assert c.outcomes and c.outcomes[0].poses[0].score == outcome.poses[0].score


def test_docking_refuses_to_start_without_prerequisites(window):
    c = window.controller
    errors = []
    c.error.connect(lambda t, m: errors.append((t, m)))
    c.run_docking(DockingParameters())
    assert errors and "Prepare the receptor" in errors[0][1] and "Define a docking box" in errors[0][1]


def test_similarity_via_controller_and_dialog(qapp, settings, tmp_path):
    """Similarity/clustering needs no receptor or Vina — only prepared ligands."""
    w = MainWindow(settings, use_viewer=False)
    w.controller.set_project(Project.create(tmp_path / "simproj"))
    c = w.controller
    errors = []
    c.error.connect(lambda t, m: errors.append((t, m)))
    try:
        assert c.add_smiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
        assert c.add_smiles("Oc1ccccc1C(=O)O", "salicylic_acid")
        assert c.add_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "caffeine")
        c.prepare_ligands(ionize="none", n_conformers=1, seed=1)
        assert wait_until(lambda: len(c.prepared) == 3, 60), errors
        assert not errors

        results = []
        c.similarity_ready.connect(results.append)
        c.compute_similarity(cutoff=0.6)
        assert wait_until(lambda: bool(results), 30), errors
        result = results[0]
        by_name = dict(zip(result.names, result.clusters))
        assert by_name["aspirin"] == by_name["salicylic_acid"] != by_name["caffeine"]

        # dialog opens and reflects a just-computed result without erroring
        w._show_similarity()
        assert w._similarity_dialog is not None and w._similarity_dialog.isVisible()
        assert not errors
    finally:
        w.job_runner.wait_for_done(10000)
        w.close()
