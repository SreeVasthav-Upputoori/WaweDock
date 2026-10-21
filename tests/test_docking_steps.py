"""Docking steps preserve inputs and report controller state without running chemistry."""
from pathlib import Path
from types import SimpleNamespace

import pytest

from mds.chem.box import DockingBox
from mds.chem.report import Report
from mds.gui.controller import Controller
from mds.gui.panels.docking_panel import DockingPanel
from mds.jobs import JobRunner
from mds.project import Project


@pytest.fixture
def panel(qapp, tmp_path, monkeypatch):
    controller = Controller(JobRunner())
    controller.set_project(Project.create(tmp_path / 'project'))
    monkeypatch.setattr(controller, 'validate_box', lambda: Report())
    widget = DockingPanel(controller)
    yield widget
    widget.close()


def test_only_selected_step_expands_and_inputs_survive_collapse(panel):
    assert panel.current_step == 0
    assert [not page.isHidden() for page in panel.step_pages] == [True, False, False, False, False]
    panel.ph.setValue(6.8)
    panel.exh.setValue(32)
    panel.step_buttons[3].click()
    assert panel.current_step == 3
    assert sum(not page.isHidden() for page in panel.step_pages) == 1
    panel.step_buttons[3].click()
    assert panel.current_step == -1
    assert all(page.isHidden() for page in panel.step_pages)
    panel.show_step(0)
    assert panel.ph.value() == 6.8
    assert panel.exh.value() == 32


def test_every_step_fits_a_narrow_sidebar(panel, qapp):
    panel.resize(360, 800)
    panel.show()
    for index in range(5):
        panel.show_step(index)
        qapp.processEvents()
        assert panel.horizontalScrollBar().maximum() == 0
        assert panel.step_buttons[index].height() >= panel.fontMetrics().lineSpacing() * 3


def test_completion_and_partial_preparation_follow_project_state(panel):
    c = panel.c
    assert panel.step_status == ['Load a structure', 'Add ligands', 'Define a box', 'Not run', 'No poses yet']
    c.source_pdb = Path('protein.pdb')
    c.structure_loaded.emit('protein', str(c.source_pdb))
    assert panel.step_status[0] == 'Not prepared'
    c.receptor = SimpleNamespace(n_atoms=20, report=Report())
    c.receptor_prepared.emit(c.receptor)
    assert '\u2713' in panel.step_labels[0].text()
    c.entries = [SimpleNamespace(name='a'), SimpleNamespace(name='b')]
    c.prepared = [SimpleNamespace(name='a')]
    c.ligands_changed.emit()
    assert panel.step_status[1] == '1/2 prepared'
    assert '\u2713' not in panel.step_labels[1].text()
    c.prepared.append(SimpleNamespace(name='b'))
    c.ligands_prepared.emit(SimpleNamespace())
    assert '\u2713' in panel.step_labels[1].text()
    c.set_box(DockingBox((0, 0, 0), (20, 20, 20)))
    assert panel.step_status[2] == 'Ready'
    assert panel.run_readiness.text().startswith('Ready to dock')
    # A structural edit invalidates receptor and box completion.
    c.receptor = None
    c.set_box(None)
    c.structure_edited.emit('protein', 'edited.pdb')
    assert panel.step_status[0] == 'Not prepared'
    assert panel.step_status[2] == 'Define a box'


def test_invalid_box_never_shows_complete(panel, monkeypatch):
    report = Report()
    report.error('Box does not contain the binding site')
    monkeypatch.setattr(panel.c, 'validate_box', lambda: report)
    panel.c.set_box(DockingBox((0, 0, 0), (20, 20, 20)))
    assert panel.step_status[2] == 'Needs attention'
    assert '\u2713' not in panel.step_labels[2].text()


def test_running_cancel_failure_and_results_review(panel):
    c = panel.c
    reviews = []
    panel.results_requested.connect(lambda: reviews.append(True))
    assert not panel.review_btn.isEnabled()
    c.docking_started.emit()
    assert panel.step_status[3] == 'Running'
    assert not panel.run_btn.isEnabled() and panel.cancel_btn.isEnabled()
    c.runner.job_cancelled.emit(1)
    assert panel.step_status[3] == 'Cancelled'
    assert panel.run_btn.isEnabled() and not panel.cancel_btn.isEnabled()
    c.docking_started.emit()
    c.runner.job_failed.emit(2, 'test failure')
    assert panel.step_status[3] == 'Failed'
    c.outcomes = [SimpleNamespace(status='ok', poses=[object(), object()])]
    c.docking_finished.emit(c.outcomes)
    assert panel.step_status[3:] == ['Completed', '2 poses available']
    panel.show_step(4)
    panel.review_btn.click()
    assert reviews == [True]
    # Partial failures remain reviewable but never get a completion checkmark.
    c.outcomes.append(SimpleNamespace(status='failed', poses=[]))
    c.docking_finished.emit(c.outcomes)
    assert panel.step_status[3] == 'Needs attention'
    assert '\u2713' not in panel.step_labels[3].text()
    assert panel.review_btn.isEnabled()
    c.new_project(c.project.root.parent / 'new-project')
    assert panel.step_status == ['Load a structure', 'Add ligands', 'Define a box', 'Not run', 'No poses yet']
    assert not panel.review_btn.isEnabled()
    assert panel.lig_list.count() == 0


def test_existing_actions_forward_all_parameters(panel, monkeypatch):
    c = panel.c
    calls = []
    monkeypatch.setattr(c, 'prepare_receptor', lambda options: calls.append(options))
    monkeypatch.setattr(c, 'prepare_ligands', lambda **options: calls.append(options))
    monkeypatch.setattr(c, 'run_docking', lambda options: calls.append(options))
    panel.ph.setValue(6.5)
    panel.chains.setText('A, B')
    panel.keep_metals.setChecked(False)
    panel._prepare_receptor()
    assert panel.current_step == 0
    assert calls[-1].ph == 6.5 and calls[-1].chains == ['A', 'B']
    assert not calls[-1].keep_metals
    panel.ionize.setCurrentIndex(1)
    panel.confs.setValue(3)
    panel.lig_seed.setValue(100)
    panel._prepare_ligands()
    assert panel.current_step == 1
    assert calls[-1] == dict(ionize='simple', n_conformers=3, seed=100)
    panel.exh.setValue(32)
    panel.modes.setValue(5)
    panel.seed.setValue(123)
    panel.scoring.setCurrentIndex(1)
    panel._run()
    assert panel.current_step == 3
    assert (calls[-1].exhaustiveness, calls[-1].num_modes, calls[-1].seed, calls[-1].scoring) == (32, 5, 123, 'vinardo')
