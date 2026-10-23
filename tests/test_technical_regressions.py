"""Technical regression checks for job lifecycle and complete result exports."""
import csv
from types import SimpleNamespace

import openpyxl
import pytest
Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend", exc_type=ImportError)

from mds.chem.box import DockingBox
from mds.chem.report import Report
from mds.chem.workflow import DockingParameters, _write_ligand_outputs, results_rows
from mds.gui.controller import Controller
from mds.jobs import JobRunner
from mds.project import Project
from test_results_workspace import sample_outcomes


@pytest.fixture
def controller(qapp, tmp_path):
    c = Controller(JobRunner())
    c.set_project(Project.create(tmp_path / 'project'))
    return c


def test_duplicate_docking_is_refused_and_started_observers_can_cancel(controller, monkeypatch):
    c = controller
    c.receptor = SimpleNamespace(pdbqt=c.project.root / 'receptor.pdbqt')
    c.prepared = [SimpleNamespace(name='ligand')]
    c.box = DockingBox((0, 0, 0), (20, 20, 20))
    monkeypatch.setattr(c, 'validate_box', lambda: Report())
    submitted, cancelled, started = [], [], []
    def submit(job, callback):
        job_id = 100 + len(submitted)
        submitted.append(job)
        c._callbacks[job_id] = callback
        c.runner._names[job_id] = job.name
        return job_id
    monkeypatch.setattr(c, '_submit', submit)
    monkeypatch.setattr(c.runner, 'cancel', cancelled.append)
    c.docking_started.connect(lambda: started.append(c.docking_running))
    c.docking_started.connect(c.cancel_docking)
    first_params = DockingParameters(seed=1)
    c.run_docking(first_params)
    run_dir = c.run_dir
    assert started == [True] and cancelled == [100]
    c.run_docking(DockingParameters(seed=2))
    assert len(submitted) == 1 and c.run_dir == run_dir and c.params is first_params
    # Cancellation frees the controller to start the next run.
    c.runner.job_cancelled.emit(100)
    assert not c.docking_running
    c.run_docking(DockingParameters(seed=3))
    assert len(submitted) == 2 and c.run_dir != run_dir
    c.runner.job_failed.emit(101, 'test failure')
    assert not c.docking_running


@pytest.mark.parametrize('reference', [[], [0.8]])
def test_all_export_formats_retain_poses_without_optional_rmsd(controller, tmp_path, reference):
    outcomes = sample_outcomes()
    outcomes[0].reference_rmsd = reference
    controller.outcomes = outcomes
    rows = results_rows(outcomes)
    assert len(rows) == 3
    assert rows[1]['RMSD to reference (A)'] == ''
    assert rows[0]['RMSD to reference (A)'] == (0.8 if reference else '')
    csv_path, excel_path, sdf_path = [tmp_path / f'results.{ext}' for ext in ('csv', 'xlsx', 'sdf')]
    controller.export_csv(csv_path)
    controller.export_excel(excel_path)
    controller.export_sdf(sdf_path)
    with csv_path.open(encoding='utf-8', newline='') as stream:
        assert len(list(csv.DictReader(stream))) == 3
    book = openpyxl.load_workbook(excel_path, read_only=True)
    assert book['Docking results'].max_row == 4
    book.close()
    molecules = list(Chem.SDMolSupplier(str(sdf_path)))
    assert len(molecules) == 3
    assert molecules[1].GetProp('RMSD_to_reference') == ''
    folder = tmp_path / 'ligand-output'
    folder.mkdir()
    _write_ligand_outputs(folder, outcomes[0])
    assert len(list(Chem.SDMolSupplier(str(folder / 'poses.sdf')))) == 2
