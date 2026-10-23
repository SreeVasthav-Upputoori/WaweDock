"""Native protein workspace stays useful when the real chemistry backend cannot load."""
import json
from pathlib import Path
import subprocess
import sys
from PySide6.QtCore import QTimer
from conftest import wait_until
import pytest
from test_protein_health_integration import window as base_window,js
from mds.core.capabilities import ChemistryStatus
from mds.project import Project
from mds.chem.workflow import DockingParameters


@pytest.fixture
def window(base_window):
    assert wait_until(lambda:base_window.controller.chemistry is not None,20)
    return base_window


def test_core_imports_and_controller_construction_do_not_import_rdkit(tmp_path):
    root=Path(__file__).resolve().parents[1]
    script=f'''
import sys
sys.path.insert(0,{str(root)!r})
from PySide6.QtWidgets import QApplication
from mds.gui.main_window import MainWindow
from mds.jobs import JobRunner
from mds.gui.controller import Controller
app=QApplication([])
c=Controller(JobRunner())
assert 'rdkit.Chem' not in sys.modules
'''
    result=subprocess.run([sys.executable,'-c',script],cwd=tmp_path,capture_output=True,text=True,timeout=20)
    assert result.returncode==0,result.stderr


def test_unavailable_chemistry_disables_tools_and_protein_health_still_runs(window):
    w=window;c=w.controller
    c._chemistry_checked(ChemistryStatus(False,'Policy blocked native dependency'))
    assert w.chemistry_notice.isVisible()
    for key in ('file.import_ligands','build.from_smiles','prepare.protein','prepare.ligands','docking.run','analyse.similarity'):
        assert not w.actions_[key].isEnabled()
    assert not w.docking_panel.prep_rec_btn.isEnabled()
    assert not w.docking_panel.prep_lig_btn.isEnabled()
    assert not w.docking_panel.run_btn.isEnabled()
    assert w.docking_panel.health_btn.isEnabled() and w.docking_panel.check_rec_btn.isEnabled()
    c.load_structure(Path('tests/data/3PTB.pdb'))
    w._show_protein_health();w.health_panel.check_btn.click()
    assert wait_until(lambda:w.health_panel.report is not None,20)
    assert w.health_panel.report.issues
    assert js(w,f'objects[{json.dumps(w._source_object)}].model.selectedAtoms({{}}).length')>1000
    c.inspect_receptor(w.docking_panel._receptor_options())
    assert wait_until(lambda:not c._receptor_check_job,20)
    assert w.docking_panel.rec_check.text().startswith('Check:')


def test_programmatic_docking_and_ligand_calls_fail_before_mutating_state(window):
    w=window;c=w.controller
    assert wait_until(lambda:c.chemistry is not None,10)
    c._chemistry_checked(ChemistryStatus(False,'Unavailable chemistry'))
    c.error.disconnect(w._show_error)
    errors=[];c.error.connect(lambda title,message:errors.append((title,message)))
    assert c.add_smiles('CCO','ethanol') is False
    c.import_ligand_file(Path('missing.sdf'))
    c.prepare_ligands(ionize='none',n_conformers=1,seed=42)
    c.run_docking(DockingParameters())
    assert len(errors)==4
    assert not c.entries and not c.prepared and not c.docking_running


def test_loading_and_editing_project_preserves_unloaded_ligands_and_runs(window,tmp_path):
    w=window;c=w.controller
    c._chemistry_checked(ChemistryStatus(False,'Policy blocked chemistry'))
    source=c.project.import_original(Path('tests/data/3PTB.pdb'))
    ligands=[{'name':'saved ligand','smiles':'CCO'}]
    c.project.state.update(source_structure=c.project.rel(source),ligands=ligands,runs=['runs/historical-run'],custom_metadata={'keep':True})
    c.project.save()
    c.load_project_state()
    assert c.source_pdb==source and not c.entries
    c.save_project()
    assert c.project.state['ligands']==ligands
    c.cleanup_structure(waters=True)
    assert c.source_pdb!=source
    c.save_project()
    saved=json.loads((c.project.root/'project.json').read_text())
    assert saved['ligands']==ligands
    assert saved['runs']==['runs/historical-run'] and saved['custom_metadata']=={'keep':True}
    c.undo_structure_edit();assert c.source_pdb==source


def test_successful_backend_probe_enables_existing_workflows(window):
    w=window
    # UI capability handling, without emulating chemistry calculations.
    w.controller._chemistry_checked(ChemistryStatus(True,version='test-status'))
    assert w.chemistry_notice.isHidden()
    assert w.actions_['file.import_ligands'].isEnabled()
    assert w.actions_['docking.run'].isEnabled()
    assert w.docking_panel.prep_lig_btn.isEnabled()


def test_probe_cancellation_finishes_capability_state(window):
    c=window.controller
    c._chemistry_job=123456
    c._callbacks[123456]=lambda result:None
    c._job_cancelled(123456)
    assert c.chemistry is not None and not c.chemistry.available
    assert 'cancelled' in c.chemistry.reason
    assert c._chemistry_job is None and 123456 not in c._callbacks
