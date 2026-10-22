"""Real application/job/viewer integration for the health diagnostics increment."""
from pathlib import Path
import json

import pytest
from conftest import wait_until
from mds.gui.main_window import MainWindow
from mds.project import Project


@pytest.fixture
def window(qapp,settings,tmp_path):
    w=MainWindow(settings,use_viewer=True)
    w.controller.set_project(Project.create(tmp_path/'project'))
    w.show()
    assert wait_until(lambda:w.viewer._ready,20)
    yield w
    w.job_runner.cancel_all()
    w.job_runner.wait_for_done(10000)
    w.close()


def js(w,expression):
    results=[]
    w.viewer.page().runJavaScript(expression,results.append)
    assert wait_until(lambda:bool(results),10)
    return results[0]


def test_health_job_navigation_preserves_appearance_and_components(window):
    w=window
    w.controller.load_structure(Path('tests/data/3PTB.pdb'))
    source=w.controller.source_pdb
    before=source.read_bytes()
    assert w.dock_health.isHidden()
    w.docking_panel.health_btn.click()
    assert w.dock_health.isVisible()
    name=json.dumps(w._source_object)
    js(w,"setBackgroundColor('#182838');setRepresentation("+name+",'sticks');true")
    original_styles=js(w,"JSON.stringify(objects["+name+"].model.selectedAtoms({}).map(a=>a.style))")
    original_atoms=js(w,"objects["+name+"].model.selectedAtoms({}).length")
    w.health_panel.check_btn.click()
    assert not w.health_panel.check_btn.isEnabled()
    assert wait_until(lambda:w.health_panel.report is not None,10)
    report=w.health_panel.report
    assert source.read_bytes()==before
    assert w.controller.source_pdb==source
    assert w.controller.receptor is None
    assert not w.controller._undo_structures
    assert w.health_panel.export_btn.isEnabled()
    assert js(w,'healthShapes.length')>0
    row=next(i for i,issue in enumerate(report.issues) if issue.atoms)
    w.health_panel.issues.setCurrentRow(row)
    assert wait_until(lambda:bool(w.controller.selection),5)
    assert set(a['atom_index'] for a in w.controller.selection[w._source_object])==set(a['atom_index'] for a in report.issues[row].atoms)
    assert js(w,'backgroundColor')=='#182838'
    assert js(w,'objects['+name+'].rep')=='sticks'
    w.viewer.clear_selection()
    w.health_panel.highlights.setChecked(False)
    assert js(w,'healthShapes.length')==0
    assert js(w,'objects['+name+'].model.selectedAtoms({}).length')==original_atoms
    assert js(w,'JSON.stringify(objects['+name+'].model.selectedAtoms({}).map(a=>a.style))')==original_styles
    w.health_panel.highlights.setChecked(True)
    assert js(w,'healthShapes.length')>0
    w.viewer.set_visible(w._source_object,False)
    assert js(w,'healthShapes.length')==0
    w.viewer.set_visible(w._source_object,True)
    assert js(w,'healthShapes.length')>0


def test_source_change_discards_in_flight_health_result(window,qapp):
    w=window
    w.controller.load_structure(Path('tests/data/3PTB.pdb'))
    w.health_panel.check()
    w.controller.load_structure(Path('tests/data/3PTB.pdb'))
    w.job_runner.wait_for_done(10000)
    for _ in range(5): qapp.processEvents()
    assert w.health_panel.report is None
    assert w.health_panel.issues.count()==0
    assert w.health_panel.check_btn.isEnabled()
    assert not w.health_panel.export_btn.isEnabled()
    assert js(w,'healthShapes.length')==0
