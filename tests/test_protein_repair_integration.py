"""Real native repair workflow: proposals, worker process, review, Undo and colours."""
import json
from pathlib import Path
from PySide6.QtCore import QTimer, Qt
from conftest import wait_until
from test_protein_health_integration import window,js
from test_protein_repair import damaged
from test_protein_loops import loop_source


def test_repair_preview_accept_undo_redo_and_viewer_styles(window,damaged):
    w=window; c=w.controller; panel=w.health_panel
    c.load_structure(damaged)
    original=c.source_pdb; original_bytes=original.read_bytes()
    old_name=w._source_object
    js(w,f"setBackgroundColor('#102030');setRepresentation({json.dumps(old_name)},'ribbon');true")
    w._show_protein_health()
    panel.repair_btn.click()
    assert wait_until(lambda:panel.repair_plan is not None,15)
    assert c.source_pdb==original
    assert panel.build_btn.isEnabled()
    ticks=[]; timer=QTimer(); timer.setInterval(20);timer.timeout.connect(lambda:ticks.append(1));timer.start()
    panel.build_btn.click()
    assert not panel.accept_btn.isEnabled()
    try: assert wait_until(lambda:panel.repair_result is not None,120),panel.details.toPlainText()
    finally: timer.stop()
    result=panel.repair_result
    assert ticks,'GUI timer must keep running during native repair'
    assert result.can_accept,result.blockers
    assert c.source_pdb==original and original.read_bytes()==original_bytes
    assert not c._undo_structures
    assert panel.accept_btn.isEnabled()
    preview=json.dumps(w._repair_preview_name)
    assert js(w,f'!!objects[{preview}]')
    assert not js(w,f'!!objects[{preview}].editable')
    index=result.added[0]['atom_index']
    assert js(w,f"objects[{preview}].model.selectedAtoms({{}})[{index}].style.stick.color")=='#00d6d6'
    # Original anchors are exposed with real parsed bonds even in ribbon view.
    assert js(w,f"!!objects[{preview}].model.selectedAtoms({{resi:32,atom:'CB'}})[0].style.stick")
    assert js(w,f"!!objects[{preview}].model.selectedAtoms({{resi:32,atom:'CA'}})[0].style.stick")
    assert '0 loop model(s) built' in panel.candidate_summary.text()
    assert 'require sequence evidence' in panel.repair_limits.text()
    from dataclasses import replace
    panel._candidate_ready(replace(result,blockers=['Severe reconstructed-atom overlap: example, 1.42 A.']))
    assert '1.42 A.' in panel.candidate_summary.text()
    assert not panel.accept_btn.isEnabled()
    panel._candidate_ready(result)
    panel.comparison.setCurrentIndex(0)
    assert not js(w,f'!!objects[{preview}]')
    assert w.project_explorer.is_visible(old_name)
    panel.comparison.setCurrentIndex(1)
    assert js(w,f'!!objects[{preview}]')
    w.dock_health.hide()
    assert not js(w,f'!!objects[{preview}]')
    assert w.project_explorer.is_visible(old_name)
    assert panel.comparison.currentIndex()==0
    w._show_protein_health()
    panel.accept_btn.click()
    assert c.source_pdb==Path(result.candidate)
    assert wait_until(lambda:panel.report is not None,15)
    assert c.receptor is None
    assert panel.undo_btn.isEnabled()
    current=json.dumps(w._source_object)
    assert js(w,'backgroundColor')=='#102030'
    assert js(w,f'objects[{current}].rep')=='ribbon'
    assert js(w,f'objects[{current}].reconstructedIndices.length')==1
    panel.undo_btn.click()
    assert c.source_pdb==original
    c.redo_structure_edit()
    assert c.source_pdb==Path(result.candidate)
    current=json.dumps(w._source_object)
    assert js(w,f'objects[{current}].reconstructedIndices.length')==1
    assert original.read_bytes()==original_bytes


def test_manual_repair_uses_existing_exact_viewer_selection(window,damaged):
    w=window; c=w.controller
    c.load_structure(damaged)
    # Atom selection expands the manual scope to this exact residue, not the chain.
    name=json.dumps(w._source_object)
    js(w,f"setSelectionMode('atom');pickSelection({name},objects[{name}].model.selectedAtoms({{resi:32,atom:'CB'}})[0],{{}});true")
    assert wait_until(lambda:bool(c.selection),5)
    w._repair_health_selection()
    panel=w.health_panel
    assert wait_until(lambda:panel.repair_plan is not None,15)
    assert panel.repair_scope.currentIndex()==1
    assert [ch.kind for ch in panel.repair_plan.changes]==['heavy_atoms']
    assert '32' in panel.repair_plan.changes[0].label


def test_cancel_native_repair_leaves_original_active(window,damaged):
    w=window; c=w.controller; p=w.health_panel
    c.load_structure(damaged); original=c.source_pdb
    p.repair_btn.click()
    assert wait_until(lambda:p.repair_plan is not None,15)
    p.build_btn.click()
    p.cancel_btn.click()
    assert wait_until(lambda:c.runner.active_count()==0,15)
    assert c.source_pdb==original
    assert not c._undo_structures
    assert p.repair_result is None
    assert not p.accept_btn.isEnabled()

def test_repair_choices_stay_readable_and_conformers_are_mutually_exclusive(window,tmp_path,qapp):
    from test_protein_health import structure_file,GLY
    from PySide6.QtWidgets import QAbstractButton
    w=window; panel=w.health_panel
    source=structure_file(tmp_path,[('A',1,'ALA',GLY+[('CB',1.5,-1.3,0,'A'),('CB',1.5,-1.3,1,'B')])])
    w.controller.load_structure(source)
    panel.repair_btn.click()
    assert wait_until(lambda:panel.repair_plan is not None,15)
    panel.proposed_changes.item(1).setCheckState(Qt.CheckState.Checked)
    assert panel.proposed_changes.item(0).checkState()==Qt.CheckState.Unchecked
    # Detach only for this layout check, keeping the real controller/plan.
    panel.setParent(None)
    panel.resize(240,650)
    panel.show()
    for _ in range(3): qapp.processEvents()
    try:
        assert panel.scroll.horizontalScrollBar().maximum()==0
        viewport=panel.scroll.viewport()
        for button in panel.findChildren(QAbstractButton):
            if not button.isVisible(): continue
            x=button.mapTo(viewport,button.rect().topLeft()).x()
            assert x>=0 and x+button.width()<=viewport.width()
    finally:
        panel.close()
        w.dock_health.setWidget(panel)


def test_loop_job_recovers_project_sequence_previews_and_undoes(window,tmp_path,loop_source):
    import gemmi
    from test_protein_loops import _loops
    w=window;c=w.controller;p=w.health_panel
    source,full,_=loop_source
    c.load_structure(full)  # immutable original with sequence evidence
    stripped=tmp_path/'previously-cleaned.pdb'
    gemmi.read_structure(str(source)).write_minimal_pdb(str(stripped))
    c.load_structure(stripped);original=c.source_pdb;original_bytes=original.read_bytes()
    w._show_protein_health();p.repair_btn.click()
    assert wait_until(lambda:p.repair_plan is not None,20)
    loop=_loops(p.repair_plan)[0]
    assert Path(loop.loop['evidence']).parent==c.project.root/'originals'
    assert Path(loop.loop['evidence'])!=original
    item=next(p.proposed_changes.item(i) for i in range(p.proposed_changes.count())
              if p.proposed_changes.item(i).data(Qt.ItemDataRole.UserRole)==loop.id)
    assert item.checkState()==Qt.CheckState.Unchecked
    item.setCheckState(Qt.CheckState.Checked)
    assert all(p.proposed_changes.item(i).checkState()==Qt.CheckState.Unchecked
               for i,ch in enumerate(p.repair_plan.changes) if ch.kind!='loop')
    ticks=[];timer=QTimer();timer.setInterval(20);timer.timeout.connect(lambda:ticks.append(1));timer.start()
    p.build_btn.click()
    try: assert wait_until(lambda:p.repair_result is not None,90),p.details.toPlainText()
    finally: timer.stop()
    result=p.repair_result
    assert result.can_accept,result.blockers
    assert len(result.loops)==1 and ticks
    assert c.source_pdb==original and original.read_bytes()==original_bytes
    assert '1 loop model(s) built' in p.candidate_summary.text()
    preview=json.dumps(w._repair_preview_name)
    assert js(w,f'objects[{preview}].reconstructedIndices.length')==len(result.added)
    assert not js(w,f'objects[{preview}].editable')
    p.accept_btn.click()
    assert c.source_pdb==Path(result.candidate)
    assert wait_until(lambda:p.report is not None,15)
    c.undo_structure_edit();assert c.source_pdb==original
    c.redo_structure_edit();assert c.source_pdb==Path(result.candidate)
    assert original.read_bytes()==original_bytes
