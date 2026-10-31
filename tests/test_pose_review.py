"""Reference locking and reversible display masks are independent of chemistry imports."""
import json
import pytest
from conftest import wait_until
from mds.render.pose_review import PoseReviewSession
from mds.gui.panels.pose_review import PoseReviewBar


class Viewer:
    def __init__(self): self.calls=[]
    def set_visible(self,*args): self.calls.append(('visible',*args))
    def add_object(self,*args,**kwargs): self.calls.append(('add',args,kwargs))
    def remove_object(self,*args): self.calls.append(('remove',*args))
    def _call(self,*args): self.calls.append(args)


def test_reference_lock_and_protein_hide_restore_original_visibility():
    viewer=Viewer(); session=PoseReviewSession(viewer)
    proteins={'raw':False,'prepared':True}; references={'ref A':True,'ref B':True}
    session.apply(proteins,references,'ref A',locked=True,hide_protein=True)
    assert ('visible','prepared',False) in viewer.calls
    assert ('visible','ref B',False) in viewer.calls
    assert not session.requested_visibility('prepared',True)
    assert session.requested_visibility('ref A',False)
    session.apply({'raw':False,'prepared':False},{'ref A':True,'ref B':False},'ref A')
    assert ('visible','prepared',True) in viewer.calls
    assert ('visible','ref B',True) in viewer.calls
    assert ('visible','raw',True) not in viewer.calls
    assert not session.enforced_visibility


def test_pose_replacement_preserves_style_and_focus_excludes_protein():
    viewer=Viewer(); session=PoseReviewSession(viewer)
    session.apply({}, {'reference':True},'reference',locked=True)
    session.show_pose('pose one'); session.show_pose('pose two')
    adds=[call for call in viewer.calls if call[0]=='add']
    assert adds[0][1][0]==adds[1][1][0]==session.pose_name
    assert all(call[2]=={'preserve_view':True,'preserve_style':True} for call in adds)
    assert len([call for call in viewer.calls if call[0]=='setLigandColor'])==1
    session.fit_view()
    assert viewer.calls[-1]==('zoomToObjects',['reference',session.pose_name])
    session.clear_pose()
    assert viewer.calls[-1]==('remove',session.pose_name)
    assert session.reference=='reference'


def test_review_bar_navigation_and_lock_state(qapp):
    bar=PoseReviewBar(); steps=[]; bar.step_requested.connect(steps.append)
    try:
        assert not bar.lock_reference.isEnabled()
        assert not bar.next.isEnabled()
        bar.set_references(['reference'])
        bar.keep_camera.setChecked(False)
        bar.lock_reference.setChecked(True)
        assert bar.keep_camera.isChecked() and not bar.keep_camera.isEnabled()
        bar.set_position(0,6)
        assert not bar.previous.isEnabled() and bar.next.isEnabled()
        bar.next.click()
        bar.set_position(5,6)
        assert bar.previous.isEnabled() and not bar.next.isEnabled()
        bar.previous.click()
        assert steps==[1,-1]
        bar.set_references([])
        assert not bar.lock_reference.isChecked() and bar.keep_camera.isEnabled()
    finally:
        bar.close()


def test_live_pose_replacement_keeps_camera_and_custom_appearance(qapp):
    from mds.render.mol3d import Mol3DViewer
    viewer=Mol3DViewer(); viewer.resize(500,400); viewer.show()
    def js(expression):
        result=[]; viewer.page().runJavaScript(expression,result.append)
        assert wait_until(lambda:bool(result),10)
        return result[0]
    try:
        assert wait_until(lambda:viewer._ready,20)
        pdb='HETATM    1  C1  BEN A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n'
        viewer.add_object('ref',pdb,'pdb','ligand')
        session=PoseReviewSession(viewer)
        session.show_pose(pdb,'pdb')
        js("setRepresentation('pose: current','ball-and-stick'); setLigandColor('pose: current','#cc33cc',true); true")
        camera=js('JSON.stringify(viewer.getView())')
        session.show_pose(pdb.replace('0.000','1.000',1),'pdb')
        assert js('JSON.stringify(viewer.getView())')==camera
        assert js("objects['pose: current'].rep")=='ball-and-stick'
        assert js("objects['pose: current'].model.selectedAtoms({})[0].style.stick.colorscheme.map.C")=='#cc33cc'
        assert js('Object.keys(objects).length')==2
        assert js("objects.ref.model.selectedAtoms({})[0].x")==0
    finally:
        viewer.close()


def test_review_bar_hint_and_reference_selection(qapp):
    from mds.gui.panels.pose_review import NO_REFERENCE_HINT
    bar=PoseReviewBar(); bar.show()
    try:
        assert bar.hint.isVisible() and bar.hint.text()==NO_REFERENCE_HINT
        assert not bar.reference.isEnabled() and not bar.reference_label.isEnabled()
        assert not bar.select_reference('STI (crystal)')
        bar.set_references(['STI (crystal)','my ligand'])
        assert not bar.hint.isVisible() and bar.reference_label.isEnabled()
        assert bar.select_reference('my ligand')
        assert bar.reference.currentText()=='my ligand' and bar.lock_reference.isChecked()
        assert not bar.keep_camera.isEnabled()
    finally:
        bar.close()


def test_review_bar_position_detail_and_shortcut_gating(qapp):
    bar=PoseReviewBar(); steps=[]; bar.step_requested.connect(steps.append)
    try:
        bar.set_position(2,6,'caffeine #3 · -5.08 kcal/mol')
        assert bar.position.text()=='Pose 3 / 6 — caffeine #3 · -5.08 kcal/mol'
        bar.set_position(-1,6,'ignored')
        assert bar.position.text()=='Pose 0 / 6'
        bar.set_position(1,6,'caffeine #2 · -5.12 kcal/mol')
        bar._step_shortcut(bar.next)          # hidden bar: keys stay inert
        assert steps==[]
        bar.show()
        bar._step_shortcut(bar.next)
        bar.set_position(5,6)
        bar._step_shortcut(bar.next)          # disabled at the last pose
        assert steps==[1]
    finally:
        bar.close()


def test_explorer_offers_pose_reference_for_ligands_only(qapp):
    from mds.gui.panels.project_explorer import ProjectExplorer
    ex=ProjectExplorer()
    try:
        ex.upsert_object('Ligands','STI (crystal)',{},kind='ligand')
        ex.upsert_object('Proteins','1XBB',{},kind='protein')
        requested=[]; ex.reference_requested.connect(requested.append)

        def reference_actions(name):
            menu=ex.build_context_menu(ex.find(name))
            return [a for a in menu.actions() if a.text()=='Use as pose reference']

        assert not reference_actions('1XBB')
        assert not [a for a in ex.build_context_menu(None).actions()
                    if a.text()=='Use as pose reference']
        actions=reference_actions('STI (crystal)')
        assert len(actions)==1
        actions[0].trigger()
        assert requested==['STI (crystal)']
    finally:
        ex.close()
