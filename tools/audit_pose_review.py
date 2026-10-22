"""Review an isolated Vina audit run without requiring RDKit.

Usage: python tools/audit_pose_review.py test/pose-review-audit
The folder contains inputs/, pose-<rank>.pdbqt, and audit-job.json from one engine run.
"""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from PySide6.QtWidgets import QApplication,QVBoxLayout,QWidget
from mds.render.mol3d import Mol3DViewer
from mds.render.pose_review import PoseReviewSession
from mds.gui.panels.pose_review import PoseReviewBar


def audit(folder:Path):
    metadata=json.loads((folder/'audit-job.json').read_text())
    app=QApplication.instance() or QApplication([])
    window=QWidget(); window.setWindowTitle('WaweDock pose review audit'); window.resize(1000,800)
    layout=QVBoxLayout(window); bar=PoseReviewBar(); viewer=Mol3DViewer()
    layout.addWidget(bar); layout.addWidget(viewer,1); window.show()
    def wait(test,timeout=20):
        deadline=time.monotonic()+timeout
        while time.monotonic()<deadline:
            app.processEvents()
            if test(): return
            time.sleep(0.01)
        raise AssertionError('Viewer audit timed out')
    def js(code):
        result=[]; viewer.page().runJavaScript(code,result.append); wait(lambda:bool(result)); return result[0]
    try:
        wait(lambda:viewer._ready)
        viewer.add_object('protein',(folder/'inputs/receptor.pdb').read_text(),'pdb','protein')
        viewer.add_object('reference',(folder/'inputs/reference.pdb').read_text(),'pdb','ligand')
        session=PoseReviewSession(viewer)
        states={'protein':True,'reference':True}
        def visibility(name,on):
            expected=session.requested_visibility(name,on)
            if expected!=on: viewer.set_visible(name,expected)
            elif name in states: states[name]=on
        viewer.visibility_changed.connect(visibility)
        def apply():
            session.apply({'protein':states['protein']},{'reference':states['reference']},bar.reference.currentText(),
                          locked=bar.lock_reference.isChecked(),hide_protein=bar.hide_protein.isChecked())
        bar.settings_changed.connect(apply)
        bar.fit_requested.connect(session.fit_view)
        bar.set_references(['reference'])
        index=[0]; poses=metadata['poses']
        def show(i):
            index[0]=i; rank=poses[i]['rank']
            session.show_pose((folder/f'pose-{rank}.pdbqt').read_text(),'pdbqt',keep_camera=bar.keep_camera.isChecked())
            bar.set_position(i,len(poses))
        bar.step_requested.connect(lambda delta:show(max(0,min(len(poses)-1,index[0]+delta))))
        show(0)
        assert js("objects['pose: current'].model.selectedAtoms({}).length")>0
        bar.lock_reference.setChecked(True); bar.hide_protein.setChecked(True); bar.fit.click()
        js("setRepresentation('pose: current','ball-and-stick'); setLigandColor('pose: current','#cc33cc',true); true")
        camera=js('JSON.stringify(viewer.getView())')
        reference=js("JSON.stringify(objects.reference.model.selectedAtoms({}).map(a=>[a.x,a.y,a.z]))")
        screenshots=folder/'screenshots'; screenshots.mkdir(exist_ok=True)
        positions=[]
        for i in range(len(poses)):
            if i: bar.next.click()
            assert js('JSON.stringify(viewer.getView())')==camera
            assert js("JSON.stringify(objects.reference.model.selectedAtoms({}).map(a=>[a.x,a.y,a.z]))")==reference
            assert js('objects.reference.visible && !objects.protein.visible')
            assert js("objects['pose: current'].rep")=='ball-and-stick'
            assert js("objects['pose: current'].model.selectedAtoms({elem:'C'}).every(a=>a.style.stick.colorscheme.map.C==='#cc33cc')")
            positions.append(js("JSON.stringify(objects['pose: current'].model.selectedAtoms({}).map(a=>[a.x,a.y,a.z]))"))
            assert js('Object.keys(objects).length')==3
            app.processEvents(); viewer.grab().save(str(screenshots/f'pose-{i+1}.png'))
        assert len(set(positions))==len(poses)
        bar.previous.click(); assert index[0]==len(poses)-2
        assert js('JSON.stringify(viewer.getView())')==camera
        viewer.set_visible('reference',False); viewer.set_visible('protein',True)
        assert js('objects.reference.visible && !objects.protein.visible')
        bar.hide_protein.setChecked(False)
        assert js('objects.protein.visible')
        assert js('Object.keys(objects).length')==3
        session.clear_pose()
        assert js('objects.reference.visible && !!objects.protein && !objects["pose: current"]')
        result={'poses_reviewed':len(poses),'unique_pose_coordinates':len(set(positions)),
                'camera_preserved':True,'reference_coordinates_preserved':True,'protein_hide_restore':True,
                'locked_visibility_enforced':True,'appearance_preserved':True,'previous_next':True,
                'limitation':'Full MainWindow/controller launch is blocked by Windows App Control rejecting RDKit. Shared review session, controls and native viewer audited.'}
        (folder/'viewer-audit.json').write_text(json.dumps(result,indent=2))
        print(json.dumps(result,indent=2),flush=True)
    finally:
        window.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_folder',type=Path)
    audit(parser.parse_args().run_folder)
