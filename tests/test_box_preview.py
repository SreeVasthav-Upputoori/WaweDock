"""Geometry checks for a positioned reference, distinct from ligand size advisories."""
from pathlib import Path
import json
import numpy as np
import pytest
from conftest import wait_until
from mds.chem.box import DockingBox, box_from_coords, receptor_atom_coords, validate_box


def test_positioned_ligand_fit_reports_counts_and_clearance():
    box=DockingBox((10,20,30),(10,10,10))
    inside=np.array([[10,20,30],[11,21,31]])
    report=validate_box(box,ligand_coords=inside,ligand_label='BEN A:1')
    assert not report.errors and not report.warnings
    assert 'all 2 heavy atoms inside' in report.render()
    assert '4.00 A' in report.render()
    outside=np.array([[10,20,30],[16,20,30],[10,13,30]])
    report=validate_box(box,ligand_coords=outside)
    assert any('2/3 heavy atom(s) outside' in w and '2.00 A' in w for w in report.warnings)
    assert not report.errors  # This pose is advisory, not a restriction on docking search.


def test_box_face_and_near_face_are_not_misclassified_as_outside():
    box=DockingBox((0,0,0),(10,10,10))
    report=validate_box(box,ligand_coords=np.array([[5,0,0],[-5,0,0],[5+1e-8,0,0]]))
    assert 'all 3 heavy atoms inside' in report.render()
    assert any('less than 2 A' in w for w in report.warnings)
    assert not any('outside' in w for w in report.warnings)


@pytest.mark.parametrize('box',[DockingBox((float('nan'),0,0),(10,10,10)),
                                DockingBox((0,0,0),(float('inf'),10,10)),
                                DockingBox((0,0,0),(0,10,10))])
def test_invalid_box_rejected(box):
    assert validate_box(box).errors


def test_unpositioned_ligand_is_not_compared_with_box_origin():
    report=validate_box(DockingBox((100,100,100),(10,10,10)),ligand_extent=12)
    assert not report.errors
    assert 'No positioned reference ligand' in report.render()
    assert any('longest dimension' in w for w in report.warnings)
    assert not any('atom(s) outside' in w for w in report.warnings)


def test_empty_receptor_and_invalid_reference_have_useful_messages(tmp_path):
    path=tmp_path/'empty.pdbqt'; path.write_text('REMARK no atoms\n')
    assert receptor_atom_coords(path).shape==(0,3)
    report=validate_box(DockingBox((0,0,0),(10,10,10)),path,ligand_coords=np.array([[np.nan,0,0]]))
    assert any('no receptor atoms' in e for e in report.errors)
    assert any('cannot assess fit' in w for w in report.warnings)
    with pytest.raises(ValueError): box_from_coords(np.array([[np.nan,0,0]]))


def test_viewer_box_status_colors_preserve_camera(qapp):
    # The standalone viewer needs no RDKit, allowing real rendering validation.
    from mds.render.mol3d import Mol3DViewer
    viewer=Mol3DViewer()
    viewer.resize(500,400); viewer.show()
    try:
        assert wait_until(lambda:viewer._ready,20)
        for status,color in [('ok','#2ecc71'),('warning','#f39c12'),('error','#e74c3c')]:
            results=[]
            script="""JSON.stringify((() => {
              const before=viewer.getView(), add=viewer.addBox, colors=[];
              viewer.addBox=function(spec) { colors.push(spec.color); return add.call(this,spec); };
              try { showBox([1,2,3],[10,12,14],STATUS); }
              finally { viewer.addBox=add; }
              return {colors,count:boxShapes.length,same:JSON.stringify(before)===JSON.stringify(viewer.getView())};
            })())""".replace('STATUS',json.dumps(status))
            viewer.page().runJavaScript(script,results.append)
            assert wait_until(lambda:bool(results),10)
            result=json.loads(results[0])
            assert result=={'colors':[color,color],'count':2,'same':True}
    finally:
        viewer.close()


def test_controller_tracks_positioned_reference_and_manual_box(qapp,tmp_path):
    from mds.gui.controller import Controller
    from mds.jobs import JobRunner
    import gemmi
    ligand=tmp_path/'BEN.pdb'
    residue=gemmi.Residue(); residue.name='BEN'; residue.seqid=gemmi.SeqId(1,' '); residue.het_flag='H'
    atom=gemmi.Atom(); atom.name='C1'; atom.element=gemmi.Element('C'); atom.pos=gemmi.Position(10,20,30)
    residue.add_atom(atom)
    structure=gemmi.Structure(); model=gemmi.Model('1'); chain=gemmi.Chain('A')
    chain.add_residue(residue); model.add_chain(chain); structure.add_model(model); structure.write_pdb(str(ligand))
    controller=Controller(JobRunner())
    controller.set_box_reference(ligand)
    controller.set_box(DockingBox((10,20,30),(10,10,10)))
    assert 'all 1 heavy atoms inside' in controller.validate_box().render()
    controller.set_box(DockingBox((100,20,30),(10,10,10)))
    assert any('1/1 heavy atom(s) outside' in w for w in controller.validate_box().warnings)
    controller.set_box_reference(None)
    assert 'No positioned reference ligand' in controller.validate_box().render()
