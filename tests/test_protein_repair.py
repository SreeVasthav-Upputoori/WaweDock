"""Real repair backends, immutable candidates and scientific acceptance gates."""
from pathlib import Path
import json
import gemmi
import pytest
from mds.chem.repair import (propose_repairs, build_repairs, fingerprint,
                             export_repaired_structure, _validate_additions)
from test_protein_health import structure_file, GLY, translated

DATA=Path('tests/data/3PTB.pdb')


@pytest.fixture
def damaged(tmp_path):
    st=gemmi.read_structure(str(DATA))
    residue=next(r for r in st[0][0] if r.name=='SER')
    residue.remove_atom('OG','*')
    path=tmp_path/'missing.pdb'
    options=gemmi.PdbWriteOptions()
    options.preserve_serial=True
    options.conect_records=True
    st.write_pdb(str(path),options)
    return path


def snapshot(path):
    st=gemmi.read_structure(str(path))
    return {(mi,ci,ri,a.name,a.altloc):(a.element.name,tuple(a.pos),a.occ,a.b_iso,a.serial)
            for mi,m in enumerate(st) for ci,c in enumerate(m) for ri,r in enumerate(c) for a in r}


def test_real_atom_repair_preserves_existing_coordinates_and_conect(damaged,tmp_path):
    before=snapshot(damaged); digest=fingerprint(damaged)
    plan=propose_repairs(damaged)
    chosen=[c.id for c in plan.changes if c.kind=='heavy_atoms']
    assert len(chosen)==1
    result=build_repairs(plan,chosen,tmp_path/'candidate')
    assert result.can_accept, result.blockers
    assert len(result.added)==1 and result.added[0]['atom']=='OG'
    after=snapshot(result.candidate)
    assert all(after[key]==value for key,value in before.items())
    assert fingerprint(damaged)==digest
    assert gemmi.read_structure(str(damaged)).conect_map==gemmi.read_structure(result.candidate).conect_map
    assert not any(i.code=='missing_atoms' for i in result.after.issues)
    report=json.loads(Path(result.candidate+'.repair.json').read_text())
    assert report['status']=='Review required'
    assert report['added'][0]['atom']=='OG'
    assert 'computational' in report['coordinate_provenance']
    export_repaired_structure(result,tmp_path/'export.cif')
    assert '_wavedock_repair_atom.' in (tmp_path/'export.cif').read_text()
    assert any(c.name.startswith('wavedock_bond_') for c in gemmi.read_structure(str(tmp_path/'export.cif')).connections)
    assert (tmp_path/'export.cif.repair.json').exists()


def test_real_hydrogen_addition_preserves_nonprotein_components(tmp_path):
    plan=propose_repairs(DATA)
    result=build_repairs(plan,['hydrogens'],tmp_path/'hydrogen')
    assert result.can_accept or result.blockers
    assert len(result.added)>1000
    assert all(a['atom'].startswith('H') for a in result.added)
    disulfides={r.seqid.num for ch in gemmi.read_structure(str(DATA))[0] for r in ch if r.name=='CYS'}
    assert not any(a['resi'] in disulfides and a['atom']=='HG' for a in result.added)
    after=snapshot(result.candidate)
    assert all(after[key]==value for key,value in snapshot(DATA).items())
    assert not any(i.code=='hydrogens' for i in result.after.issues)
    assert 'WAVEDOCK REPAIR' in Path(result.candidate).read_text()
    repeat=propose_repairs(result.candidate)
    assert not any(c.kind=='hydrogens' for c in repeat.changes)


def test_manual_scope_and_stale_preview(damaged,tmp_path):
    assert not any(c.kind=='heavy_atoms' for c in propose_repairs(damaged,selected_residues={(0,0,0)}).changes)
    plan=propose_repairs(damaged)
    damaged.write_text(damaged.read_text()+'REMARK changed\n')
    with pytest.raises(ValueError,match='stale'):
        build_repairs(plan,[c.id for c in plan.changes if c.kind=='heavy_atoms'],tmp_path/'candidate')


def test_conformer_choice_requires_one_choice_and_preserves_others(tmp_path):
    # Shared backbone and explicit, separate CB conformers.
    atoms=GLY+[('CB',1.5,-1.3,0,'A'),('CB',1.5,-1.3,1,'B')]
    path=structure_file(tmp_path,[('A',1,'ALA',atoms),('Z',1,'GLY',translated(GLY,20))])
    plan=propose_repairs(path)
    choices=[c.id for c in plan.changes if c.kind=='altloc']
    assert len(choices)==2
    with pytest.raises(ValueError,match='one conformer'):
        build_repairs(plan,choices,tmp_path/'invalid')
    result=build_repairs(plan,[choices[0]],tmp_path/'good')
    assert result.can_accept,result.blockers
    assert len(result.removed)==1
    res=gemmi.read_structure(result.candidate)[0][0][0]
    assert all(a.altloc=='\0' for a in res)
    assert [r for r in snapshot(path) if r[1]==1]==[r for r in snapshot(result.candidate) if r[1]==1]


def test_missing_backbone_anchors_and_unresolved_gaps_are_not_auto_rebuilt(tmp_path):
    path=structure_file(tmp_path,[('A',1,'ALA',[a for a in GLY if a[0]!='N']),('A',9,'GLY',translated(GLY,30))])
    plan=propose_repairs(path)
    assert not any(c.kind in ('heavy_atoms','hydrogens') for c in plan.changes)
    assert any('advanced modelling' in n for n in plan.notes)


def test_reconstructed_geometry_gate_rejects_bad_bonds_and_contacts():
    st=gemmi.Structure(); m=gemmi.Model('1'); c=gemmi.Chain('A'); r=gemmi.Residue()
    r.name='ALA';r.seqid=gemmi.SeqId(1,' ')
    for name,xyz in [('CA',(0,0,0)),('CB',(8,0,0)),('O',(8.1,0,0))]:
        a=gemmi.Atom();a.name=name;a.element=gemmi.Element(name[0]);a.pos=gemmi.Position(*xyz);r.add_atom(a)
    c.add_residue(r);m.add_chain(c);st.add_model(m)
    blockers=_validate_additions(st,{(0,0,'CB')},[((0,0,'CA'),(0,0,'CB'))])
    assert any('bond' in b for b in blockers)
    assert any('overlap' in b for b in blockers)


def test_no_unapproved_export_or_cancellation(damaged,tmp_path):
    plan=propose_repairs(damaged)
    with pytest.raises(ValueError,match='Select valid'):
        build_repairs(plan,[],tmp_path/'none')
    def cancel(): raise RuntimeError('cancelled')
    with pytest.raises(RuntimeError,match='cancelled'):
        build_repairs(plan,[c.id for c in plan.changes if c.kind=='heavy_atoms'],tmp_path/'none',check_cancelled=cancel)

def test_mmcif_repairs_and_provenance_survive_export_and_reimport(damaged,tmp_path):
    from mds.chem.repair import reconstructed_atom_indices
    cif=tmp_path/'source.cif'
    gemmi.read_structure(str(damaged)).make_mmcif_document().write_file(str(cif))
    source_doc=gemmi.cif.read_file(str(cif))
    source_doc.sole_block().set_mmcif_category('_chem_comp_bond.',{
        'comp_id':['BEN'],'atom_id_1':['N1'],'atom_id_2':['C1'],'value_order':['SING']})
    source_doc.sole_block().set_pair('_audit.creation_method','test-source-metadata')
    source_doc.write_file(str(cif))
    plan=propose_repairs(cif)
    result=build_repairs(plan,[c.id for c in plan.changes if c.kind=='heavy_atoms'],tmp_path/'candidate-cif')
    assert result.can_accept,result.blockers
    saved_block=gemmi.cif.read_file(result.candidate).sole_block()
    assert saved_block.get_mmcif_category('_chem_comp_bond.')==source_doc.sole_block().get_mmcif_category('_chem_comp_bond.')
    assert saved_block.find_value('_audit.creation_method')=='test-source-metadata'
    assert reconstructed_atom_indices(Path(result.candidate))==[result.added[0]['atom_index']]
    export_repaired_structure(result,tmp_path/'reimport.pdb')
    assert gemmi.read_structure(str(tmp_path/'reimport.pdb')).conect_map
    # File annotations remain meaningful even when an external importer omits the sidecar.
    (tmp_path/'reimport.pdb.repair.json').unlink()
    assert reconstructed_atom_indices(tmp_path/'reimport.pdb')==[result.added[0]['atom_index']]


def test_later_hydrogen_revision_keeps_prior_reconstructed_atom_provenance(damaged,tmp_path):
    from mds.chem.repair import reconstructed_atom_indices
    plan=propose_repairs(damaged)
    first=build_repairs(plan,[c.id for c in plan.changes if c.kind=='heavy_atoms'],tmp_path/'first')
    second_plan=propose_repairs(Path(first.candidate))
    second=build_repairs(second_plan,['hydrogens'],tmp_path/'second')
    assert len(second.inherited)==1 and second.inherited[0]['atom']=='OG'
    assert len(reconstructed_atom_indices(Path(second.candidate)))==len(second.added)+1
    assert all(a['atom']!='HG' for a in second.added if a['resn']=='CYS')

def test_geometry_gate_rejects_linear_reconstructed_sidechain_angle():
    st=gemmi.Structure(); m=gemmi.Model('1'); c=gemmi.Chain('A'); r=gemmi.Residue()
    r.name='LEU';r.seqid=gemmi.SeqId(1,' ')
    for i,name in enumerate(('CA','CB','CG')):
        a=gemmi.Atom();a.name=name;a.element=gemmi.Element('C');a.pos=gemmi.Position(i*1.53,0,0);r.add_atom(a)
    c.add_residue(r);m.add_chain(c);st.add_model(m)
    assert any('angle' in text for text in _validate_additions(st,{(0,0,'CG')},
        [((0,0,'CA'),(0,0,'CB')),((0,0,'CB'),(0,0,'CG'))]))

def test_single_conformer_label_resolution_is_a_real_reversible_change(tmp_path):
    atoms=[(name,x,y,z,'A') if name=='CA' else (name,x,y,z) for name,x,y,z in GLY]
    path=structure_file(tmp_path,[('A',1,'GLY',atoms)])
    plan=propose_repairs(path)
    result=build_repairs(plan,[c.id for c in plan.changes if c.kind=='altloc'],tmp_path/'single-alt')
    assert result.can_accept,result.blockers
    assert not result.added and not result.removed
    assert len(result.normalised)==1
    assert not any(i.code=='alternate_conformations' for i in result.after.issues)

def test_modified_candidate_cannot_be_activated(qapp,damaged,tmp_path):
    from mds.gui.controller import Controller
    from mds.jobs import JobRunner
    from mds.project import Project
    controller=Controller(JobRunner())
    controller.set_project(Project.create(tmp_path/'project'))
    controller.load_structure(damaged)
    original=controller.source_pdb
    plan=propose_repairs(original)
    result=build_repairs(plan,[c.id for c in plan.changes if c.kind=='heavy_atoms'],tmp_path/'stage')
    candidate=Path(result.candidate)
    candidate.write_text(candidate.read_text()+'REMARK externally modified\n')
    errors=[]
    controller.error.connect(lambda title,message:errors.append(message))
    assert not controller.accept_protein_repair(result)
    assert controller.source_pdb==original and not controller._undo_structures
    assert errors and 'changed' in errors[-1]
