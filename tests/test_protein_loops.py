"""Real sequence-backed loop reconstruction and conservative eligibility checks."""
from pathlib import Path
import json
import gemmi
import pytest
from mds.chem.edit import atom_records,edit_atoms
from mds.chem.repair import propose_repairs,build_repairs,fingerprint,export_repaired_structure
from mds.chem.loops import propose_loops,validate_loops,_positions
from test_protein_health import structure_file,GLY,translated


@pytest.fixture
def loop_source(tmp_path):
    """Remove one observed residue from a real trypsin backbone (known sequence)."""
    native=gemmi.read_structure('tests/data/3PTB.pdb')
    # Five contiguous experimentally observed residues with complete anchors.
    residues=[r.clone() for r in list(native[0][0])[15:20]]
    st=gemmi.Structure();m=gemmi.Model('1');ch=gemmi.Chain('A')
    for i,r in enumerate(residues):
        r.seqid=gemmi.SeqId(i+1,' ');r.label_seq=i+1;r.subchain='A';ch.add_residue(r)
    m.add_chain(ch);st.add_model(m);st.setup_entities()
    st.entities[0].full_sequence=[r.name for r in residues]
    # Preserve the real non-protein ligand/water/metal content as fixed context.
    for r in native[0][0]:
        if r.het_flag=='H': st[0][0].add_residue(r.clone())
    full=tmp_path/'native.pdb';st.write_pdb(str(full))
    target=tmp_path/'gap.pdb'
    edit_atoms(full,target,[a for a in atom_records(st) if a['resi']==3 and a['resn']==residues[2].name])
    return target,full,residues[2].name


def _loops(plan): return [c for c in plan.changes if c.kind=='loop']


def test_sequence_evidence_and_manual_loop_scope(loop_source):
    source,full,name=loop_source
    plan=propose_repairs(source)
    loops=_loops(plan)
    assert len(loops)==1 and not loops[0].recommended
    assert loops[0].loop['numbers']==[3] and loops[0].loop['names']==[name]
    assert not _loops(propose_repairs(source,selected_residues={(0,0,0)}))
    assert len(_loops(propose_repairs(source,selected_residues={(0,0,1)})))==1
    st=gemmi.read_structure(str(source));st.write_minimal_pdb(str(source))
    assert not _loops(propose_repairs(source))
    recovered=_loops(propose_repairs(source,sequence_references=[full]))
    assert len(recovered)==1 and recovered[0].loop['evidence']==str(full)


def test_no_sequence_ambiguous_mapping_termini_or_explicit_ter_are_not_modelled(tmp_path):
    path=structure_file(tmp_path,[('A',1,'GLY',GLY),('A',3,'GLY',translated(GLY,8))])
    assert not _loops(propose_repairs(path))
    assert _positions(gemmi.read_structure(str(path))[0][0],['GLY']*5) is None
    original=gemmi.read_structure(str(path));original.setup_entities();original.entities[0].full_sequence=['GLY','ALA','GLY']
    original.write_pdb(str(path));st=gemmi.read_structure(str(path))
    assert len(propose_loops(path,st,[],())[0])==1
    assert not propose_loops(path,st,[(0,1)],())[0]
    st.entities[0].full_sequence=['SER','GLY','ALA','GLY','LYS'];st.write_pdb(str(path))
    assert len(_loops(propose_repairs(path)))==1  # internal gap only, no terminal extension


@pytest.mark.parametrize("suffix",[".pdb",".cif"])
def test_real_loop_model_preserves_original_atoms_and_exports(loop_source,tmp_path,suffix):
    source,_,name=loop_source
    if suffix==".cif":
        converted=source.with_suffix(suffix)
        gemmi.read_structure(str(source)).make_mmcif_document().write_file(str(converted))
        source=converted
    plan=propose_repairs(source);changes=_loops(plan);original=source.read_bytes()
    result=build_repairs(plan,[c.id for c in changes],tmp_path/'candidate')
    assert result.loops and result.modelling['trials']
    assert result.can_accept,result.blockers
    assert {a['resn'] for a in result.added}=={name}
    assert {a['resi'] for a in result.added}=={3}
    before=gemmi.read_structure(str(source));after=gemmi.read_structure(result.candidate)
    def existing(st):
        return {(mi,ch.name,r.seqid.num,r.seqid.icode,r.name,a.name,a.altloc):(tuple(a.pos),a.serial,a.occ,a.b_iso)
                for mi,m in enumerate(st) for ch in m for r in ch for a in r if r.seqid.num!=3 or r.name!=name}
    assert existing(before)==existing(after)
    assert source.read_bytes()==original
    assert not any(i.code=='backbone_connectivity' for i in result.after.issues)
    assert len(result.modelling['geometry'])==3  # two peptide junctions plus one phi/psi
    assert all(not t.get('error') for t in result.modelling['trials'])
    report=json.loads(Path(result.candidate+'.repair.json').read_text())
    assert report['parameters']['loop_reconstruction']
    assert 'experimental' in result.warnings[0]
    export_repaired_structure(result,tmp_path/'loop-export.cif')
    assert '_wavedock_repair_atom.' in (tmp_path/'loop-export.cif').read_text()
    if suffix=='.cif':
        export_repaired_structure(result,tmp_path/'loop-export.pdb')
        assert 'REMARK 950 ADDED' in (tmp_path/'loop-export.pdb').read_text()
    # Contacts with halogenated ligands and metals must not disappear from screening.
    for element in ('Cl','Zn'):
        context=after.clone();res=gemmi.Residue();res.name=element.upper();res.seqid=gemmi.SeqId(9999,' ');res.het_flag='H'
        atom=gemmi.Atom();atom.name=element.upper();atom.element=gemmi.Element(element)
        position=next(r for r in context[0][0] if r.seqid.num==3 and r.name==name).find_atom('CA','*').pos
        atom.pos=gemmi.Position(position.x+0.1,position.y,position.z);res.add_atom(atom);context[0][0].add_residue(res)
        findings,_=validate_loops(context,result.loops)
        assert any('overlap' in b for b in findings)
        if element=='Zn': assert any('specialist coordination' in b for b in findings)
    # Severe closure damage must be rejected independently of minimisation.
    modeled=next(r for r in after[0][0] if r.seqid.num==3 and r.name==name)
    modeled.find_atom('N','*').pos.x+=8
    blockers,_=validate_loops(after,result.loops)
    assert any('closure' in b or 'bond' in b for b in blockers)


def test_changed_sequence_and_mixed_operations_cannot_build(loop_source,tmp_path):
    source,full,_=loop_source
    st=gemmi.read_structure(str(source));st.write_minimal_pdb(str(source))
    plan=propose_repairs(source,sequence_references=[full]);choices=_loops(plan)
    full.write_text(full.read_text()+'REMARK changed\n')
    with pytest.raises(ValueError,match='evidence changed'):
        build_repairs(plan,[c.id for c in choices],tmp_path/'candidate')


def test_mixed_loop_and_atom_repairs_require_separate_stages(loop_source,tmp_path):
    source,_,_=loop_source
    st=gemmi.read_structure(str(source));st[0][0][0].remove_atom('O','*');st.write_pdb(str(source))
    plan=propose_repairs(source)
    loop=_loops(plan)[0];heavy=next(c for c in plan.changes if c.kind=='heavy_atoms')
    with pytest.raises(ValueError,match='Build loops separately'):
        build_repairs(plan,[loop.id,heavy.id],tmp_path/'candidate')


def test_long_loops_and_mismatched_original_coordinates_are_not_guessed(tmp_path,loop_source):
    path=structure_file(tmp_path,[('A',1,'GLY',GLY),('A',11,'SER',translated(GLY+[('CB',1.45,-1.4,0),('OG',1.45,-2.8,0)],15))],sequence=['GLY']+['ALA']*9+['SER'])
    plan=propose_repairs(path)
    assert not _loops(plan) and any('exceed' in n for n in plan.notes)
    source,full,_=loop_source
    st=gemmi.read_structure(str(source));st[0][0][0].find_atom('CA','*').pos.x+=0.1
    st.write_minimal_pdb(str(source))
    assert not _loops(propose_repairs(source,sequence_references=[full]))
