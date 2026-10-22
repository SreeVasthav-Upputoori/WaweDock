"""Scientific boundaries for protein diagnostics, independent of RDKit."""
from pathlib import Path
import json

import gemmi
import pytest

from mds.chem.health import check_protein


def structure_file(tmp_path, residues, *, sequence=None, suffix='.pdb'):
    structure = gemmi.Structure()
    model = gemmi.Model('1')
    chains = {}
    for chain_id, number, name, atoms in residues:
        if chain_id not in chains: chains[chain_id] = gemmi.Chain(chain_id)
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(str(number))
        residue.het_flag = 'A'
        for atom_name, x,y,z,*alt in atoms:
            atom = gemmi.Atom()
            atom.name = atom_name
            atom.element = gemmi.Element(atom_name[0])
            atom.pos = gemmi.Position(x,y,z)
            if alt: atom.altloc = alt[0]
            residue.add_atom(atom)
        chains[chain_id].add_residue(residue)
    for chain in chains.values(): model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    if sequence:
        structure.entities[0].full_sequence = sequence
    path = tmp_path / ('protein'+suffix)
    if suffix == '.cif': structure.make_mmcif_document().write_file(str(path))
    else: structure.write_pdb(str(path))
    return path


GLY = [('N',0,0,0),('CA',1.45,0,0),('C',2.0,1.42,0),('O',1.5,2.54,0)]


def translated(atoms, dx):
    return [(name,x+dx,y,z,*rest) for name,x,y,z,*rest in atoms]


def test_legitimate_termini_and_no_seqres_numbering_gap(tmp_path):
    path = structure_file(tmp_path,[('A',1,'GLY',GLY),('A',90,'GLY',translated(GLY,20)),
                                    ('B',1,'GLY',translated(GLY,40))])
    original = path.read_bytes()
    report = check_protein(path)
    assert path.read_bytes() == original
    assert not any(i.code in ('missing_atoms','unmodelled_sequence') for i in report.issues)
    breaks = [i for i in report.issues if i.code == 'backbone_connectivity']
    assert len(breaks) == 1
    assert 'A:1-90' == breaks[0].location
    assert not any(i.code=='steric_overlap' for i in report.issues)


def test_oxt_at_chain_terminus_is_not_a_clash(tmp_path):
    atoms = GLY + [('OXT',3.2,1.45,0)]
    report = check_protein(structure_file(tmp_path,[('A',1,'GLY',atoms)]))
    assert not any(i.code in ('steric_overlap','missing_atoms') for i in report.issues)


def test_missing_sidechain_and_conformer_specific_completeness(tmp_path):
    atoms = GLY + [('CB',1.4,-1.4,0,'A')]
    atoms = [(name,x,y,z,'B') if name=='CA' else (name,x,y,z,*rest) for name,x,y,z,*rest in atoms]
    report = check_protein(structure_file(tmp_path,[('A','10A','ALA',atoms)]))
    missing = [i for i in report.issues if i.code=='missing_atoms']
    assert len(missing)==2
    assert any("Conformer B" in i.message and "CB" in i.message for i in missing)
    assert any("Conformer A" in i.message and "CA" in i.message for i in missing)
    assert '10A' in missing[0].location.upper()
    assert any(i.code=='alternate_conformations' for i in report.issues)
    assert all('atom_index' in a for i in report.issues for a in i.atoms)


def test_entity_sequence_supports_unmodelled_terminal_region(tmp_path):
    path = structure_file(tmp_path,[('A',7,'GLY',GLY)], sequence=['GLY','ALA'],suffix='.cif')
    report = check_protein(path)
    missing = [i for i in report.issues if i.code=='unmodelled_sequence']
    assert len(missing)==1
    assert missing[0].severity=='Note'
    assert 'sequence terminus' in missing[0].message


def test_clash_and_backbone_geometry_screens(tmp_path):
    path = structure_file(tmp_path,[('A',1,'GLY',GLY),('B',1,'GLY',translated(GLY,0.1))])
    report = check_protein(path)
    assert any(i.code=='steric_overlap' for i in report.issues)
    bad = GLY.copy(); bad[1]=('CA',8,0,0)
    report = check_protein(structure_file(tmp_path,[('A',1,'GLY',bad)]))
    assert any(i.code=='bond_geometry' for i in report.issues)


def test_export_discloses_checks_and_unsupported_validation(tmp_path):
    report = check_protein(Path('tests/data/3PTB.pdb'))
    report.export(tmp_path/'report.json')
    data = json.loads((tmp_path/'report.json').read_text())
    assert len(data['source_sha256'])==64
    assert any('Ramachandran' in item for item in data['limitations'])
    assert not any(i.code=='steric_overlap' for i in report.issues)
    assert any(i.code=='structural_component' for i in report.issues)


def test_cancel_stops_computation(tmp_path):
    path=structure_file(tmp_path,[('A',1,'GLY',GLY)])
    def cancel(): raise RuntimeError('cancelled')
    with pytest.raises(RuntimeError,match='cancelled'):
        check_protein(path,check_cancelled=cancel)

def test_ambiguous_seqres_alignment_does_not_invent_missing_regions(tmp_path):
    path=structure_file(tmp_path,[('A',5,'GLY',GLY)],sequence=['GLY','GLY'])
    assert 'SEQRES' in path.read_text()
    report=check_protein(path)
    assert not any(i.code=='unmodelled_sequence' for i in report.issues)
    assert any('no unambiguous' in text for text in report.limitations)


def test_ter_separated_segments_are_not_backbone_breaks(tmp_path):
    path=structure_file(tmp_path,[('A',1,'GLY',GLY),('A',2,'GLY',translated(GLY,30))])
    text=path.read_text()
    lines=text.splitlines()
    atom_lines=[i for i,line in enumerate(lines) if line.startswith('ATOM')]
    lines.insert(atom_lines[4],'TER')
    path.write_text('\n'.join(lines)+'\n')
    report=check_protein(path)
    assert not any(i.code=='backbone_connectivity' for i in report.issues)
