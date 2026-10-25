"""Chain termini classified as true sequence ends or truncation artefacts."""
from pathlib import Path

from mds.chem.health import build_assessment, check_protein
from test_protein_health import GLY, structure_file, translated

DATA = Path(__file__).parent / 'data' / '3PTB.pdb'


def termini(report):
    return [i for i in report.issues if i.code == 'artificial_terminus']


def chain(names):
    """Backbone-only residues with the given names, spaced along x."""
    return [('A', number, name, translated(GLY, 3.3 * i))
            for i, (number, name) in enumerate(names)]


def test_truncated_chain_flags_both_ends(tmp_path):
    # Entity sequence ALA-GLY-PHE-GLY-ALA; only residues 2-4 are modelled, so
    # both ends are truncations. Distinct names keep the alignment unambiguous.
    path = structure_file(tmp_path, chain([(2, 'GLY'), (3, 'PHE'), (4, 'GLY')]),
                          sequence=['ALA', 'GLY', 'PHE', 'GLY', 'ALA'], suffix='.cif')
    found = termini(check_protein(path))
    assert len(found) == 2
    assert sorted(i.message[0] for i in found) == ['C', 'N']
    assert all(i.severity == 'Warning' and i.atoms for i in found)
    assert all('1 entity-sequence residue' in i.message for i in found)


def test_fully_modelled_chain_has_no_artificial_terminus(tmp_path):
    path = structure_file(tmp_path, chain([(1, 'GLY'), (2, 'PHE'), (3, 'GLY')]),
                          sequence=['GLY', 'PHE', 'GLY'], suffix='.cif')
    assert termini(check_protein(path)) == []


def test_only_one_end_truncated(tmp_path):
    # Residues 1-3 of a four-residue entity: the C-terminus alone is artificial.
    path = structure_file(tmp_path, chain([(1, 'GLY'), (2, 'PHE'), (3, 'GLY')]),
                          sequence=['GLY', 'PHE', 'GLY', 'ALA'], suffix='.cif')
    found = termini(check_protein(path))
    assert len(found) == 1 and found[0].message.startswith('C-terminus')


def test_complete_real_structure_reports_no_truncation():
    # 3PTB models its full 223-residue entity, so neither end is an artefact.
    report = check_protein(DATA)
    assert termini(report) == []
    assert report.assessment['metrics']['artificial_termini'] == 0


def test_termini_are_reported_but_do_not_penalise_the_score(tmp_path):
    path = structure_file(tmp_path, chain([(2, 'GLY'), (3, 'PHE'), (4, 'GLY')]),
                          sequence=['ALA', 'GLY', 'PHE', 'GLY', 'ALA'], suffix='.cif')
    report = check_protein(path)
    assert report.assessment['metrics']['artificial_termini'] == 2
    # Truncation is a preparation concern, not a geometry defect: no penalty.
    residues = report.assessment['metrics']['residues']
    before = report.assessment['score']
    report.issues = [i for i in report.issues if i.code != 'artificial_terminus']
    assert build_assessment(report, residues)['score'] == before


def test_coverage_and_limitation_document_the_method():
    report = check_protein(DATA)
    assert any('termini' in line.lower() for line in report.coverage)
    assert any('deposited entity sequence' in line for line in report.limitations)
