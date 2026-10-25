"""Tests for cysteine disulfide-bond detection."""
from pathlib import Path

import gemmi

from mds.chem.disulfides import find_disulfides
from mds.chem.health import check_protein

DATA = Path(__file__).parent / 'data' / '3PTB.pdb'


def _structure(sg_positions):
    """Build a model with one CYS per given SG (x, y, z) position."""
    st = gemmi.Structure()
    model = gemmi.Model('1')
    chain = gemmi.Chain('A')
    for i, (x, y, z) in enumerate(sg_positions, start=1):
        res = gemmi.Residue()
        res.name = 'CYS'
        res.seqid = gemmi.SeqId(str(i))
        atom = gemmi.Atom()
        atom.name = 'SG'
        atom.element = gemmi.Element('S')
        atom.pos = gemmi.Position(x, y, z)
        res.add_atom(atom)
        chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    return st


def test_bonded_pair_within_window():
    bonds, free = find_disulfides(_structure([(0, 0, 0), (2.05, 0, 0)]))
    assert len(bonds) == 1 and not free
    assert bonds[0].distance == 2.05


def test_distant_cysteines_are_free_not_bonded():
    bonds, free = find_disulfides(_structure([(0, 0, 0), (6.0, 0, 0)]))
    assert not bonds
    assert len(free) == 2


def test_distance_window_edges():
    # Just inside vs just outside the 1.8-2.5 A window.
    assert len(find_disulfides(_structure([(0, 0, 0), (2.4, 0, 0)]))[0]) == 1
    assert len(find_disulfides(_structure([(0, 0, 0), (2.6, 0, 0)]))[0]) == 0


def test_real_structure_trypsin_has_six_bridges():
    bonds, free = find_disulfides(gemmi.read_structure(str(DATA)))
    assert len(bonds) == 6
    assert all(1.8 <= b.distance <= 2.5 for b in bonds)
    assert free == []


def test_report_carries_disulfides_and_issues():
    report = check_protein(DATA)
    assert len(report.disulfides) == 6
    ss_issues = [i for i in report.issues if i.code == 'disulfide_bond']
    assert len(ss_issues) == 6
    assert all(issue.atoms for issue in ss_issues)  # highlightable in the viewer
    assert any('Disulfide' in line for line in report.coverage)
