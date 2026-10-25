"""Per-residue protonation states surfaced from receptor preparation."""
import gemmi
import pytest

from mds.chem.receptor import PROTONATION_STATES, collect_protonation


def structure(residues):
    st = gemmi.Structure()
    model = gemmi.Model('1')
    chain = gemmi.Chain('A')
    for i, (number, name) in enumerate(residues, start=1):
        res = gemmi.Residue()
        res.name = name
        res.seqid = gemmi.SeqId(str(number))
        atom = gemmi.Atom()
        atom.name = 'CA'
        atom.element = gemmi.Element('C')
        atom.pos = gemmi.Position(float(i), 0.0, 0.0)
        res.add_atom(atom)
        chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    return st


def test_each_state_maps_to_its_standard_residue():
    found = collect_protonation(structure(
        [(1, 'HID'), (2, 'HIE'), (3, 'HIP'), (4, 'ASH'), (5, 'GLH'), (6, 'LYN'), (7, 'CYX')]))
    assert [p.state for p in found] == ['HID', 'HIE', 'HIP', 'ASH', 'GLH', 'LYN', 'CYX']
    assert [p.residue for p in found] == ['HIS', 'HIS', 'HIS', 'ASP', 'GLU', 'LYS', 'CYS']
    assert all(p.meaning for p in found)


def test_only_charge_changing_states_are_marked():
    found = {p.state: p.alters_charge for p in collect_protonation(structure(
        [(1, 'HID'), (2, 'HIE'), (3, 'HIP'), (4, 'ASH'), (5, 'GLH'), (6, 'LYN'), (7, 'CYX')]))}
    # These change the formal charge and shift binding-site electrostatics.
    assert found['HIP'] and found['ASH'] and found['GLH'] and found['LYN']
    # Tautomer and disulfide assignments do not.
    assert not found['HID'] and not found['HIE'] and not found['CYX']


@pytest.mark.parametrize('name,state', [('NHID', 'HID'), ('CHIP', 'HIP'), ('NCYX', 'CYX')])
def test_amber_terminal_prefixes_are_recognised(name, state):
    found = collect_protonation(structure([(1, name)]))
    assert len(found) == 1 and found[0].state == state


def test_standard_residues_are_not_reported():
    assert collect_protonation(structure([(1, 'GLY'), (2, 'HIS'), (3, 'CYS'), (4, 'ASP')])) == []


def test_location_is_recorded_for_review():
    found = collect_protonation(structure([(57, 'HIP')]))
    assert (found[0].chain, found[0].seqid) == ('A', 57)


def test_state_table_is_self_consistent():
    for state, (standard, meaning, alters) in PROTONATION_STATES.items():
        assert len(state) == 3 and len(standard) == 3
        assert isinstance(meaning, str) and meaning
        assert isinstance(alters, bool)
