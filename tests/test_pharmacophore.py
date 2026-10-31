"""Interaction pharmacophores: derivation from classified contacts, and pose matching.

These run without Open Babel: the ProLIF path is exercised with a synthetic protein molecule
built directly as RDKit, as in ``test_interactions.py``.
"""

import pytest

Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend",
                           exc_type=ImportError)
from rdkit.Chem import AllChem

from mds.chem import pharmacophore as ph
from mds.chem.interactions import InteractionContact, compute_interactions


def _embedded(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    return mol


def _linear(smiles, position):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (position + i * 0.1, 0, 0))
    mol.AddConformer(conf)
    return mol


def _translated(mol, dx):
    moved = Chem.Mol(mol)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + dx, pos.y, pos.z))
    return moved


@pytest.fixture
def phenol():
    mol = _embedded("c1ccccc1O")
    oxygen = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O")
    ring = tuple(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    return mol, oxygen, ring


def test_derivation_types_features_from_the_ligand_side(phenol):
    mol, oxygen, ring = phenol
    pharmacophore = ph.derive([
        InteractionContact("ASP25:A", "HBDonor", "x", (oxygen,)),
        InteractionContact("THR26:A", "HBAcceptor", "x", (oxygen,)),
        InteractionContact("PHE53:A", "PiStacking", "x", ring),
    ], mol, source="phenol")
    assert pharmacophore.counts() == {"Aromatic": 1, "HBA": 1, "HBD": 1}
    assert pharmacophore.source == "phenol"
    assert ph.CAVEAT in pharmacophore.to_json()
    aromatic = next(f for f in pharmacophore.features if f.kind == "Aromatic")
    assert aromatic.label == "Aromatic ring"
    # the ring feature sits at the ring centroid, not on any one atom
    conf = mol.GetConformer()
    centroid = [sum(conf.GetAtomPosition(i).x for i in ring) / len(ring),
                sum(conf.GetAtomPosition(i).y for i in ring) / len(ring)]
    assert aromatic.x == pytest.approx(centroid[0]) and aromatic.y == pytest.approx(centroid[1])
    # polar features are tighter than bulk hydrophobic/aromatic contact
    assert next(f for f in pharmacophore.features if f.kind == "HBD").tolerance < aromatic.tolerance


def test_van_der_waals_is_excluded_and_repeats_merge(phenol):
    mol, oxygen, ring = phenol
    pharmacophore = ph.derive([
        InteractionContact("ILE50:A", "VdWContact", "x", (oxygen,)),
        InteractionContact("PHE53:A", "PiStacking", "x", ring),
        InteractionContact("TYR59:A", "PiStacking", "x", ring),
        InteractionContact("TRP60:A", "HBDonor", "x", ()),        # no indices reported
    ], mol)
    assert len(pharmacophore) == 1
    assert pharmacophore.features[0].kind == "Aromatic"
    # one aromatic requirement naming both partners, not two features
    assert pharmacophore.features[0].residues == ("PHE53:A", "TYR59:A")


def test_matching_requires_both_position_and_chemistry(phenol):
    mol, oxygen, ring = phenol
    pharmacophore = ph.derive([
        InteractionContact("ASP25:A", "HBDonor", "x", (oxygen,)),
        InteractionContact("PHE53:A", "PiStacking", "x", ring),
    ], mol)

    exact = ph.match(pharmacophore, mol)
    assert exact.n_matched == 2 and exact.recovery == 1.0 and not exact.missed()
    assert all(m.distance == pytest.approx(0.0, abs=1e-6) for m in exact.matches)

    moved = ph.match(pharmacophore, _translated(mol, 5.0))
    assert moved.n_matched == 0 and moved.recovery == 0.0
    assert {f.kind for f in moved.missed()} == {"HBD", "Aromatic"}
    assert all(m.distance is not None for m in moved.matches)   # nearest atom still reported

    # Benzene holds the ring in place but has no donor: right position, wrong chemistry.
    benzene = Chem.Mol(mol)
    editable = Chem.RWMol(benzene)
    for index in sorted([oxygen] + [n.GetIdx() for n in benzene.GetAtomWithIdx(oxygen).GetNeighbors()
                                    if n.GetSymbol() == "H"], reverse=True):
        editable.RemoveAtom(index)
    stripped = editable.GetMol()
    Chem.SanitizeMol(stripped)
    partial = ph.match(pharmacophore, stripped)
    assert [f.kind for f in partial.missed()] == ["HBD"]
    assert partial.n_matched == 1 and partial.recovery == pytest.approx(0.5)


def test_tolerance_boundary_is_inclusive_and_overridable(phenol):
    mol, oxygen, _ring = phenol
    contacts = [InteractionContact("ASP25:A", "HBDonor", "x", (oxygen,))]
    tight = ph.derive(contacts, mol, tolerances={"HBD": 0.5})
    assert tight.features[0].tolerance == 0.5
    assert ph.match(tight, _translated(mol, 0.5)).n_matched == 1     # exactly on the boundary
    assert ph.match(tight, _translated(mol, 0.6)).n_matched == 0


def test_metal_direction_follows_prolif_naming():
    """MetalDonor means the ligand IS the metal; MetalAcceptor means it chelates one."""
    zinc = _linear("[Zn+2]", 2)
    ion = ph.derive([InteractionContact("LYS1:A", "MetalDonor", "x", (0,))], zinc)
    assert ion.features[0].kind == "MetalIon"
    assert ph.match(ion, zinc).n_matched == 1
    # an oxygen at the same point is not a metal ion, however well placed
    assert ph.match(ion, _linear("O", 2)).n_matched == 0

    water = _linear("O", 2)
    binder = ph.derive([InteractionContact("ZN1:A", "MetalAcceptor", "x", (0,))], water)
    assert binder.features[0].kind == "MetalBinder"
    assert ph.match(binder, water).n_matched == 1
    assert ph.match(binder, _linear("[Zn+2]", 2)).n_matched == 0


def test_derivation_rejects_mismatched_pose_and_missing_conformer(phenol):
    mol, oxygen, _ring = phenol
    contacts = [InteractionContact("ASP25:A", "HBDonor", "x", (mol.GetNumAtoms() + 3,))]
    with pytest.raises(ValueError, match="do not match"):
        ph.derive(contacts, mol)
    with pytest.raises(ValueError, match="no conformer"):
        ph.derive([InteractionContact("ASP25:A", "HBDonor", "x", (oxygen,))],
                  Chem.MolFromSmiles("CCO"))
    with pytest.raises(ValueError, match="no conformer"):
        ph.match(ph.derive([], mol), Chem.MolFromSmiles("CCO"))


def test_empty_pharmacophore_reports_no_recovery(phenol):
    mol, _oxygen, _ring = phenol
    empty = ph.derive([], mol)
    assert len(empty) == 0 and empty.counts() == {}
    result = ph.match(empty, mol)
    assert result.recovery == 0.0 and result.n_matched == 0


def test_pharmacophore_from_a_real_prolif_run():
    """End-to-end: ProLIF classifies, the contacts carry atom indices, features follow."""
    protein = _linear("N", 0)
    for atom in protein.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("LYS")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)
    ligand = _linear("[Zn+2]", 2)

    result = compute_interactions(protein, ligand, diagram=False)
    assert result.contacts, "interaction contacts should be reported alongside the residue table"
    assert {c.interaction for c in result.contacts} >= {"MetalDonor"}
    assert all(c.ligand_atoms and max(c.ligand_atoms) < ligand.GetNumAtoms()
               for c in result.contacts)

    pharmacophore = ph.derive(result.contacts, ligand, source="zinc pose")
    assert pharmacophore.counts() == {"MetalIon": 1}   # the van der Waals contact is not a feature
    assert ph.match(pharmacophore, ligand).recovery == 1.0


def test_ring_matching_compares_centroids_not_ring_atoms(phenol):
    """A superposed ring must match with its full tolerance to spare, not 0.2 A."""
    mol, _oxygen, ring = phenol
    pharmacophore = ph.derive([InteractionContact("PHE53:A", "PiStacking", "x", ring)], mol)
    feature = pharmacophore.features[0]

    conformer = mol.GetConformer()
    atom_distances = [
        ((conformer.GetAtomPosition(i).x - feature.x) ** 2
         + (conformer.GetAtomPosition(i).y - feature.y) ** 2
         + (conformer.GetAtomPosition(i).z - feature.z) ** 2) ** 0.5 for i in ring]
    assert min(atom_distances) > 1.0, "ring atoms sit well off the centroid, as expected"

    exact = ph.match(pharmacophore, mol)
    assert exact.matches[0].distance == pytest.approx(0.0, abs=1e-6)

    # Tighten below the centroid-to-atom distance: centroid matching still succeeds.
    tight = ph.derive([InteractionContact("PHE53:A", "PiStacking", "x", ring)], mol,
                      tolerances={"Aromatic": 0.3})
    assert ph.match(tight, mol).n_matched == 1
