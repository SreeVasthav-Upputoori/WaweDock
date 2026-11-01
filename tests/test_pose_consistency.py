"""Pose consistency: the six comparison axes, the consensus, and the verdicts.

The interaction axis is driven with hand-built contacts rather than a ProLIF run, because
``compare`` takes classified contacts as an argument precisely so it stays testable without
the optional Open Babel backend; ``test_interactions.py`` covers the real classification.
"""

import math

import pytest

Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend",
                           exc_type=ImportError)
from rdkit.Chem import AllChem

from mds.chem import pose_consistency as pc
from mds.chem.interactions import InteractionContact
from mds.chem.poses import mcs_rmsd


def _embedded(smiles, seed=42):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _translated(mol, dx=0.0, dy=0.0, dz=0.0):
    moved = Chem.Mol(mol)
    conformer = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        position = conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i, (position.x + dx, position.y + dy, position.z + dz))
    return moved


def _rotated(mol, degrees, axis=(0.0, 0.0, 1.0)):
    """Rotate about the molecule's own centroid, so only the orientation changes."""
    import numpy as np

    turned = Chem.Mol(mol)
    conformer = turned.GetConformer()
    xyz = np.array(conformer.GetPositions())
    centre = xyz.mean(axis=0)
    k = np.array(axis, dtype=float)
    k /= np.linalg.norm(k)
    theta = math.radians(degrees)
    cross = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    rotation = np.eye(3) + math.sin(theta) * cross + (1 - math.cos(theta)) * (cross @ cross)
    for i, point in enumerate((xyz - centre) @ rotation.T + centre):
        conformer.SetAtomPosition(i, tuple(float(v) for v in point))
    return turned


@pytest.fixture
def ligand():
    return _embedded("c1ccccc1CCNC(=O)C")        # anilide-like, clearly elongated


def _contacts(*pairs):
    return [InteractionContact(residue, interaction, interaction, (0,))
            for residue, interaction in pairs]


# ------------------------------------------------------------------------- geometry
def test_identical_pose_is_conserved_on_every_axis(ligand):
    result = pc.compare(ligand, Chem.Mol(ligand), "self")
    assert result.rmsd == 0.0
    assert result.mcs_rmsd == 0.0
    assert result.centroid_shift == 0.0
    assert result.axis_angle == 0.0
    assert result.volume_overlap == pytest.approx(100.0)
    assert result.feature_recovery == 1.0
    assert not result.flipped
    assert result.consensus == 1.0
    assert result.verdict == pc.CONSERVED


def test_displaced_pose_loses_geometry_and_shape_but_keeps_usr(ligand):
    moved = _translated(ligand, dx=9.0)
    result = pc.compare(ligand, moved, "displaced")
    assert result.rmsd > 8.0
    assert result.centroid_shift == pytest.approx(9.0, abs=1e-6)
    # A 9 A shift of an ~11 A molecule still leaves its tails brushing past each other, so a
    # little volume overlap survives; what matters is that almost none of it does.
    assert result.volume_overlap < 15.0
    # USR is built from internal distances only, so a pure translation cannot change it.
    assert result.usr == pytest.approx(1.0, abs=1e-3)
    assert result.verdict in (pc.DIFFERENT, pc.OUTLIER)
    assert result.consensus < pc.VERDICT_THRESHOLDS["partial"]


def test_mcs_rmsd_compares_a_different_ligand_on_the_shared_scaffold(ligand):
    """A substituted analogue has no whole-molecule RMSD but still has an MCS one."""
    analogue = _embedded("c1ccccc1CCNC(=O)CCl")
    # Put the analogue's benzene where the reference's is, so the scaffolds coincide.
    import numpy as np
    ref_ring = [a.GetIdx() for a in ligand.GetAtoms() if a.GetIsAromatic()]
    alt_ring = [a.GetIdx() for a in analogue.GetAtoms() if a.GetIsAromatic()]
    offset = (np.array(ligand.GetConformer().GetPositions())[ref_ring].mean(axis=0)
              - np.array(analogue.GetConformer().GetPositions())[alt_ring].mean(axis=0))
    analogue = _translated(analogue, *offset)

    result = pc.compare(ligand, analogue, "analogue")
    assert result.rmsd is None, "whole-molecule RMSD is undefined across different chemistry"
    assert result.same_chemistry is False
    assert result.mcs_rmsd is not None and result.mcs_atoms >= 6
    assert 0 < result.mcs_fraction <= 1.0
    assert result.best_rmsd == result.mcs_rmsd
    assert any("shared scaffold" in note for note in result.notes)


def test_mcs_rmsd_is_not_realigned(ligand):
    """Moving a pose must raise its MCS RMSD; an aligning implementation would not notice."""
    near, _atoms, _fraction = mcs_rmsd(Chem.Mol(ligand), ligand)
    far, _atoms, _fraction = mcs_rmsd(_translated(ligand, dx=5.0), ligand)
    assert near == pytest.approx(0.0, abs=1e-6)
    assert far == pytest.approx(5.0, abs=1e-3)


# ---------------------------------------------------------------------- orientation
def test_axis_angle_is_unsigned_and_a_flip_is_reported_separately(ligand):
    turned = _rotated(ligand, 180.0, axis=(0, 0, 1))
    result = pc.compare(ligand, turned, "turned")
    # A half turn maps the long axis onto its own reflection, so it lands back on nearly the
    # same line; the residual few degrees are the axis's component along the rotation axis.
    assert result.axis_angle < 15.0
    # ... but the molecule is threaded the other way round, which is what the flip test sees.
    assert result.flipped
    assert result.verdict == pc.FLIPPED


def test_flipped_verdict_outranks_a_high_consensus(ligand):
    """A reversed pose keeping most contacts must not be reported as conserved."""
    assert pc.verdict(0.95, flipped=True, ifp_tanimoto=0.8) == pc.FLIPPED
    assert pc.verdict(0.95, flipped=False, ifp_tanimoto=0.8) == pc.CONSERVED
    # Nothing in common at all is an outlier, flipped or not.
    assert pc.verdict(0.05, flipped=True, ifp_tanimoto=0.0) == pc.OUTLIER


def test_isotropic_ligand_reports_no_axis_angle():
    """Carbon tetrachloride has no long axis; a number there would be noise, not data."""
    ball = _embedded("ClC(Cl)(Cl)Cl")
    assert pc.long_axis(ball) is None
    result = pc.compare(ball, _rotated(ball, 55.0), "rotated ball")
    assert result.axis_angle is None
    assert any("isotropic" in note for note in result.notes)
    # The remaining axes still produce a verdict rather than failing.
    assert result.verdict in (pc.CONSERVED, pc.PARTIAL, pc.DIFFERENT, pc.OUTLIER)


# ---------------------------------------------------------------------- interactions
def test_interaction_tanimoto_splits_kept_lost_and_gained(ligand):
    reference = _contacts(("ASP25:A", "HBDonor"), ("PHE53:A", "PiStacking"),
                          ("ILE50:A", "Hydrophobic"))
    pose = _contacts(("ASP25:A", "HBDonor"), ("TRP60:A", "Hydrophobic"))
    result = pc.compare(ligand, Chem.Mol(ligand), "same place, other contacts",
                        reference_contacts=reference, pose_contacts=pose)
    assert result.ifp_tanimoto == pytest.approx(0.25)     # 1 shared of 4 in the union
    assert result.conserved == ("ASP25:A (HBDonor)",)
    assert set(result.missing) == {"PHE53:A (PiStacking)", "ILE50:A (Hydrophobic)"}
    assert result.gained == ("TRP60:A (Hydrophobic)",)


def test_van_der_waals_contacts_are_excluded_from_the_fingerprint():
    """They are reported across most of the interface, so they would flatter every pose."""
    reference = pc.interaction_set(_contacts(("ASP25:A", "HBDonor"), ("ILE50:A", "VdWContact")))
    pose = pc.interaction_set(_contacts(("ILE50:A", "VdWContact")))
    assert reference == {("ASP25:A", "HBDonor")}
    assert pose == frozenset()
    assert pc.ifp_similarity(reference, pose)["ifp_tanimoto"] == 0.0


def test_missing_contacts_disable_the_axis_rather_than_scoring_zero(ligand):
    result = pc.compare(ligand, Chem.Mol(ligand), "no backend")
    assert result.ifp_tanimoto is None
    assert any("Interaction fingerprint unavailable" in note for note in result.notes)
    # The consensus is taken over the other axes only, so a perfect pose still scores 1.0.
    assert result.consensus == 1.0


# --------------------------------------------------------------------- pharmacophore
def test_feature_recovery_needs_position_and_the_right_kind(ligand):
    matched, total, missed = pc.feature_recovery(Chem.Mol(ligand), ligand)
    assert total > 0 and matched == total and missed == ()

    moved = _translated(ligand, dx=10.0)
    matched, total, missed = pc.feature_recovery(moved, ligand)
    assert matched == 0 and total > 0 and missed

    # Benzene sits where the reference's ring does but offers no donor or acceptor.
    benzene = _embedded("c1ccccc1")
    matched, total, _missed = pc.feature_recovery(benzene, ligand)
    assert 0 <= matched < total


def test_feature_tolerance_is_honoured(ligand):
    nudged = _translated(ligand, dx=1.0)
    tight = pc.compare(ligand, nudged, "nudged",
                       options=pc.ConsistencyOptions(feature_tolerance=0.5))
    loose = pc.compare(ligand, nudged, "nudged",
                       options=pc.ConsistencyOptions(feature_tolerance=2.0))
    assert tight.matched_features < loose.matched_features


# -------------------------------------------------------------------------- scoring
def test_consensus_renormalises_over_the_axes_actually_computed():
    assert pc.consensus({"geometry": 1.0, "interactions": None, "features": None,
                         "shape": None, "orientation": None}) == 1.0
    assert pc.consensus({"geometry": 1.0, "interactions": 0.0}) == pytest.approx(0.5)
    assert pc.consensus({k: None for k in pc.METHOD_WEIGHTS}) == 0.0
    assert sum(pc.METHOD_WEIGHTS.values()) == pytest.approx(1.0)


def test_switching_axes_off_skips_their_work(ligand):
    bare = pc.compare(ligand, Chem.Mol(ligand), "bare",
                      options=pc.ConsistencyOptions(shape=False, features=False))
    assert bare.usr is None and bare.volume_overlap is None
    assert bare.feature_recovery is None
    assert bare.centroid_shift == 0.0           # orientation is always computed; it is cheap


def test_comparison_requires_coordinates(ligand):
    flat = Chem.MolFromSmiles("CCO")
    with pytest.raises(ValueError, match="3D coordinates"):
        pc.compare(ligand, flat, "no conformer")
    with pytest.raises(ValueError, match="3D coordinates"):
        pc.compare(flat, ligand, "no conformer")


def test_rows_flattens_results_for_export(ligand):
    reference = _contacts(("ASP25:A", "HBDonor"), ("PHE53:A", "PiStacking"))
    pose = _contacts(("ASP25:A", "HBDonor"))
    result = pc.compare(ligand, Chem.Mol(ligand), "pose A",
                        reference_contacts=reference, pose_contacts=pose)
    row = pc.rows([result])[0]
    assert row["Pose"] == "pose A" and row["Verdict"]
    assert row["Conserved contacts"] == "ASP25:A (HBDonor)"
    assert row["Lost contacts"] == "PHE53:A (PiStacking)"
    assert row["New contacts"] == ""
    assert all(value is not None for value in row.values())


def test_consistency_from_a_real_prolif_run():
    """End to end: ProLIF classifies both poses and the fingerprint axis separates them."""
    from mds.chem.interactions import compute_interactions

    def linear(smiles, start):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            conformer.SetAtomPosition(i, (start + i * 0.1, 0.0, 0.0))
        mol.AddConformer(conformer)
        return mol

    protein = linear("N", 0)
    for atom in protein.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("LYS")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)

    reference = linear("O", 2.0)
    reference_contacts = compute_interactions(protein, reference, diagram=False).contacts
    assert {c.interaction for c in reference_contacts} >= {"HBAcceptor"}

    near = pc.compare(reference, linear("O", 2.2), "near", reference_contacts=reference_contacts,
                      pose_contacts=compute_interactions(protein, linear("O", 2.2),
                                                         diagram=False).contacts)
    far = pc.compare(reference, linear("O", 9.0), "far", reference_contacts=reference_contacts,
                     pose_contacts=compute_interactions(protein, linear("O", 9.0),
                                                        diagram=False).contacts)
    assert near.ifp_tanimoto == 1.0 and near.verdict == pc.CONSERVED
    assert far.ifp_tanimoto == 0.0 and far.verdict == pc.OUTLIER
    assert far.missing == ("LYS1:A (HBAcceptor)",)
    # The van der Waals contact ProLIF also reports is excluded, so a Tanimoto of exactly 1.0
    # here means the real interaction matched -- not that a permissive contact padded both sets.
    assert near.conserved == ("LYS1:A (HBAcceptor)",)
