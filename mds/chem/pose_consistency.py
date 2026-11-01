"""Pose consistency: how closely docked poses reproduce a reference binding mode.

What this answers
-----------------
Docking gives a score per pose, and a score alone does not say whether a pose *binds the way
the reference does*. This module compares each docked pose with one trusted reference pose
(a co-crystal ligand, or a pose the user has decided to trust) along six independent axes:

===================  ==========================================================================
Geometry             in-place heavy-atom RMSD (same ligand) and MCS RMSD (different ligands)
Interactions         Tanimoto over the set of classified protein contacts, plus which
                     reference contacts are kept, lost or newly gained
Pharmacophore        how many of the reference's ligand-side features the pose reproduces
                     within a distance tolerance, feature kind by feature kind
Orientation          centre-of-mass shift and the angle between the two long axes
Shape (placement)    volume overlap of the two poses where they actually sit
Shape (conformer)    USR similarity, which is position-independent by construction
===================  ==========================================================================

Nothing here is realigned. Every geometric number is measured between the coordinates as
docked, because a pose that only matches the reference *after* being moved has not reproduced
the binding mode.

The consensus score
-------------------
:data:`METHOD_WEIGHTS` combines the axes into one number, and :data:`VERDICT_THRESHOLDS` turns
that into a verdict. Those weights and thresholds are **conventions, not a fitted or validated
model**: no activity data was used to choose them, and a different reasonable choice would move
borderline poses between verdicts. They exist to order a long pose list, and the individual
columns — not the consensus — are what should be quoted. See :data:`CAVEAT`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from rdkit import Chem

    from mds.chem.interactions import InteractionContact

import math
from dataclasses import dataclass, field

import numpy as np

from mds.chem.pharmacophore import FEATURE_LABELS, feature_points
from mds.chem.poses import centroid_shift, mcs_rmsd, symmetric_rmsd

CAVEAT = ("Pose consistency compares each pose with one reference pose. Every number is "
          "calculated from predicted coordinates, so it measures agreement with that "
          "reference, not correctness. The consensus score and the verdicts use the weights "
          "and cut-offs shown in the tooltips: they are conventions for ordering a pose list, "
          "not a validated model of binding, and the individual columns are what should be "
          "quoted.")

# Relative contribution of each axis to the consensus score. Only the axes that could actually
# be computed for a pose take part, and the weights are renormalised over those.
METHOD_WEIGHTS = {
    "geometry": 0.25,       # RMSD, or MCS RMSD for a different ligand
    "interactions": 0.25,   # interaction fingerprint Tanimoto
    "features": 0.20,       # pharmacophore feature recovery
    "shape": 0.15,          # volume overlap where the poses sit
    "orientation": 0.15,    # centroid shift and long-axis agreement
}

VERDICT_THRESHOLDS = {"conserved": 0.75, "partial": 0.40, "different": 0.15}

CONSERVED = "Conserved"
PARTIAL = "Partially conserved"
DIFFERENT = "Different"
FLIPPED = "Flipped"
OUTLIER = "Outlier"

# Van der Waals contacts are left out of the interaction fingerprint for the same reason they
# are left out of a pharmacophore: ProLIF reports them across most of the interface, so they
# raise every Tanimoto towards 1 without distinguishing binding modes.
IGNORED_INTERACTIONS = frozenset({"VdWContact"})


@dataclass(frozen=True)
class ConsistencyOptions:
    """Which axes to compute, and the cut-offs the verdict uses."""

    interactions: bool = True       # needs ProLIF + Open Babel; contacts are supplied by the caller
    shape: bool = True              # volume overlap is the slowest axis on large ligands
    features: bool = True
    rmsd_conserved: float = 2.0     # Å; the usual redocking success criterion
    rmsd_different: float = 4.0     # Å
    feature_tolerance: float = 2.0  # Å between a reference feature and the pose's own
    flip_angle: float = 120.0       # degrees between long axes
    mcs_timeout: int = 10           # seconds
    grid_spacing: float = 0.6       # Å, volume-overlap grid


@dataclass(frozen=True)
class PoseConsistency:
    """One pose measured against the reference. ``None`` means "not computed", never zero."""

    label: str
    same_chemistry: bool = False
    rmsd: float | None = None
    mcs_rmsd: float | None = None
    mcs_atoms: int = 0
    mcs_fraction: float = 0.0
    ifp_tanimoto: float | None = None
    conserved: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    gained: tuple[str, ...] = ()
    feature_recovery: float | None = None
    matched_features: int = 0
    reference_features: int = 0
    missed_features: tuple[str, ...] = ()
    centroid_shift: float | None = None
    axis_angle: float | None = None
    flipped: bool = False
    usr: float | None = None
    volume_overlap: float | None = None
    consensus: float = 0.0
    verdict: str = OUTLIER
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def best_rmsd(self) -> float | None:
        """The RMSD to quote: the whole-molecule one when the ligands match, else the MCS one."""
        return self.rmsd if self.rmsd is not None else self.mcs_rmsd


# --------------------------------------------------------------------------- orientation
def _heavy_coordinates(mol: Chem.Mol) -> np.ndarray:
    positions = np.array(mol.GetConformer().GetPositions())
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return positions[heavy] if heavy else positions


def long_axis(mol: Chem.Mol, *, anisotropy: float = 1.15) -> np.ndarray | None:
    """Unit vector along the ligand's longest extent, or None when that axis is not defined.

    The axis comes from the eigenvector of the heavy-atom covariance with the largest
    eigenvalue. For a ligand that is nearly as wide as it is long the two leading eigenvalues
    are almost equal, the eigenvector they span is arbitrary, and any angle measured from it is
    noise — so ``None`` is returned instead of a meaningless number. ``anisotropy`` is the
    ratio the leading eigenvalue must exceed.
    """
    coords = _heavy_coordinates(mol)
    if len(coords) < 3:
        return None
    centred = coords - coords.mean(axis=0)
    values, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    if values[1] <= 0 or values[0] / values[1] < anisotropy:
        return None
    return vectors[:, 0]


def axis_angle(probe: Chem.Mol, reference: Chem.Mol) -> float | None:
    """Angle in degrees between the two long axes, or None when either is undefined.

    A principal axis has no sign, so the angle is folded into 0-90°: two ligands lying along the
    same line agree whichever way round they are threaded. End-to-end reversal is a separate
    question and :func:`_reversal_angle` answers it from the atom order, not from the axis.
    """
    a, b = long_axis(probe), long_axis(reference)
    if a is None or b is None:
        return None
    return round(math.degrees(math.acos(min(1.0, abs(float(np.dot(a, b)))))), 1)


def _reversal_angle(probe: Chem.Mol, reference: Chem.Mol) -> float | None:
    """Angle between the two ligands' first-to-last heavy-atom vectors, 0-180°.

    This is what detects a flip: the long axis is unsigned, so it cannot tell a pose that is
    threaded head-to-tail from one that is tail-to-head. It needs the same atom ordering in
    both molecules, so it is only meaningful for the same ligand.
    """
    p, r = _heavy_coordinates(probe), _heavy_coordinates(reference)
    if len(p) < 2 or len(p) != len(r):
        return None
    pv, rv = p[-1] - p[0], r[-1] - r[0]
    if np.linalg.norm(pv) < 1e-6 or np.linalg.norm(rv) < 1e-6:
        return None
    cosine = float(np.dot(pv, rv) / (np.linalg.norm(pv) * np.linalg.norm(rv)))
    return round(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))), 1)


def orientation_score(shift: float | None, angle: float | None) -> float | None:
    """Agreement in placement and direction, 0-1. Both terms decay smoothly, neither is a cliff."""
    terms = []
    if shift is not None:
        terms.append(math.exp(-0.3 * shift))
    if angle is not None:
        terms.append(math.exp(-0.015 * angle))
    if not terms:
        return None
    return round(float(np.prod(terms)), 3)


# -------------------------------------------------------------------------------- shape
def _usr_descriptor(coords: np.ndarray) -> np.ndarray | None:
    """The 12 USR moments (Ballester & Richards 2007) of a heavy-atom cloud."""
    if len(coords) < 4:
        return None
    centroid = coords.mean(axis=0)
    from_centroid = np.linalg.norm(coords - centroid, axis=1)
    closest = coords[int(np.argmin(from_centroid))]
    farthest = coords[int(np.argmax(from_centroid))]
    from_farthest = np.linalg.norm(coords - farthest, axis=1)
    beyond = coords[int(np.argmax(from_farthest))]

    moments = []
    for reference in (centroid, closest, farthest, beyond):
        d = np.linalg.norm(coords - reference, axis=1)
        mean, sd = float(d.mean()), float(d.std())
        skew = float(np.mean(((d - mean) / sd) ** 3)) if sd > 1e-9 else 0.0
        moments += [mean, sd, skew]
    return np.array(moments)


def usr_similarity(probe: Chem.Mol, reference: Chem.Mol) -> float | None:
    """USR shape similarity, 0-1.

    USR is built from internal distances only, so it is invariant to where the pose sits: it
    compares the *conformers*, not their placement. It is reported separately for that reason
    and is deliberately kept out of the consensus score, where :func:`volume_overlap` carries
    the shape term instead.
    """
    a = _usr_descriptor(_heavy_coordinates(probe))
    b = _usr_descriptor(_heavy_coordinates(reference))
    if a is None or b is None:
        return None
    return round(1.0 / (1.0 + float(np.abs(a - b).sum()) / 12.0), 3)


def volume_overlap(probe: Chem.Mol, reference: Chem.Mol, *, spacing: float = 0.6,
                   max_points: int = 4_000_000) -> float | None:
    """Jaccard overlap of the two poses' van der Waals volumes, as a percentage.

    Computed on a grid over the union's bounding box with each atom's own vdW radius. Because
    the coordinates are left where they are, this measures whether the two poses occupy the
    same space — not merely whether they are the same shape.
    """
    from rdkit import Chem

    table = Chem.GetPeriodicTable()

    def spheres(mol: Chem.Mol):
        positions = np.array(mol.GetConformer().GetPositions())
        keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        if not keep:
            return None, None
        radii = np.array([table.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in keep])
        return positions[keep], radii

    p_xyz, p_r = spheres(probe)
    r_xyz, r_r = spheres(reference)
    if p_xyz is None or r_xyz is None:
        return None

    pad = float(max(p_r.max(), r_r.max()))
    low = np.minimum(p_xyz.min(axis=0), r_xyz.min(axis=0)) - pad
    high = np.maximum(p_xyz.max(axis=0), r_xyz.max(axis=0)) + pad
    while np.prod(np.ceil((high - low) / spacing) + 1) > max_points:
        spacing *= 2.0
    axes = [np.arange(low[i], high[i] + spacing, spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)

    def occupied(xyz: np.ndarray, radii: np.ndarray) -> np.ndarray:
        inside = np.zeros(len(grid), dtype=bool)
        for centre, radius in zip(xyz, radii):
            inside |= np.sum((grid - centre) ** 2, axis=1) <= radius * radius
        return inside

    a, b = occupied(p_xyz, p_r), occupied(r_xyz, r_r)
    union = int(np.count_nonzero(a | b))
    if not union:
        return None
    return round(int(np.count_nonzero(a & b)) / union * 100, 1)


# ------------------------------------------------------------------------- interactions
def interaction_set(contacts: list[InteractionContact]) -> frozenset[tuple[str, str]]:
    """The ``(residue, interaction)`` pairs a pose makes, which is the fingerprint compared."""
    return frozenset((c.residue, c.interaction) for c in contacts
                     if c.interaction not in IGNORED_INTERACTIONS)


def _describe(pairs) -> tuple[str, ...]:
    return tuple(sorted(f"{residue} ({interaction})" for residue, interaction in pairs))


def ifp_similarity(reference: frozenset, pose: frozenset) -> dict:
    """Tanimoto over the two interaction sets, with what was kept, lost and gained."""
    union = reference | pose
    return {
        "ifp_tanimoto": round(len(reference & pose) / len(union), 3) if union else None,
        "conserved": _describe(reference & pose),
        "missing": _describe(reference - pose),
        "gained": _describe(pose - reference),
    }


# ------------------------------------------------------------------------- pharmacophore
def feature_recovery(probe: Chem.Mol, reference: Chem.Mol, *,
                     tolerance: float = 2.0) -> tuple[int, int, tuple[str, ...]]:
    """How many of the reference's own pharmacophore features the pose places nearby.

    Features come from the same RDKit definitions the interaction pharmacophore uses
    (:func:`mds.chem.pharmacophore.feature_points`), and are matched only within their own kind,
    so a donor can never be satisfied by a ring. Each reference feature is paired with the
    closest unused pose feature of that kind, shortest distance first — a greedy assignment
    rather than an optimal one, which keeps this dependency-free and is equivalent whenever the
    matches are unambiguous, the only case in which the number is worth reading.

    Returns ``(matched, reference_total, labels_of_the_misses)``.
    """
    ref_points = feature_points(reference)
    pose_points = feature_points(probe)
    matched, total, missed = 0, 0, []
    for kind, wanted in ref_points.items():
        total += len(wanted)
        available = list(pose_points.get(kind, ()))
        pairs = sorted(((math.dist(w, c), i, j) for i, w in enumerate(wanted)
                        for j, c in enumerate(available) if math.dist(w, c) <= tolerance))
        used_reference: set[int] = set()
        used_pose: set[int] = set()
        for _distance, i, j in pairs:
            if i not in used_reference and j not in used_pose:
                used_reference.add(i)
                used_pose.add(j)
        matched += len(used_reference)
        if len(used_reference) < len(wanted):
            missed.append(FEATURE_LABELS.get(kind, kind))
    return matched, total, tuple(sorted(missed))


# ------------------------------------------------------------------------------ scoring
def _geometry_score(rmsd: float | None, options: ConsistencyOptions) -> float | None:
    if rmsd is None:
        return None
    if rmsd <= options.rmsd_conserved:
        return 1.0
    if rmsd >= options.rmsd_different:
        return 0.0
    span = options.rmsd_different - options.rmsd_conserved
    return round(1.0 - (rmsd - options.rmsd_conserved) / span, 3) if span > 0 else 0.0


def consensus(axes: dict[str, float | None]) -> float:
    """Weighted mean of whichever axes were computed, renormalised over their weights."""
    weight = sum(METHOD_WEIGHTS[k] for k, v in axes.items() if v is not None)
    if weight <= 0:
        return 0.0
    total = sum(METHOD_WEIGHTS[k] * v for k, v in axes.items() if v is not None)
    return round(total / weight, 3)


def verdict(score: float, *, flipped: bool, ifp_tanimoto: float | None) -> str:
    """Name the binding mode. A flip is called out separately because the score cannot show it.

    A pose that is reversed end to end can still keep most of the reference's contacts and
    score respectably; reporting that as conserved would hide the one thing a modeller most
    needs to see, so the flip takes precedence whenever any contact survives.
    """
    if flipped and (ifp_tanimoto is None or ifp_tanimoto > 0):
        return FLIPPED
    if score >= VERDICT_THRESHOLDS["conserved"]:
        return CONSERVED
    if score >= VERDICT_THRESHOLDS["partial"]:
        return PARTIAL
    if score >= VERDICT_THRESHOLDS["different"]:
        return DIFFERENT
    return OUTLIER


# ------------------------------------------------------------------------------ compare
def compare(reference: Chem.Mol, pose: Chem.Mol, label: str, *,
            options: ConsistencyOptions | None = None,
            reference_contacts: list[InteractionContact] | None = None,
            pose_contacts: list[InteractionContact] | None = None) -> PoseConsistency:
    """Measure one pose against the reference.

    Classified interactions are passed in rather than computed here, because they need the
    receptor and the optional Open Babel backend; when they are absent that axis is reported as
    unavailable and the consensus is taken over the remaining ones.
    """
    options = options or ConsistencyOptions()
    if not reference.GetNumConformers() or not pose.GetNumConformers():
        raise ValueError("Pose consistency needs 3D coordinates for both the pose and the reference")
    notes: list[str] = []

    rmsd = None
    try:
        rmsd = round(symmetric_rmsd(pose, reference), 3)
    except ValueError:
        notes.append("Different chemistry from the reference: compared on the shared scaffold "
                     "(MCS RMSD), not whole-molecule RMSD.")
    try:
        mcs, mcs_atoms, mcs_fraction = mcs_rmsd(pose, reference, timeout=options.mcs_timeout)
    except Exception as exc:  # noqa: BLE001 — one unmatchable ligand must not stop the batch
        mcs, mcs_atoms, mcs_fraction = None, 0, 0.0
        notes.append(f"Common-substructure search failed: {exc}")
    if rmsd is None and mcs is None:
        notes.append("No common substructure of three or more atoms: the geometric axes could "
                     "not be computed.")

    ifp: dict = {"ifp_tanimoto": None, "conserved": (), "missing": (), "gained": ()}
    if options.interactions:
        if reference_contacts is None or pose_contacts is None:
            notes.append("Interaction fingerprint unavailable (ProLIF/Open Babel, or no receptor).")
        else:
            ifp = ifp_similarity(interaction_set(reference_contacts), interaction_set(pose_contacts))

    matched = total = 0
    missed: tuple[str, ...] = ()
    recovery = None
    if options.features:
        matched, total, missed = feature_recovery(pose, reference,
                                                  tolerance=options.feature_tolerance)
        recovery = round(matched / total, 3) if total else None
        if total == 0:
            notes.append("The reference has no recognised pharmacophore features.")

    shift = round(centroid_shift(pose, reference), 3)
    angle = axis_angle(pose, reference)
    if angle is None:
        notes.append("Long axis undefined for one of the two ligands (too isotropic), so no "
                     "axis angle or flip test.")
    reversal = _reversal_angle(pose, reference)
    flipped = reversal is not None and reversal > options.flip_angle

    usr = volume = None
    if options.shape:
        usr = usr_similarity(pose, reference)
        volume = volume_overlap(pose, reference, spacing=options.grid_spacing)

    axes = {
        "geometry": _geometry_score(rmsd if rmsd is not None else mcs, options),
        "interactions": ifp["ifp_tanimoto"],
        "features": recovery,
        "shape": None if volume is None else round(volume / 100.0, 3),
        "orientation": orientation_score(shift, angle),
    }
    score = consensus(axes)
    return PoseConsistency(
        label=label, same_chemistry=rmsd is not None, rmsd=rmsd,
        mcs_rmsd=mcs, mcs_atoms=mcs_atoms, mcs_fraction=mcs_fraction,
        ifp_tanimoto=ifp["ifp_tanimoto"], conserved=ifp["conserved"],
        missing=ifp["missing"], gained=ifp["gained"],
        feature_recovery=recovery, matched_features=matched, reference_features=total,
        missed_features=missed, centroid_shift=shift, axis_angle=angle, flipped=flipped,
        usr=usr, volume_overlap=volume, consensus=score,
        verdict=verdict(score, flipped=flipped, ifp_tanimoto=ifp["ifp_tanimoto"]),
        notes=tuple(notes),
    )


CSV_COLUMNS = [
    ("label", "Pose"), ("verdict", "Verdict"), ("consensus", "Consensus"),
    ("rmsd", "RMSD (A)"), ("mcs_rmsd", "MCS RMSD (A)"), ("mcs_atoms", "MCS atoms"),
    ("mcs_fraction", "MCS fraction of reference"), ("ifp_tanimoto", "Interaction Tanimoto"),
    ("feature_recovery", "Feature recovery"), ("matched_features", "Features matched"),
    ("reference_features", "Reference features"), ("centroid_shift", "Centroid shift (A)"),
    ("axis_angle", "Long-axis angle (deg)"), ("flipped", "Flipped"),
    ("usr", "USR (conformer only)"), ("volume_overlap", "Volume overlap (%)"),
]


def rows(results: list[PoseConsistency]) -> list[dict]:
    """Flatten results for CSV export, with the kept/lost/gained contacts spelled out."""
    out = []
    for r in results:
        row = {header: getattr(r, name) for name, header in CSV_COLUMNS}
        row["Conserved contacts"] = "; ".join(r.conserved)
        row["Lost contacts"] = "; ".join(r.missing)
        row["New contacts"] = "; ".join(r.gained)
        row["Notes"] = " ".join(r.notes)
        out.append({k: ("" if v is None else v) for k, v in row.items()})
    return out
