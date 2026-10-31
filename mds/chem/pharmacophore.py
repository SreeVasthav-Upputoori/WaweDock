"""Interaction pharmacophores derived from a docked pose, and pose matching against them.

What this does
--------------
A docked pose plus its ProLIF classification says *which ligand atoms* make *which kind* of
interaction with *which residue*. Turning that into a pharmacophore means keeping the chemistry
and the position but dropping the specific ligand: an H-bond donor was needed near this point,
an aromatic ring near that one. Other poses — of the same ligand or a different one — can then be
asked how much of that pattern they reproduce.

The intended use is a reference-led one: build the pharmacophore from a crystal ligand's pose
(see ``mds.render.pose_review`` for the reference-locking workflow), then score docked poses by
how much of the known binding mode they recover. Feature recovery is reported as a count and a
fraction, never as an affinity or a probability of activity.

Why not a statistical pharmacophore model
-----------------------------------------
Commercial suites fit pharmacophore hypotheses over ligand sets and validate them against actives
and decoys. That needs labelled activity data this program does not have, so nothing here is
fitted, weighted or trained: a feature exists because an interaction was classified at those
coordinates, and matching is a geometric test against an explicit tolerance. The derived
pharmacophore inherits every limitation of the pose it came from — see :data:`CAVEAT`.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit import Chem

    from mds.chem.interactions import InteractionContact

CAVEAT = ("A pharmacophore derived here describes one docked or crystal pose: the features are "
          "placed where that pose's interactions were classified. It is a hypothesis about what "
          "a ligand needs at those positions, not a validated activity model, and feature "
          "recovery is a geometric score, not a predicted affinity.")

# Ligand-side feature type for each ProLIF interaction. The ligand's own role decides the type:
# "HBDonor" means the LIGAND donates, so the ligand needs a donor there.  Van der Waals contacts
# are deliberately absent -- they are too permissive to constrain anything, and ProLIF reports
# them for most of the interface.
FEATURE_OF_INTERACTION = {
    "HBDonor": "HBD", "HBAcceptor": "HBA",
    "XBDonor": "Halogen", "XBAcceptor": "HBA",
    "Hydrophobic": "Hydrophobic",
    "PiStacking": "Aromatic", "EdgeToFace": "Aromatic", "FaceToFace": "Aromatic",
    "PiCation": "Aromatic",      # ligand ring, protein cation
    "CationPi": "Cation",        # ligand cation, protein ring
    "Cationic": "Cation", "Anionic": "Anion",
    # ProLIF's names are from the ligand's side: MetalDonor means the LIGAND is the metal,
    # MetalAcceptor means the ligand chelates a metal in the protein.  They need opposite
    # chemistry of a candidate ligand, so they are separate feature types.
    "MetalDonor": "MetalIon", "MetalAcceptor": "MetalBinder",
}

FEATURE_LABELS = {
    "HBD": "H-bond donor", "HBA": "H-bond acceptor", "Hydrophobic": "Hydrophobic",
    "Aromatic": "Aromatic ring", "Cation": "Positive ionisable",
    "Anion": "Negative ionisable", "Halogen": "Halogen",
    "MetalIon": "Metal ion", "MetalBinder": "Metal-coordinating",
}

# Matching radius per feature type, in angstrom. Directional polar features get a tighter
# tolerance than bulk hydrophobic/aromatic contact, whose position is a centroid of several atoms
# and therefore less precisely defined. These are conventional starting values, not fitted ones.
DEFAULT_TOLERANCES = {
    "HBD": 1.2, "HBA": 1.2, "Cation": 1.3, "Anion": 1.3, "Halogen": 1.2,
    "MetalIon": 1.2, "MetalBinder": 1.2,
    "Hydrophobic": 1.8, "Aromatic": 1.6,
}

# RDKit BaseFeatures families that can satisfy each feature type, plus halogens by element:
# RDKit's feature definitions have no halogen-bond family.
_SATISFIED_BY = {
    "HBD": ("Donor",), "HBA": ("Acceptor",),
    "Cation": ("PosIonizable",), "Anion": ("NegIonizable",),
    "Hydrophobic": ("Hydrophobe", "LumpedHydrophobe"),
    "Aromatic": ("Aromatic",),
    "MetalBinder": ("Acceptor",),   # plus any N/O/S by element, see _feature_points
    # Matched by element instead of an RDKit feature family: its definitions have no
    # halogen-bond family, and no family marks an atom as being a metal.
    "Halogen": (), "MetalIon": (),
}
_HALOGENS = {9, 17, 35, 53}
# Everything heavier than helium that is not one of these is treated as a metal. Listing the
# non-metals and metalloids is shorter and more stable than listing every metal.
_NON_METALS = {1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86}
# Metal coordination uses a lone pair, which is broader than RDKit's H-bond Acceptor
# definition: it excludes thioethers, and gives water no features at all.
_LONE_PAIR = {7, 8, 16}


@dataclass(frozen=True)
class PharmacophoreFeature:
    kind: str                    # key of FEATURE_LABELS
    x: float
    y: float
    z: float
    tolerance: float             # matching radius, angstrom
    residues: tuple[str, ...]    # partnering protein residues, e.g. ("ASP25:A",)

    @property
    def label(self) -> str:
        return FEATURE_LABELS.get(self.kind, self.kind)

    @property
    def center(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass
class Pharmacophore:
    features: list[PharmacophoreFeature]
    source: str = ""             # what the features were derived from, for provenance
    caveat: str = CAVEAT

    def __len__(self) -> int:
        return len(self.features)

    def counts(self) -> dict[str, int]:
        """How many features of each kind, for a summary line."""
        out: dict[str, int] = {}
        for feature in self.features:
            out[feature.kind] = out.get(feature.kind, 0) + 1
        return out

    def to_json(self) -> str:
        return json.dumps({"source": self.source, "caveat": self.caveat,
                           "features": [asdict(f) for f in self.features]}, indent=2)


@dataclass(frozen=True)
class FeatureMatch:
    feature: PharmacophoreFeature
    matched: bool
    distance: float | None = None     # to the nearest satisfying atom, when one was found


@dataclass
class MatchResult:
    matches: list[FeatureMatch] = field(default_factory=list)

    @property
    def n_matched(self) -> int:
        return sum(1 for m in self.matches if m.matched)

    @property
    def recovery(self) -> float:
        """Fraction of the pharmacophore's features this pose satisfies (0.0 with no features)."""
        return self.n_matched / len(self.matches) if self.matches else 0.0

    def missed(self) -> list[PharmacophoreFeature]:
        return [m.feature for m in self.matches if not m.matched]


def derive(contacts: list[InteractionContact], pose: Chem.Mol, *,
           source: str = "", tolerances: dict[str, float] | None = None) -> Pharmacophore:
    """Build a pharmacophore from one pose's classified interactions.

    Each contact becomes a feature at the centroid of the ligand atoms responsible for it.
    Contacts of the same feature type on the same atoms are merged into one feature that lists
    every partnering residue, so a ring stacking against two residues is one aromatic
    requirement rather than two. Contacts whose interaction type has no ligand-side feature
    (van der Waals) or whose atom indices ProLIF did not report are skipped.
    """
    radii = {**DEFAULT_TOLERANCES, **(tolerances or {})}
    conformer = pose.GetConformer() if pose.GetNumConformers() else None
    if conformer is None:
        raise ValueError("The pose has no conformer, so pharmacophore features have no position")

    merged: dict[tuple[str, tuple[int, ...]], list[str]] = {}
    for contact in contacts:
        kind = FEATURE_OF_INTERACTION.get(contact.interaction)
        if kind is None or not contact.ligand_atoms:
            continue
        atoms = tuple(sorted(contact.ligand_atoms))
        if max(atoms) >= pose.GetNumAtoms():
            raise ValueError(f"Interaction references atom {max(atoms)}, outside the pose's "
                             f"{pose.GetNumAtoms()} atoms — pose and interactions do not match")
        residues = merged.setdefault((kind, atoms), [])
        if contact.residue not in residues:
            residues.append(contact.residue)

    features = []
    for (kind, atoms), residues in merged.items():
        positions = [conformer.GetAtomPosition(i) for i in atoms]
        n = len(positions)
        features.append(PharmacophoreFeature(
            kind=kind,
            x=sum(p.x for p in positions) / n,
            y=sum(p.y for p in positions) / n,
            z=sum(p.z for p in positions) / n,
            tolerance=radii.get(kind, 1.5),
            residues=tuple(residues),
        ))
    # Stable, readable order: by feature type, then by position.
    features.sort(key=lambda f: (f.kind, round(f.x, 3), round(f.y, 3), round(f.z, 3)))
    return Pharmacophore(features, source=source)


def match(pharmacophore: Pharmacophore, pose: Chem.Mol) -> MatchResult:
    """Score how much of ``pharmacophore`` this pose reproduces.

    A feature is satisfied when the pose offers compatible chemistry within the feature's
    tolerance. Compatible chemistry comes from RDKit's own feature definitions, so a donor
    requirement is not satisfied by an atom that merely happens to be nearby.

    Positions are compared feature to feature, not feature to atom: a ring's requirement sits at
    a ring centroid, and the candidate ring's centroid is what has to reach it. Measuring to the
    candidate's individual ring atoms instead would leave a perfectly superposed six-ring 1.4 A
    adrift of its own centroid and barely inside tolerance.
    """
    if not pose.GetNumConformers():
        raise ValueError("The pose has no conformer, so it cannot be matched")

    points = _feature_points(pose)
    matches = []
    for feature in pharmacophore.features:
        distances = [math.dist(feature.center, point) for point in points.get(feature.kind, ())]
        best = min(distances) if distances else None
        matches.append(FeatureMatch(feature, best is not None and best <= feature.tolerance, best))
    return MatchResult(matches)


Point = tuple[float, float, float]


def _feature_points(mol: Chem.Mol) -> dict[str, tuple[Point, ...]]:
    """Positions in ``mol`` that can satisfy each pharmacophore feature type.

    RDKit places a multi-atom feature (an aromatic ring, a lumped hydrophobe) at its centroid,
    which is what a derived feature of that kind is compared against. Element-matched types have
    no RDKit family, so they contribute their atom positions.
    """
    import os.path

    from rdkit import Chem, RDConfig  # noqa: PLC0415 — chemistry import stays lazy
    from rdkit.Chem import ChemicalFeatures

    factory = ChemicalFeatures.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
    by_family: dict[str, list[Point]] = {}
    for feature in factory.GetFeaturesForMol(Chem.Mol(mol)):
        position = feature.GetPos()
        by_family.setdefault(feature.GetFamily(), []).append((position.x, position.y, position.z))

    conformer = mol.GetConformer()

    def atoms(predicate) -> list[Point]:
        out = []
        for atom in mol.GetAtoms():
            if predicate(atom.GetAtomicNum()):
                position = conformer.GetAtomPosition(atom.GetIdx())
                out.append((position.x, position.y, position.z))
        return out

    out: dict[str, tuple[Point, ...]] = {}
    for kind, families in _SATISFIED_BY.items():
        points: list[Point] = [p for family in families for p in by_family.get(family, ())]
        if kind == "Halogen":
            points += atoms(lambda n: n in _HALOGENS)
        elif kind == "MetalIon":
            points += atoms(lambda n: n > 2 and n not in _NON_METALS)
        elif kind == "MetalBinder":
            points += atoms(lambda n: n in _LONE_PAIR)
        out[kind] = tuple(points)
    return out
