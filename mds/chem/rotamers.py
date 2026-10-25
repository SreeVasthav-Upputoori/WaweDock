"""Side-chain rotamer screening from staggered chi-torsion geometry (read-only).

For every residue with rotatable sp3 side-chain bonds, chi angles are computed
exactly with Gemmi and compared to the three ideal staggered positions
(gauche+, trans, gauche- ~ +60/180/-60). A residue is an ``outlier`` when any
screened chi sits far from all three (an eclipsed/strained conformation).

This is a transparent geometric screen, NOT a rotamer-library or MolProbity
percentile analysis. Only sp3-sp3 chi angles are screened; planar or aromatic
terminal chi angles (carboxylate, amide, aromatic ring, guanidinium) are
skipped because they are not staggered and would be misread. No coordinates
are changed.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import gemmi

# Screened sp3-sp3 chi angles only (chi1 first). Planar/aromatic terminal chi
# angles are intentionally omitted: ASP/ASN/HIS/PHE/TYR/TRP chi2, GLN/GLU chi3,
# ARG chi4+. Proline is excluded (ring-constrained, not freely rotameric).
CHI_ATOMS = {
    'ARG': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'NE')],
    'ASN': [('N', 'CA', 'CB', 'CG')],
    'ASP': [('N', 'CA', 'CB', 'CG')],
    'CYS': [('N', 'CA', 'CB', 'SG')],
    'GLN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
    'GLU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
    'HIS': [('N', 'CA', 'CB', 'CG')],
    'ILE': [('N', 'CA', 'CB', 'CG1'), ('CA', 'CB', 'CG1', 'CD1')],
    'LEU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'LYS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'CE'), ('CG', 'CD', 'CE', 'NZ')],
    'MET': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'SD'), ('CB', 'CG', 'SD', 'CE')],
    'PHE': [('N', 'CA', 'CB', 'CG')],
    'SER': [('N', 'CA', 'CB', 'OG')],
    'THR': [('N', 'CA', 'CB', 'OG1')],
    'TRP': [('N', 'CA', 'CB', 'CG')],
    'TYR': [('N', 'CA', 'CB', 'CG')],
    'VAL': [('N', 'CA', 'CB', 'CG1')],
}
_STAGGERED = (-60.0, 60.0, 180.0)
_FAVORED_TOL = 30.0   # within this of a staggered minimum -> clearly rotameric
_ALLOWED_TOL = 40.0   # beyond this of every minimum -> eclipsed/strained


@dataclass
class RotamerPoint:
    chain: str
    seqid: int
    icode: str
    resname: str
    chis: list        # measured chi angles, degrees
    classification: str  # 'favored' | 'allowed' | 'outlier'


def _wrap(radians: float) -> float:
    degrees = math.degrees(radians)
    while degrees <= -180:
        degrees += 360
    while degrees > 180:
        degrees -= 360
    return degrees


def _atom(res, name):
    best = None
    for atom in res:
        if atom.name != name or atom.is_hydrogen():
            continue
        if atom.altloc in ('\0', ' '):
            return atom.pos
        if best is None:
            best = atom.pos
    return best


def _deviation(angle: float) -> float:
    """Smallest angular distance (degrees) from any staggered minimum."""
    def diff(a, b):
        d = abs(a - b) % 360
        return min(d, 360 - d)
    return min(diff(angle, s) for s in _STAGGERED)


def classify(deviations) -> str:
    worst = max(deviations)
    if worst > _ALLOWED_TOL:
        return 'outlier'
    if worst > _FAVORED_TOL:
        return 'allowed'
    return 'favored'


def compute_rotamers(st, *, check_cancelled: Callable = lambda: None) -> list[RotamerPoint]:
    """Rotamer classification for model-1 residues with screened chi angles."""
    points: list[RotamerPoint] = []
    if not len(st):
        return points
    for chain in st[0]:
        for res in chain:
            quartets = CHI_ATOMS.get(res.name)
            if not quartets:
                continue
            check_cancelled()
            chis, deviations = [], []
            for quartet in quartets:
                positions = [_atom(res, name) for name in quartet]
                if any(p is None for p in positions):
                    break  # incomplete side chain; stop in chi order
                angle = _wrap(gemmi.calculate_dihedral(*positions))
                chis.append(round(angle, 1))
                deviations.append(_deviation(angle))
            if not deviations:
                continue
            points.append(RotamerPoint(chain.name or '', res.seqid.num, (res.seqid.icode or '').strip(),
                                       res.name, chis, classify(deviations)))
    return points


def summarize(points: list[RotamerPoint]) -> dict:
    total = len(points)
    counts = {k: sum(p.classification == k for p in points) for k in ('favored', 'allowed', 'outlier')}
    percent = {k: round(100.0 * v / total, 1) if total else 0.0 for k, v in counts.items()}
    return {'total': total, 'counts': counts, 'percent': percent,
            'method': 'Staggered sp3 chi-torsion screen (+/-60/180); not a rotamer library or MolProbity contours.'}
