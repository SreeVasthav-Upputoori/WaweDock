"""Ramachandran backbone-torsion screening and plotting.

phi/psi angles are computed exactly from input coordinates with Gemmi. Residues
are then classified with a transparent, documented region screen into
``favored`` / ``allowed`` / ``outlier``.

This region screen is a simplification: it uses fixed phi/psi boxes for the
general, glycine and proline cases. It is NOT a MolProbity percentile-contour
analysis, and the boxes are deliberately generous so an ``outlier`` flags a
genuinely unusual backbone conformation for review, not a borderline one. No
coordinates are changed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import gemmi

# Region boxes as (phi_min, phi_max, psi_min, psi_max), degrees in (-180, 180].
# Favored boxes are subsets of the allowed boxes. psi wrap near +/-180 is
# covered with explicit boxes rather than modular arithmetic for clarity.
_GENERAL_ALLOWED = [
    (-175, -30, -100, 30),   # right-handed helix + bridge, broad
    (-180, -30, 60, 180),    # beta / polyproline-II / extended
    (-180, -45, -180, -120), # beta basin wrapped past +/-180
    (25, 100, -10, 95),      # left-handed helix region
]
_GENERAL_FAVORED = [
    (-160, -40, -75, 5),     # alpha-R core (right-handed helix), broad
    (-180, -45, 95, 180),    # beta / PPII core
    (35, 85, 15, 80),        # alpha-L core
]
# Proline: phi is restricted by the ring (~-63); two psi basins.
_PROLINE_ALLOWED = [
    (-100, -30, -80, 20),
    (-100, -30, 100, 180),
    (-100, -30, -180, -150),
]
_PROLINE_FAVORED = [
    (-85, -45, -60, 5),
    (-85, -45, 120, 180),
]


def _in(boxes, phi, psi):
    return any(a <= phi <= b and c <= psi <= d for a, b, c, d in boxes)


def classify(phi: float, psi: float, kind: str) -> str:
    """Return 'favored', 'allowed' or 'outlier' for one residue.

    kind is 'glycine', 'proline' or 'general'. Glycine is near-symmetric in
    phi/psi, so it is accepted if either (phi, psi) or its point reflection
    falls in the general regions.
    """
    if kind == 'proline':
        allowed, favored = _PROLINE_ALLOWED, _PROLINE_FAVORED
    elif kind == 'glycine':
        if _in(_GENERAL_FAVORED, phi, psi) or _in(_GENERAL_FAVORED, -phi, -psi):
            return 'favored'
        if _in(_GENERAL_ALLOWED, phi, psi) or _in(_GENERAL_ALLOWED, -phi, -psi):
            return 'allowed'
        return 'outlier'
    else:
        allowed, favored = _GENERAL_ALLOWED, _GENERAL_FAVORED
    if _in(favored, phi, psi):
        return 'favored'
    if _in(allowed, phi, psi):
        return 'allowed'
    return 'outlier'


@dataclass
class RamaPoint:
    chain: str
    seqid: int
    icode: str
    resname: str
    phi: float
    psi: float
    residue_class: str   # 'general' | 'glycine' | 'proline'
    classification: str  # 'favored' | 'allowed' | 'outlier'


def _primary(res, name):
    """Position of the named backbone atom, preferring the blank/first altloc."""
    best = None
    for atom in res:
        if atom.name != name or atom.is_hydrogen():
            continue
        if atom.altloc in ('\0', ' '):
            return atom.pos
        if best is None:
            best = atom.pos
    return best


def _residue_kind(name: str) -> str:
    if name == 'GLY':
        return 'glycine'
    if name == 'PRO':
        return 'proline'
    return 'general'


def compute_ramachandran(st, *, check_cancelled: Callable = lambda: None) -> list[RamaPoint]:
    """phi/psi for every amino-acid residue of model 1 that has both defined.

    phi needs the previous residue's C bonded to this N; psi needs the next
    residue's N bonded to this C (peptide bond <= 1.8 A, same chain). Residues
    at a break or chain end therefore have no entry.
    """
    points: list[RamaPoint] = []
    if not len(st):
        return points
    for chain in st[0]:
        residues = [res for res in chain
                    if (info := gemmi.find_tabulated_residue(res.name)) and info.is_amino_acid()]
        for i, res in enumerate(residues):
            check_cancelled()
            n, ca, c = _primary(res, 'N'), _primary(res, 'CA'), _primary(res, 'C')
            if n is None or ca is None or c is None:
                continue
            prev_c = _primary(residues[i - 1], 'C') if i > 0 else None
            next_n = _primary(residues[i + 1], 'N') if i + 1 < len(residues) else None
            if prev_c is None or next_n is None:
                continue
            # Require real peptide bonds so angles are not computed across gaps.
            if prev_c.dist(n) > 1.8 or c.dist(next_n) > 1.8:
                continue
            phi = _wrap(gemmi.calculate_dihedral(prev_c, n, ca, c))
            psi = _wrap(gemmi.calculate_dihedral(n, ca, c, next_n))
            kind = _residue_kind(res.name)
            points.append(RamaPoint(chain.name or '', res.seqid.num, (res.seqid.icode or '').strip(),
                                    res.name, round(phi, 1), round(psi, 1), kind, classify(phi, psi, kind)))
    return points


def _wrap(radians: float) -> float:
    import math
    degrees = math.degrees(radians)
    while degrees <= -180:
        degrees += 360
    while degrees > 180:
        degrees -= 360
    return degrees


def summarize(points: list[RamaPoint]) -> dict:
    """Counts and percentages by classification over all computed residues."""
    total = len(points)
    counts = {k: sum(p.classification == k for p in points) for k in ('favored', 'allowed', 'outlier')}
    percent = {k: round(100.0 * v / total, 1) if total else 0.0 for k, v in counts.items()}
    return {'total': total, 'counts': counts, 'percent': percent,
            'method': 'Region-based phi/psi screen (general/glycine/proline boxes); not MolProbity contours.'}


# --- Plot -------------------------------------------------------------------
# A self-contained SVG scatter on a white panel so it reads in any UI theme.
_W = 360          # plot square in px
_PAD = 38         # margin for axis labels
_FILL = {'favored': '#cfe8d0', 'allowed': '#fae6c2'}
_DOT = {'favored': '#2e7d32', 'allowed': '#e8820c', 'outlier': '#d32f2f'}


def _xy(phi, psi):
    x = _PAD + (phi + 180) / 360.0 * _W
    y = _PAD + (180 - psi) / 360.0 * _W   # psi increases upward
    return x, y


def _box_rect(box, fill):
    a, b, c, d = box
    x0, y0 = _xy(a, d)   # top-left: min phi, max psi
    x1, y1 = _xy(b, c)
    return f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="{fill}" />'


def render_svg(points: list[RamaPoint]) -> str:
    """Return a standalone SVG Ramachandran plot for the given residues."""
    full = _W + 2 * _PAD
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{full}" height="{full+26}" '
             f'viewBox="0 0 {full} {full+26}" font-family="sans-serif" font-size="11">',
             f'<rect x="0" y="0" width="{full}" height="{full+26}" fill="#ffffff"/>']
    # Allowed regions under favored regions, for every residue type (union view).
    for box in _GENERAL_ALLOWED + _PROLINE_ALLOWED:
        parts.append(_box_rect(box, _FILL['allowed']))
    for box in _GENERAL_FAVORED + _PROLINE_FAVORED:
        parts.append(_box_rect(box, _FILL['favored']))
    # Axes frame and zero lines.
    parts.append(f'<rect x="{_PAD}" y="{_PAD}" width="{_W}" height="{_W}" fill="none" stroke="#888" stroke-width="1"/>')
    for v in (-90, 0, 90):
        x, _ = _xy(v, 0)
        _, y = _xy(0, v)
        parts.append(f'<line x1="{x:.1f}" y1="{_PAD}" x2="{x:.1f}" y2="{_PAD+_W}" stroke="#ddd"/>')
        parts.append(f'<line x1="{_PAD}" y1="{y:.1f}" x2="{_PAD+_W}" y2="{y:.1f}" stroke="#ddd"/>')
    for v in (-180, -90, 0, 90, 180):
        x, _ = _xy(v, 0)
        _, y = _xy(0, v)
        parts.append(f'<text x="{x:.1f}" y="{_PAD+_W+14:.1f}" text-anchor="middle" fill="#555">{v}</text>')
        parts.append(f'<text x="{_PAD-6:.1f}" y="{y+4:.1f}" text-anchor="end" fill="#555">{v}</text>')
    parts.append(f'<text x="{_PAD+_W/2:.1f}" y="{full-2:.1f}" text-anchor="middle" fill="#333">phi (degrees)</text>')
    parts.append(f'<text x="12" y="{_PAD+_W/2:.1f}" text-anchor="middle" fill="#333" '
                 f'transform="rotate(-90 12 {_PAD+_W/2:.1f})">psi (degrees)</text>')
    # Residue points; outliers drawn last so they stay visible.
    for p in sorted(points, key=lambda q: q.classification != 'outlier'):
        x, y = _xy(p.phi, p.psi)
        r = 3.2 if p.classification == 'outlier' else 2.2
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{_DOT[p.classification]}" '
                     f'fill-opacity="0.85"><title>{p.chain}:{p.seqid}{p.icode} {p.resname} '
                     f'phi={p.phi:.0f} psi={p.psi:.0f} ({p.classification})</title></circle>')
    # Legend.
    lx = _PAD
    for label in ('favored', 'allowed', 'outlier'):
        parts.append(f'<circle cx="{lx+4:.1f}" cy="{full+18:.1f}" r="3.5" fill="{_DOT[label]}"/>')
        parts.append(f'<text x="{lx+12:.1f}" y="{full+22:.1f}" fill="#333">{label}</text>')
        lx += 70
    parts.append('</svg>')
    return ''.join(parts)
