"""Ramachandran backbone-torsion screening and plotting.

phi/psi angles are computed exactly from input coordinates with Gemmi. Residues
are then classified with a smooth density model of the backbone phi/psi
distribution: a Gaussian mixture whose basins are the right-handed helix, beta
sheet, polyproline-II and left-handed helix regions (with glycine treated
symmetrically and proline restricted). The same density draws the contour
background of the plot, so the shaded favored/allowed regions have the familiar
curved Ramachandran shapes.

This analytic density is a transparent approximation of the Ramachandran
distribution for screening and drawing; it is NOT the MolProbity empirical
percentile contour data. No coordinates are changed.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import gemmi

# Gaussian basins: (mean_phi, mean_psi, sd_phi, sd_psi, weight), degrees.
_GENERAL = [
    (-63.0, -43.0, 15.0, 17.0, 1.00),   # alpha-R, right-handed helix
    (-125.0, 128.0, 26.0, 28.0, 0.95),  # beta sheet
    (-68.0, 148.0, 14.0, 20.0, 0.50),   # polyproline-II
    (57.0, 47.0, 14.0, 18.0, 0.18),     # alpha-L, left-handed helix
]
# Proline phi is restricted by the ring (~-63); two psi basins.
_PROLINE = [
    (-63.0, -35.0, 12.0, 18.0, 0.90),
    (-63.0, 150.0, 12.0, 20.0, 0.90),
]


def _mirror(basins):
    return [(-a, -b, sa, sb, w) for (a, b, sa, sb, w) in basins]


# Glycine is near-symmetric in phi/psi, so mirror the general basins.
_GLYCINE = _GENERAL + _mirror(_GENERAL)
_BASINS = {'general': _GENERAL, 'glycine': _GLYCINE, 'proline': _PROLINE}

# Density thresholds as a fraction of each class's peak density, calibrated so a
# well-refined structure is ~90% favored with ~1-2% outliers (PROCHECK-like).
_FAVORED_FRAC = 0.0075
_ALLOWED_FRAC = 0.0002


def _circular(delta: float) -> float:
    delta = abs(delta) % 360.0
    return min(delta, 360.0 - delta)


def _density(phi: float, psi: float, basins) -> float:
    # phi is treated linearly (no basin lies near +/-180, and far-positive phi
    # is genuinely rare), while psi stays periodic so the beta basin wraps
    # correctly across +/-180.
    total = 0.0
    for mphi, mpsi, sphi, spsi, weight in basins:
        dp = (phi - mphi) / sphi
        dq = _circular(psi - mpsi) / spsi
        total += weight * math.exp(-0.5 * (dp * dp + dq * dq))
    return total


def _peak(basins) -> float:
    best = 0.0
    for phi in range(-180, 180, 2):
        for psi in range(-180, 180, 2):
            best = max(best, _density(float(phi), float(psi), basins))
    return best


_PEAK = {kind: _peak(basins) for kind, basins in _BASINS.items()}


def classify(phi: float, psi: float, kind: str = 'general') -> str:
    """Return 'favored', 'allowed' or 'outlier' from the density model.

    kind is 'glycine', 'proline' or 'general'.
    """
    key = kind if kind in _BASINS else 'general'
    relative = _density(phi, psi, _BASINS[key]) / _PEAK[key]
    if relative >= _FAVORED_FRAC:
        return 'favored'
    if relative >= _ALLOWED_FRAC:
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


def _wrap(radians: float) -> float:
    degrees = math.degrees(radians)
    while degrees <= -180:
        degrees += 360
    while degrees > 180:
        degrees -= 360
    return degrees


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
    """phi/psi for every amino-acid residue of model 1 that has both defined."""
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
            if prev_c.dist(n) > 1.8 or c.dist(next_n) > 1.8:
                continue
            phi = _wrap(gemmi.calculate_dihedral(prev_c, n, ca, c))
            psi = _wrap(gemmi.calculate_dihedral(n, ca, c, next_n))
            kind = _residue_kind(res.name)
            points.append(RamaPoint(chain.name or '', res.seqid.num, (res.seqid.icode or '').strip(),
                                    res.name, round(phi, 1), round(psi, 1), kind, classify(phi, psi, kind)))
    return points


def summarize(points: list[RamaPoint]) -> dict:
    """Counts and percentages by classification over all computed residues."""
    total = len(points)
    counts = {k: sum(p.classification == k for p in points) for k in ('favored', 'allowed', 'outlier')}
    percent = {k: round(100.0 * v / total, 1) if total else 0.0 for k, v in counts.items()}
    return {'total': total, 'counts': counts, 'percent': percent,
            'method': 'Gaussian-mixture density model of the Ramachandran distribution; not MolProbity contours.'}


# --- Plot -------------------------------------------------------------------
# Self-contained SVG on a white panel so it reads in any UI theme. The contour
# background is the general-case density; residue dots are coloured by their own
# per-type classification.
_W = 360          # plot square in px
_PAD = 38         # margin for axis labels
_STEP = 3         # contour grid resolution, degrees
_FILL = {'favored': '#cfe8d0', 'allowed': '#fae6c2'}
_DOT = {'favored': '#2e7d32', 'allowed': '#e8820c', 'outlier': '#d32f2f'}


def _xy(phi, psi):
    x = _PAD + (phi + 180) / 360.0 * _W
    y = _PAD + (180 - psi) / 360.0 * _W   # psi increases upward
    return x, y


def _contour_spans():
    """Horizontal (phi0, phi1, psi0, psi1, band) runs of the general contour."""
    basins, peak = _GENERAL, _PEAK['general']
    spans = []
    for row in range(-180, 180, _STEP):
        bands = []
        for col in range(-180, 180, _STEP):
            rel = _density(col + _STEP / 2.0, row + _STEP / 2.0, basins) / peak
            bands.append('favored' if rel >= _FAVORED_FRAC
                         else 'allowed' if rel >= _ALLOWED_FRAC else None)
        c = 0
        while c < len(bands):
            if bands[c] is None:
                c += 1
                continue
            band, start = bands[c], c
            while c < len(bands) and bands[c] == band:
                c += 1
            spans.append((-180 + start * _STEP, -180 + c * _STEP, row, row + _STEP, band))
    return spans


def render_svg(points: list[RamaPoint]) -> str:
    """Return a standalone SVG Ramachandran plot for the given residues."""
    full = _W + 2 * _PAD
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{full}" height="{full+26}" '
             f'viewBox="0 0 {full} {full+26}" font-family="sans-serif" font-size="11">',
             f'<rect x="0" y="0" width="{full}" height="{full+26}" fill="#ffffff"/>']
    # Smooth contour background (allowed under favored, by construction of bands).
    for phi0, phi1, psi0, psi1, band in _contour_spans():
        x0, y0 = _xy(phi0, psi1)
        x1, y1 = _xy(phi1, psi0)
        parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.2f}" height="{y1-y0:.2f}" '
                     f'fill="{_FILL[band]}" shape-rendering="crispEdges"/>')
    # Frame, zero lines, ticks.
    parts.append(f'<rect x="{_PAD}" y="{_PAD}" width="{_W}" height="{_W}" fill="none" stroke="#888" stroke-width="1"/>')
    for v in (-90, 0, 90):
        x, _ = _xy(v, 0)
        _, y = _xy(0, v)
        parts.append(f'<line x1="{x:.1f}" y1="{_PAD}" x2="{x:.1f}" y2="{_PAD+_W}" stroke="#bbb" stroke-dasharray="2,2"/>')
        parts.append(f'<line x1="{_PAD}" y1="{y:.1f}" x2="{_PAD+_W}" y2="{y:.1f}" stroke="#bbb" stroke-dasharray="2,2"/>')
    for v in (-180, -90, 0, 90, 180):
        x, _ = _xy(v, 0)
        _, y = _xy(0, v)
        parts.append(f'<text x="{x:.1f}" y="{_PAD+_W+14:.1f}" text-anchor="middle" fill="#555">{v}</text>')
        parts.append(f'<text x="{_PAD-6:.1f}" y="{y+4:.1f}" text-anchor="end" fill="#555">{v}</text>')
    parts.append(f'<text x="{_PAD+_W/2:.1f}" y="{full-2:.1f}" text-anchor="middle" fill="#333">phi (degrees)</text>')
    parts.append(f'<text x="12" y="{_PAD+_W/2:.1f}" text-anchor="middle" fill="#333" '
                 f'transform="rotate(-90 12 {_PAD+_W/2:.1f})">psi (degrees)</text>')
    # Residue points; outliers last so they stay visible.
    for p in sorted(points, key=lambda q: q.classification != 'outlier'):
        x, y = _xy(p.phi, p.psi)
        r = 3.2 if p.classification == 'outlier' else 2.2
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{_DOT[p.classification]}" '
                     f'stroke="#00000022" stroke-width="0.4"><title>{p.chain}:{p.seqid}{p.icode} {p.resname} '
                     f'phi={p.phi:.0f} psi={p.psi:.0f} ({p.classification})</title></circle>')
    # Legend.
    lx = _PAD
    for label in ('favored', 'allowed', 'outlier'):
        parts.append(f'<circle cx="{lx+4:.1f}" cy="{full+18:.1f}" r="3.5" fill="{_DOT[label]}"/>')
        parts.append(f'<text x="{lx+12:.1f}" y="{full+22:.1f}" fill="#333">{label}</text>')
        lx += 70
    parts.append('</svg>')
    return ''.join(parts)
