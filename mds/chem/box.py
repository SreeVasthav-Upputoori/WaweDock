"""Docking search box: construction from ligand / residues / whole protein, and validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

from mds.chem.report import Report

# Vina's search efficiency degrades for large volumes; these are advisory thresholds.
RECOMMENDED_MAX_EDGE = 30.0
RECOMMENDED_MAX_VOLUME = 27_000.0


@dataclass(frozen=True)
class DockingBox:
    center: tuple[float, float, float]
    size: tuple[float, float, float]

    @property
    def volume(self) -> float:
        return float(np.prod(self.size))

    @property
    def min_corner(self) -> np.ndarray:
        return np.array(self.center) - np.array(self.size) / 2

    @property
    def max_corner(self) -> np.ndarray:
        return np.array(self.center) + np.array(self.size) / 2

    def contains(self, xyz: np.ndarray) -> np.ndarray:
        return np.all((xyz >= self.min_corner) & (xyz <= self.max_corner), axis=1)

    def as_dict(self) -> dict[str, float]:
        return {**{f"center_{a}": c for a, c in zip("xyz", self.center)},
                **{f"size_{a}": s for a, s in zip("xyz", self.size)}}


def box_from_coords(coords: np.ndarray, padding: float = 8.0, cubic: bool = False,
                    min_edge: float = 15.0) -> DockingBox:
    """Box that encloses ``coords`` with ``padding`` Å on every side."""
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or len(coords) == 0:
        raise ValueError("coords must be an (N, 3) array with N >= 1")
    lo, hi = coords.min(axis=0), coords.max(axis=0)
    size = np.maximum(hi - lo + 2 * padding, min_edge)
    if cubic:
        size[:] = size.max()
    centre = (lo + hi) / 2
    return DockingBox(tuple(round(float(c), 3) for c in centre),
                      tuple(round(float(s), 3) for s in size))


def heavy_atom_coords(path: Path, *, chain: str | None = None,
                      residues: set[tuple[str, int]] | None = None) -> np.ndarray:
    """Heavy-atom coordinates from a PDB/mmCIF, optionally restricted to ``(chain, number)`` residues."""
    st = gemmi.read_structure(str(path))
    pts = []
    for ch in st[0]:
        if chain and ch.name != chain:
            continue
        for res in ch:
            if residues is not None and (ch.name, res.seqid.num) not in residues:
                continue
            if res.name in ("HOH", "DOD", "WAT"):
                continue
            pts += [[a.pos.x, a.pos.y, a.pos.z] for a in res if not a.is_hydrogen()]
    if not pts:
        raise ValueError("No atoms matched the selection")
    return np.array(pts)


def box_from_ligand_file(path: Path, padding: float = 8.0) -> DockingBox:
    """Box around a (co-crystal) ligand file: its heavy-atom extent plus padding."""
    return box_from_coords(heavy_atom_coords(Path(path)), padding)


def box_from_residues(receptor_pdb: Path, residues: set[tuple[str, int]],
                      padding: float = 6.0) -> DockingBox:
    return box_from_coords(heavy_atom_coords(Path(receptor_pdb), residues=residues), padding)


def blind_box(receptor_pdb: Path, padding: float = 4.0) -> DockingBox:
    """Box around the whole protein ('blind docking'); very large, needs high exhaustiveness."""
    return box_from_coords(heavy_atom_coords(Path(receptor_pdb)), padding, min_edge=0.0)


def validate_box(box: DockingBox, receptor_pdbqt: Path | None = None,
                 ligand_extent: float | None = None) -> Report:
    """Check the box for physical sense and report advisories (never silently altered)."""
    rep = Report()
    if any(s <= 0 for s in box.size):
        rep.error(f"Box dimensions must be positive, got {box.size}.")
        return rep
    if any(s > 126 for s in box.size):
        rep.error("Box edge exceeds 126 A, the practical maximum for Vina grid maps.")
    if max(box.size) > RECOMMENDED_MAX_EDGE or box.volume > RECOMMENDED_MAX_VOLUME:
        rep.warn(f"Large search space ({box.size[0]:.0f} x {box.size[1]:.0f} x {box.size[2]:.0f} A, "
                 f"{box.volume:,.0f} A^3). Sampling becomes less reliable; raise exhaustiveness "
                 "(e.g. 32+) or define a tighter box around the binding site.")
    if ligand_extent is not None and min(box.size) < ligand_extent:
        rep.warn(f"A box edge ({min(box.size):.1f} A) is shorter than the ligand's longest "
                 f"dimension ({ligand_extent:.1f} A); the ligand may be unable to fit in some orientations.")
    if receptor_pdbqt is not None:
        xyz = receptor_atom_coords(receptor_pdbqt)
        inside = int(box.contains(xyz).sum())
        if inside == 0:
            rep.error("The box contains no receptor atoms; it does not overlap the protein at all.")
        else:
            rep.info(f"{inside} receptor atoms lie inside the box.")
            if inside < 30:
                rep.warn("Very few receptor atoms inside the box; check the box centre.")
    return rep


def receptor_atom_coords(pdbqt: Path) -> np.ndarray:
    pts = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
           for l in Path(pdbqt).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    return np.array(pts)


def ligand_extent(coords: np.ndarray) -> float:
    """Largest pairwise distance between ligand atoms."""
    coords = np.asarray(coords)
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    return float(d.max())

