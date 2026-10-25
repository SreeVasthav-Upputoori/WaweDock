"""Disulfide-bond detection from cysteine SG-SG distance (read-only).

A disulfide is called when two cysteine SG atoms lie within a covalent S-S
distance window. This is a geometric assignment from input coordinates only;
no bonds or labels are added to the structure, and oxidation-state residue
names are not required or trusted.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

# Covalent S-S is ~2.05 A; accept a generous window for experimental models.
_MIN, _MAX = 1.8, 2.5
# Names a disulfide cysteine may carry in experimental or prepared inputs.
_CYS_NAMES = ('CYS', 'CYX', 'CYM')


@dataclass
class Disulfide:
    chain_a: str
    seqid_a: int
    icode_a: str
    chain_b: str
    seqid_b: int
    icode_b: str
    distance: float


def _sg(res):
    """SG position of a cysteine, preferring the blank/first altloc."""
    best = None
    for atom in res:
        if atom.name == 'SG' and not atom.is_hydrogen():
            if atom.altloc in ('\0', ' '):
                return atom.pos
            if best is None:
                best = atom.pos
    return best


def find_disulfides(st, *, check_cancelled: Callable = lambda: None):
    """Return (disulfides, free_cysteines) for model 1.

    free_cysteines lists cysteine SG atoms not within bonding distance of any
    other SG, as dicts {chain, seqid, icode}.
    """
    sgs = []  # (chain, seqid, icode, pos)
    if not len(st):
        return [], []
    for chain in st[0]:
        for res in chain:
            if res.name in _CYS_NAMES:
                pos = _sg(res)
                if pos is not None:
                    sgs.append((chain.name or '', res.seqid.num, (res.seqid.icode or '').strip(), pos))
    bonds: list[Disulfide] = []
    bonded: set[int] = set()
    for i in range(len(sgs)):
        check_cancelled()
        for j in range(i + 1, len(sgs)):
            distance = sgs[i][3].dist(sgs[j][3])
            if _MIN <= distance <= _MAX:
                a, b = sgs[i], sgs[j]
                bonds.append(Disulfide(a[0], a[1], a[2], b[0], b[1], b[2], round(distance, 2)))
                bonded.add(i)
                bonded.add(j)
    free = [{'chain': sgs[k][0], 'seqid': sgs[k][1], 'icode': sgs[k][2]}
            for k in range(len(sgs)) if k not in bonded]
    return bonds, free
