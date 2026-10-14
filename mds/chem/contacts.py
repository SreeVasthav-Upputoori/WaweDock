"""Distance-based protein-ligand contact analysis (deliberately conservative).

Everything reported here is a *geometric candidate* derived from distances only.  Angles,
protonation, water mediation and energetics are not evaluated, so none of these entries
is a verified interaction.  Use a dedicated tool (ProLIF/PLIP) and experimental evidence
before interpreting them chemically.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

CONTACT_CUTOFF = 4.5
HBOND_MAX = 3.5          # heavy-atom N/O ... N/O distance
SALT_BRIDGE_MAX = 4.0    # charged group centres
HYDROPHOBIC_MAX = 4.0    # C...C between non-polar carbons

_POS = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ")}
_NEG = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CAVEAT = ("Distance-based candidates only; angles, protonation and water-mediation were not "
          "evaluated. These are calculated geometric features, not verified interactions.")


@dataclass
class Contact:
    residue: str
    min_distance: float
    candidates: list[str] = field(default_factory=list)


@dataclass
class _ProteinAtoms:
    xyz: np.ndarray
    residue: list[str]
    resname: list[str]
    name: list[str]
    element: list[str]
    polar_carbon: np.ndarray


def load_protein_atoms(pdb: Path) -> _ProteinAtoms:
    st = gemmi.read_structure(str(pdb))
    xyz, residue, resname, name, element = [], [], [], [], []
    for ch in st[0]:
        for res in ch:
            for a in res:
                if a.is_hydrogen():
                    continue
                xyz.append([a.pos.x, a.pos.y, a.pos.z])
                residue.append(f"{res.name}{res.seqid.num}:{ch.name}")
                resname.append(res.name)
                name.append(a.name)
                element.append(a.element.name)
    # carbons bonded to N/O/S within 1.7 A are treated as polar (not hydrophobic)
    arr = np.array(xyz)
    polar_c = np.zeros(len(arr), dtype=bool)
    hetero = np.array([e in ("N", "O", "S") for e in element])
    for i, e in enumerate(element):
        if e == "C":
            d = np.linalg.norm(arr[hetero] - arr[i], axis=1) if hetero.any() else np.array([9.9])
            polar_c[i] = bool((d < 1.7).any())
    return _ProteinAtoms(arr, residue, resname, name, element, polar_c)


def pose_contacts(pose: Chem.Mol, protein: _ProteinAtoms, cutoff: float = CONTACT_CUTOFF) -> list[Contact]:
    lig = Chem.RemoveHs(pose, sanitize=False)
    lxyz = np.array(lig.GetConformer().GetPositions())
    latoms = list(lig.GetAtoms())
    l_el = [a.GetSymbol() for a in latoms]
    l_charge = np.array([a.GetFormalCharge() for a in latoms])
    l_polar_c = np.array([a.GetSymbol() == "C" and any(n.GetSymbol() in ("N", "O", "S")
                                                      for n in a.GetNeighbors()) for a in latoms])

    d = np.linalg.norm(protein.xyz[:, None, :] - lxyz[None, :, :], axis=-1)   # (P, L)
    contacts: dict[str, Contact] = {}
    for pi in np.where(d.min(axis=1) <= cutoff)[0]:
        res = protein.residue[pi]
        c = contacts.setdefault(res, Contact(res, float("inf")))
        c.min_distance = min(c.min_distance, float(d[pi].min()))
        for li in np.where(d[pi] <= cutoff)[0]:
            dist = float(d[pi, li])
            pe, le = protein.element[pi], l_el[li]
            if pe in ("N", "O") and le in ("N", "O") and dist <= HBOND_MAX:
                _add(c, f"H-bond candidate: {protein.name[pi]}...{le}{li + 1} {dist:.2f} A")
            if pe == "C" and le == "C" and dist <= HYDROPHOBIC_MAX and \
                    not protein.polar_carbon[pi] and not l_polar_c[li]:
                _add(c, "hydrophobic contact candidate")
            key = (protein.resname[pi], protein.name[pi])
            if dist <= SALT_BRIDGE_MAX and ((key in _POS and l_charge[li] < 0) or
                                            (key in _NEG and l_charge[li] > 0)):
                _add(c, f"salt-bridge candidate ({protein.name[pi]}...{le}{li + 1} {dist:.2f} A)")
    return sorted(contacts.values(), key=lambda x: x.min_distance)


def _add(contact: Contact, text: str) -> None:
    if text not in contact.candidates:
        contact.candidates.append(text)
