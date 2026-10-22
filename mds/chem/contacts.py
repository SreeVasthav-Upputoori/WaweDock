"""Distance-based protein-ligand contact analysis (deliberately conservative).

Everything reported here is a *geometric candidate* derived from distances only.  Angles,
protonation, water mediation and energetics are not evaluated, so none of these entries
is a verified interaction.  Use a dedicated tool (ProLIF/PLIP) and experimental evidence
before interpreting them chemically.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import gemmi
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
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
    arr = np.asarray(xyz,dtype=float).reshape(-1,3)
    polar_c = _polar_carbon_flags(arr,element)
    return _ProteinAtoms(arr, residue, resname, name, element, polar_c)


def _polar_carbon_flags(arr:np.ndarray,element:list[str]) -> np.ndarray:
    polar_c = np.zeros(len(arr), dtype=bool)
    # Spatial buckets avoid a full protein-wide distance scan for every carbon.
    width=1.7
    grid={}
    for i,e in enumerate(element):
        if e in ("N","O","S"):
            key=tuple(np.floor(arr[i]/width).astype(int))
            grid.setdefault(key,[]).append(arr[i])
    for i,e in enumerate(element):
        if e!="C": continue
        cell=np.floor(arr[i]/width).astype(int)
        hit=False
        for x in (-1,0,1):
            for y in (-1,0,1):
                for z in (-1,0,1):
                    for near in grid.get((cell[0]+x,cell[1]+y,cell[2]+z),()):
                        if np.dot(near-arr[i],near-arr[i]) < width*width:
                            hit=True; break
                    if hit: break
                if hit: break
            if hit: break
        polar_c[i]=hit
    return polar_c



def cached_protein_atoms(path:Path) -> _ProteinAtoms:
    path=Path(path).resolve()
    stat=path.stat()
    return _cached_protein_atoms(str(path),stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)


@lru_cache(maxsize=2)
def _cached_protein_atoms(path:str,size:int,modified:int,changed:int) -> _ProteinAtoms:
    return load_protein_atoms(Path(path))


def pose_contacts(pose: Chem.Mol, protein: _ProteinAtoms, cutoff: float = CONTACT_CUTOFF) -> list[Contact]:
    from rdkit import Chem
    lig = Chem.RemoveHs(pose, sanitize=False)
    lxyz = np.array(lig.GetConformer().GetPositions())
    latoms = list(lig.GetAtoms())
    l_el = [a.GetSymbol() for a in latoms]
    l_charge = np.array([a.GetFormalCharge() for a in latoms])
    l_polar_c = np.array([a.GetSymbol() == "C" and any(n.GetSymbol() in ("N", "O", "S")
                                                      for n in a.GetNeighbors()) for a in latoms])

    contacts: dict[str, Contact] = {}
    for start,distances in contact_distance_blocks(protein.xyz,lxyz):
        for local in np.where(distances.min(axis=1) <= cutoff)[0]:
            pi=start+local
            res = protein.residue[pi]
            c = contacts.setdefault(res, Contact(res, float("inf")))
            c.min_distance = min(c.min_distance, float(distances[local].min()))
            for li in np.where(distances[local] <= cutoff)[0]:
                dist = float(distances[local, li])
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


def contact_distance_blocks(protein_xyz:np.ndarray,ligand_xyz:np.ndarray,block_size:int=256):
    """Same Euclidean distances as the full matrix, with bounded protein rows."""
    for start in range(0,len(protein_xyz),block_size):
        delta=protein_xyz[start:start+block_size,None,:]-ligand_xyz[None,:,:]
        yield start,np.linalg.norm(delta,axis=-1)


def _add(contact: Contact, text: str) -> None:
    if text not in contact.candidates:
        contact.candidates.append(text)
