"""Docked-pose handling: PDBQT -> RDKit molecules, symmetry-aware RMSD, SDF export."""

from __future__ import annotations

import logging
import urllib.request
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from mds.chem.report import PreparationError
from mds.engines.base import Pose

log = logging.getLogger("mds.chem.poses")


def poses_to_mols(poses: list[Pose], ligand_pdbqt: str) -> list[Chem.Mol]:
    """Convert pose PDBQT blocks to RDKit molecules with correct bond orders.

    Bond orders/aromaticity come from the ``REMARK SMILES`` records Meeko wrote into the
    input ligand; the heavy-atom coordinates come from the docked pose.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate

    header = [l for l in ligand_pdbqt.splitlines()
              if l.startswith(("REMARK SMILES", "REMARK H PARENT"))]
    mols: list[Chem.Mol] = []
    for pose in poses:
        body = [l for l in pose.pdbqt.splitlines() if not l.startswith(("REMARK SMILES", "REMARK H PARENT"))]
        text = "MODEL 1\n" + "\n".join(header + body) + "\nENDMDL\n"
        try:
            pm = PDBQTMolecule(text, is_dlg=False, skip_typing=True)
            created = RDKitMolCreate.from_pdbqt_mol(pm)
            mol = created[0] if created else None
        except Exception as exc:  # noqa: BLE001
            log.warning("Pose %d could not be converted to an RDKit molecule: %s", pose.rank, exc)
            mol = None
        if mol is None:
            raise PreparationError(f"Could not rebuild molecule for pose {pose.rank}")
        mol.SetProp("_Name", f"pose_{pose.rank}")
        mol.SetDoubleProp("Score_kcal_per_mol", pose.score)
        mols.append(mol)
    return mols


def _graph_only(mol: Chem.Mol) -> Chem.Mol:
    """Heavy-atom skeleton with identical-bond types, for bond-order-agnostic matching."""
    m = Chem.RemoveHs(mol, sanitize=False)
    rw = Chem.RWMol()
    for a in m.GetAtoms():
        rw.AddAtom(Chem.Atom(a.GetAtomicNum()))
    for b in m.GetBonds():
        rw.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    return rw.GetMol()


def symmetric_rmsd(probe: Chem.Mol, reference: Chem.Mol, max_matches: int = 20000) -> float:
    """Heavy-atom RMSD without re-alignment (as for docking validation).

    All graph automorphisms are tried, so symmetric groups (phenyl rings, carboxylates)
    are not penalised for an equivalent atom labelling.
    """
    p_graph, r_graph = _graph_only(probe), _graph_only(reference)
    if p_graph.GetNumAtoms() != r_graph.GetNumAtoms():
        raise ValueError(f"Heavy-atom counts differ ({p_graph.GetNumAtoms()} vs "
                         f"{r_graph.GetNumAtoms()}); cannot compute RMSD")
    p_xyz = np.array(Chem.RemoveHs(probe, sanitize=False).GetConformer().GetPositions())
    r_xyz = np.array(Chem.RemoveHs(reference, sanitize=False).GetConformer().GetPositions())
    matches = p_graph.GetSubstructMatches(r_graph, uniquify=False, useChirality=False,
                                          maxMatches=max_matches)
    if not matches:
        raise ValueError("Probe and reference have different connectivity; cannot compute RMSD")
    best = min(float(np.sqrt(np.mean(np.sum((p_xyz[list(m)] - r_xyz) ** 2, axis=1))))
               for m in matches)
    return best


def centroid_shift(probe: Chem.Mol, reference: Chem.Mol) -> float:
    p = np.array(Chem.RemoveHs(probe, sanitize=False).GetConformer().GetPositions()).mean(axis=0)
    r = np.array(Chem.RemoveHs(reference, sanitize=False).GetConformer().GetPositions()).mean(axis=0)
    return float(np.linalg.norm(p - r))


# ------------------------------------------------------------ reference ligands
def fetch_ccd_template(resname: str, cache_dir: Path | None = None) -> Chem.Mol:
    """Download the ideal-geometry definition (bond orders!) of a PDB ligand from RCSB."""
    resname = resname.upper()
    cache = (cache_dir / f"{resname}_ideal.sdf") if cache_dir else None
    if cache and cache.exists():
        text = cache.read_text()
    else:
        url = f"https://files.rcsb.org/ligands/download/{resname}_ideal.sdf"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:  # noqa: S310 — fixed https host
                text = resp.read().decode()
        except Exception as exc:  # noqa: BLE001
            raise PreparationError(f"Could not download definition of ligand {resname}: {exc}") from exc
        if cache:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(text)
    mol = Chem.MolFromMolBlock(text, removeHs=False)
    if mol is not None:
        mol = Chem.RemoveAllHs(mol)
    if mol is None:
        raise PreparationError(f"RCSB ligand definition for {resname} could not be parsed")
    return mol


def reference_from_pdb_ligand(ligand_pdb: Path, template: Chem.Mol) -> Chem.Mol:
    """Crystal-coordinate ligand with correct bond orders, using ``template`` for chemistry."""
    raw = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=True, sanitize=False, proximityBonding=True)
    if raw is None:
        raise PreparationError(f"Could not read ligand coordinates from {ligand_pdb.name}")
    try:
        ref = AllChem.AssignBondOrdersFromTemplate(template, raw)
    except Exception as exc:  # noqa: BLE001
        raise PreparationError(f"Ligand {ligand_pdb.name} does not match its chemical definition: "
                               f"{exc}") from exc
    return ref


# ---------------------------------------------------------------------- export
def write_poses_sdf(path: Path, mols: list[Chem.Mol], properties: list[dict]) -> None:
    writer = Chem.SDWriter(str(path))
    try:
        for mol, props in zip(mols, properties):
            m = Chem.Mol(mol)
            for key, value in props.items():
                m.SetProp(str(key), str(value))
            writer.write(m)
    finally:
        writer.close()

