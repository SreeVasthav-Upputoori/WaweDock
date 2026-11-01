"""Docked-pose handling: PDBQT -> RDKit molecules, symmetry-aware RMSD, SDF export."""

from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from rdkit import Chem

import logging
import urllib.request
from pathlib import Path

import numpy as np
# AllChem umbrella is avoided (its rdShapeHelpers DLL is blocked by Smart App Control on
# some Windows machines); AssignBondOrdersFromTemplate is pure Python so we inline a copy.

from mds.chem.report import PreparationError
from mds.engines.base import Pose


def assign_bond_orders_from_template(refmol: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
    """Copy of ``rdkit.Chem.AllChem.AssignBondOrdersFromTemplate`` (pure Python).

    Transfers bond orders from ``refmol`` onto ``mol`` after a substructure match.
    """
    from rdkit.Chem.rdchem import BondType
    from rdkit import Chem
    refmol2 = Chem.Mol(refmol)
    mol2 = Chem.Mol(mol)
    matching = mol2.GetSubstructMatch(refmol2)
    if not matching:
        for b in mol2.GetBonds():
            if b.GetBondType() != BondType.SINGLE:
                b.SetBondType(BondType.SINGLE)
                b.SetIsAromatic(False)
        for b in refmol2.GetBonds():
            b.SetBondType(BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in refmol2.GetAtoms():
            a.SetFormalCharge(0)
        for a in mol2.GetAtoms():
            a.SetFormalCharge(0)
        matching = mol2.GetSubstructMatch(refmol2)
    if not matching:
        raise ValueError("No matching found")
    for b in refmol.GetBonds():
        a1 = matching[b.GetBeginAtomIdx()]
        a2 = matching[b.GetEndAtomIdx()]
        b2 = mol2.GetBondBetweenAtoms(a1, a2)
        if b2 is None:
            raise ValueError("Template bond not found in mol")
        b2.SetBondType(b.GetBondType())
        b2.SetIsAromatic(b.GetIsAromatic())
    for a in refmol.GetAtoms():
        mol2.GetAtomWithIdx(matching[a.GetIdx()]).SetFormalCharge(a.GetFormalCharge())
    try:
        Chem.SanitizeMol(mol2)
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"Could not sanitise after bond-order transfer: {exc}") from exc
    if hasattr(mol, "__sssAtoms"):
        mol2.__sssAtoms = mol.__sssAtoms  # type: ignore[attr-defined]
    return mol2

log = logging.getLogger("mds.chem.poses")


def poses_to_mols(poses: list[Pose], ligand_pdbqt: str) -> list[Chem.Mol]:
    """Convert pose PDBQT blocks to RDKit molecules with correct bond orders.

    Bond orders/aromaticity come from the ``REMARK SMILES`` records Meeko wrote into the
    input ligand; the heavy-atom coordinates come from the docked pose.
    """
    from rdkit import Chem
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
    from rdkit import Chem
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
    from rdkit import Chem
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
    from rdkit import Chem
    p = np.array(Chem.RemoveHs(probe, sanitize=False).GetConformer().GetPositions()).mean(axis=0)
    r = np.array(Chem.RemoveHs(reference, sanitize=False).GetConformer().GetPositions()).mean(axis=0)
    return float(np.linalg.norm(p - r))


# ------------------------------------------------------------ reference ligands
def fetch_ccd_template(resname: str, cache_dir: Path | None = None) -> Chem.Mol:
    """Download the ideal-geometry definition (bond orders!) of a PDB ligand from RCSB."""
    from rdkit import Chem
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
    from rdkit import Chem
    raw = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=True, sanitize=False, proximityBonding=True)
    if raw is None:
        raise PreparationError(f"Could not read ligand coordinates from {ligand_pdb.name}")
    try:
        ref = assign_bond_orders_from_template(template, raw)
    except Exception as exc:  # noqa: BLE001
        raise PreparationError(f"Ligand {ligand_pdb.name} does not match its chemical definition: "
                               f"{exc}") from exc
    return ref


# ---------------------------------------------------------------------- export
def write_poses_sdf(path: Path, mols: list[Chem.Mol], properties: list[dict]) -> None:
    from rdkit import Chem
    writer = Chem.SDWriter(str(path))
    try:
        for mol, props in zip(mols, properties):
            m = Chem.Mol(mol)
            for key, value in props.items():
                m.SetProp(str(key), str(value))
            writer.write(m)
    finally:
        writer.close()


def mcs_rmsd(probe: Chem.Mol, reference: Chem.Mol, *, timeout: int = 10,
             max_matches: int = 200) -> tuple[float | None, int, float]:
    """In-place RMSD over the maximum common substructure of two *different* ligands.

    Returns ``(rmsd, mcs_heavy_atoms, fraction_of_reference)``. Unlike
    :func:`symmetric_rmsd` this does not require identical chemistry, so a docked analogue can
    be compared with a reference ligand on the scaffold they share. Like it, the molecules are
    never re-superposed: the RMSD is measured between the coordinates as docked, which is what
    a binding-mode comparison needs. ``rmsd`` is None when no common substructure of at least
    three atoms is found.
    """
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    p = Chem.RemoveHs(probe, sanitize=False)
    r = Chem.RemoveHs(reference, sanitize=False)
    # Skipping sanitisation keeps odd docked chemistry readable, but leaves ring info unset,
    # which the ring-aware MCS query then trips over. FastFindRings fills it in cheaply.
    for m in (p, r):
        Chem.FastFindRings(m)
    result = rdFMCS.FindMCS([r, p], atomCompare=rdFMCS.AtomCompare.CompareAny,
                            bondCompare=rdFMCS.BondCompare.CompareAny,
                            ringMatchesRingOnly=True, completeRingsOnly=True, timeout=timeout)
    if result.canceled or result.numAtoms < 3:
        # Fall back to element-matched atoms without the ring constraints, which is cheaper.
        result = rdFMCS.FindMCS([r, p], atomCompare=rdFMCS.AtomCompare.CompareElements,
                                bondCompare=rdFMCS.BondCompare.CompareAny, timeout=timeout)
    fraction = round(result.numAtoms / max(r.GetNumAtoms(), 1), 3)
    if result.numAtoms < 3:
        return None, int(result.numAtoms), fraction

    query = getattr(result, "queryMol", None) or Chem.MolFromSmarts(result.smartsString)
    if query is None:
        return None, int(result.numAtoms), fraction
    ref_match = r.GetSubstructMatch(query)
    probe_matches = p.GetSubstructMatches(query, uniquify=False, maxMatches=max_matches)
    if not ref_match or not probe_matches:
        return None, int(result.numAtoms), fraction

    r_xyz = np.array(r.GetConformer().GetPositions())[list(ref_match)]
    p_all = np.array(p.GetConformer().GetPositions())
    # Every embedding of the MCS in the probe is tried, so a symmetric scaffold is not
    # penalised for an arbitrary atom labelling; the reference needs only one embedding
    # because the probe's enumeration already spans every relative mapping.
    best = min(float(np.sqrt(np.mean(np.sum((p_all[list(m)] - r_xyz) ** 2, axis=1))))
               for m in probe_matches)
    return round(best, 3), int(result.numAtoms), fraction
