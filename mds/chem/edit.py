"""Structural edits a researcher commonly needs before docking: remove residues, drop
chains, keep only a region.  All edits read a PDB/mmCIF file with gemmi, apply the
change, and write a new file next to the original.  The input file is never modified.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi


def atom_records(st: gemmi.Structure) -> list[dict]:
    """Stable atom identities shared by the viewer and structural edit operations."""
    records = []
    for mi, model in enumerate(st):
        for ci, chain in enumerate(model):
            for ri, res in enumerate(chain):
                for atom in res:
                    records.append({
                        "atom_index": len(records), "model": mi, "chain": chain.name,
                        "chain_index": ci, "residue_index": ri,
                        "resi": res.seqid.num, "icode": res.seqid.icode.strip(),
                        "resn": res.name, "atom": atom.name,
                        "altloc": atom.altloc.replace("\0", ""),
                        "x": atom.pos.x, "y": atom.pos.y, "z": atom.pos.z,
                    })
    return records


def selected_atom_indices(st: gemmi.Structure, selections: list[dict]) -> set[int]:
    """Resolve exact viewport atoms, or legacy residue selectors, without serial ambiguity."""
    records = atom_records(st)
    selected = set()
    for item in selections:
        if "atom_index" in item:
            index = int(item["atom_index"])
            if index < 0 or index >= len(records):
                raise ValueError("The atom selection is stale; select the atoms again.")
            record = records[index]
            for key in ("model", "chain", "chain_index", "residue_index", "resi", "icode", "resn", "atom", "altloc"):
                if key in item and item[key] != record[key]:
                    raise ValueError("The atom selection no longer matches the structure.")
            selected.add(index)
        else:
            for record in records:
                if record["chain"] == item["chain"] and record["resi"] == int(item["resi"]):
                    if all(key not in item or item[key] == record[key]
                           for key in ("model", "icode", "resn")):
                        selected.add(record["atom_index"])
    return selected


def edit_atoms(src: Path, dst: Path, selections: list[dict], *, keep: bool = False) -> int:
    """Delete/retain exact atoms, pruning empty containers. Return the removed atom count."""
    st = gemmi.read_structure(str(src))
    selected = selected_atom_indices(st, selections)
    index, removed = 0, 0
    for model in st:
        empty_chains = []
        for ci, chain in enumerate(model):
            empty_residues = []
            for ri, res in enumerate(chain):
                victims = []
                for ai, _ in enumerate(res):
                    if (index in selected) != keep:
                        victims.append(ai)
                    index += 1
                for ai in reversed(victims):
                    del res[ai]
                removed += len(victims)
                if not len(res):
                    empty_residues.append(ri)
            for ri in reversed(empty_residues):
                del chain[ri]
            if not len(chain):
                empty_chains.append(ci)
        for ci in reversed(empty_chains):
            del model[ci]
    # CONECT records can reference deleted atoms; regenerate connectivity on display.
    st.clear_conect()
    _write_pdb(st, dst)
    return removed


def cleanup_selection(src: Path, *, waters=False, hetero=False, ligands=False,
                      remove_chains=(), retain=()) -> list[dict]:
    """Find cleanup victims while protecting every explicitly retained atom."""
    st = gemmi.read_structure(str(src))
    protected = selected_atom_indices(st, list(retain))
    victims, index = [], 0
    records = atom_records(st)
    for model in st:
        for chain in model:
            for res in chain:
                info = gemmi.find_tabulated_residue(res.name)
                water = res.is_water()
                ligand = not water and not (info and info.is_amino_acid()) and len(res) > 1
                remove = ((waters and water) or (hetero and res.het_flag == "H")
                          or (ligands and ligand) or chain.name in remove_chains)
                for _ in res:
                    if remove and index not in protected:
                        victims.append(records[index])
                    index += 1
    return victims


@dataclass(frozen=True)
class Residue:
    chain: str
    seqid: int


def _parse_residues(residues: list[Residue] | list[tuple[str, int]] | list[dict]) -> set[tuple[str, int]]:
    """Normalise any supported form of residue list into ``{(chain, seqid)}``."""
    out: set[tuple[str, int]] = set()
    for item in residues:
        if isinstance(item, Residue):
            out.add((item.chain, int(item.seqid)))
        elif isinstance(item, dict):
            out.add((str(item["chain"]), int(item["resi"])))
        elif isinstance(item, (tuple, list)) and len(item) == 2:
            out.add((str(item[0]), int(item[1])))
        else:
            raise TypeError(f"Cannot interpret residue: {item!r}")
    return out


def delete_residues(src: Path, dst: Path, residues) -> int:
    """Write ``src`` to ``dst`` with the given residues removed.  Returns # removed."""
    victims = _parse_residues(residues)
    st = gemmi.read_structure(str(src))
    removed = 0
    for model in st:
        for ch in list(model):
            drop_indices: list[int] = []
            for i, res in enumerate(ch):
                if (ch.name, res.seqid.num) in victims:
                    drop_indices.append(i)
                    removed += 1
            for i in reversed(drop_indices):
                del ch[i]
            if len(ch) == 0:
                model.remove_chain(ch.name)
    _write_pdb(st, dst)
    return removed


def keep_only_residues(src: Path, dst: Path, residues) -> int:
    """Write ``dst`` containing ONLY the residues given (opposite of ``delete_residues``)."""
    keepers = _parse_residues(residues)
    st = gemmi.read_structure(str(src))
    removed = 0
    for model in st:
        for ch in list(model):
            drop_indices: list[int] = []
            for i, res in enumerate(ch):
                if (ch.name, res.seqid.num) not in keepers:
                    drop_indices.append(i)
                    removed += 1
            for i in reversed(drop_indices):
                del ch[i]
            if len(ch) == 0:
                model.remove_chain(ch.name)
    _write_pdb(st, dst)
    return removed


def _write_pdb(st: gemmi.Structure, dst: Path) -> None:
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.suffix.lower() in (".cif", ".mmcif"):
        st.make_mmcif_document().write_file(str(dst))
    else:
        st.write_minimal_pdb(str(dst))


def summarise_selection(residues) -> str:
    """Human-readable summary like ``'3 waters, 1 HETATM (STI A:201), 2 residues in chain B'``."""
    parsed = _parse_residues(residues)
    if not parsed:
        return "0 residues"
    chains: dict[str, int] = {}
    for ch, _ in parsed:
        chains[ch] = chains.get(ch, 0) + 1
    parts = [f"{n} residue(s) in chain {ch}" for ch, n in sorted(chains.items())]
    return "; ".join(parts) if parts else "0 residues"
