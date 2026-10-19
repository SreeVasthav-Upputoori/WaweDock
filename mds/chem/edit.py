"""Structural edits a researcher commonly needs before docking: remove residues, drop
chains, keep only a region.  All edits read a PDB/mmCIF file with gemmi, apply the
change, and write a new file next to the original.  The input file is never modified.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi


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
    # Minimal PDB preserves atom records including HETATM entries; adequate for docking prep.
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
