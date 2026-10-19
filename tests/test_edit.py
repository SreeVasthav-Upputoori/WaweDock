"""Unit tests for structural edits (delete/keep residues)."""
from __future__ import annotations

from pathlib import Path
import gemmi
import pytest
from mds.chem.edit import delete_residues, keep_only_residues

DATA = Path(__file__).parent / "data"


def _count_residues(path: Path) -> dict[str, int]:
    st = gemmi.read_structure(str(path))
    out: dict[str, int] = {}
    for ch in st[0]:
        out[ch.name] = out.get(ch.name, 0) + len(ch)
    return out


def test_delete_removes_only_the_named_residues(tmp_path):
    dst = tmp_path / "out.pdb"
    n = delete_residues(DATA / "3PTB.pdb", dst, [("A", 1), {"chain": "A", "resi": 480}])
    assert n == 2 and dst.exists()
    st = gemmi.read_structure(str(dst))
    for ch in st[0]:
        for res in ch:
            assert (ch.name, res.seqid.num) not in {("A", 1), ("A", 480)}


def test_delete_non_matching_residues_is_a_no_op(tmp_path):
    dst = tmp_path / "out.pdb"
    original = _count_residues(DATA / "3PTB.pdb")
    n = delete_residues(DATA / "3PTB.pdb", dst, [("Z", 999)])
    assert n == 0
    assert _count_residues(dst) == original


def test_keep_only_shrinks_to_the_named_residues(tmp_path):
    dst = tmp_path / "out.pdb"
    n = keep_only_residues(DATA / "3PTB.pdb", dst, [("A", 57), ("A", 102), ("A", 195)])
    assert n > 0
    st = gemmi.read_structure(str(dst))
    kept = {(ch.name, r.seqid.num) for ch in st[0] for r in ch}
    assert kept == {("A", 57), ("A", 102), ("A", 195)}


def test_delete_empty_list_passes_through(tmp_path):
    dst = tmp_path / "out.pdb"
    n = delete_residues(DATA / "3PTB.pdb", dst, [])
    assert n == 0 and dst.exists()


def test_delete_rejects_bad_input(tmp_path):
    with pytest.raises(TypeError):
        delete_residues(DATA / "3PTB.pdb", tmp_path / "x.pdb", ["chain:resi"])
