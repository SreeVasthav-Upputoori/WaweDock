"""Selection pipeline: viewer event -> controller -> delete/keep edit -> reload."""
from __future__ import annotations

from pathlib import Path

import pytest
from conftest import wait_until

from mds.gui.main_window import MainWindow
from mds.project import Project

DATA = Path(__file__).parent / "data"


@pytest.fixture
def window(qapp, settings, tmp_path):
    w = MainWindow(settings, use_viewer=False)
    w.controller.set_project(Project.create(tmp_path / "proj"))
    yield w
    w.close()


def test_selection_signal_round_trip(window):
    c = window.controller
    received = []
    c.selection_changed.connect(lambda d, s: received.append((d, s)))

    c.set_selection({"3PTB": [{"chain": "A", "resi": 15}]}, "1 residue(s)")
    assert c.selection == {"3PTB": [{"chain": "A", "resi": 15}]}
    assert received[-1][1] == "1 residue(s)"

    c.clear_selection()
    assert c.selection == {}
    assert received[-1] == ({}, "")


def test_delete_requires_a_loaded_structure(window):
    c = window.controller
    errors = []
    c.error.connect(lambda t, m: errors.append((t, m)))
    c.set_selection({"foo": [{"chain": "A", "resi": 1}]}, "")
    c.delete_selected_residues()
    assert errors and "Select residues on the loaded protein" in errors[0][1]


def test_delete_selected_residues_rewrites_pdb(window):
    c = window.controller
    c.load_structure(DATA / "3PTB.pdb")
    assert c.source_pdb is not None
    name = c.source_pdb.stem
    before_atoms = _count_atoms(c.source_pdb)

    c.set_selection({name: [{"chain": "A", "resi": 1}, {"chain": "A", "resi": 2}]},
                    "2 residue(s)")
    assert wait_until(lambda: True, 0.05)   # flush event loop
    c.delete_selected_residues()

    # new file written
    assert c.source_pdb.name.endswith("_edited.pdb")
    after_atoms = _count_atoms(c.source_pdb)
    assert after_atoms < before_atoms
    # selection cleared
    assert c.selection == {}
    # receptor state reset (user must re-prepare)
    assert c.receptor is None


def test_keep_only_shrinks_structure(window):
    c = window.controller
    c.load_structure(DATA / "3PTB.pdb")
    name = c.source_pdb.stem
    c.set_selection({name: [{"chain": "A", "resi": 57}, {"chain": "A", "resi": 102},
                            {"chain": "A", "resi": 195}]}, "")
    c.keep_only_selected_residues()
    assert c.source_pdb.name.endswith("_edited.pdb")
    assert _count_atoms(c.source_pdb) < 100   # catalytic triad only


def _count_atoms(path: Path) -> int:
    return sum(1 for l in path.read_text().splitlines() if l.startswith(("ATOM", "HETATM")))
