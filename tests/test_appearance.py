"""Appearance panel: object dropdown, style buttons, visibility, remove."""

from __future__ import annotations

from pathlib import Path

import pytest
from conftest import wait_until

from mds.chem.receptor import ReceptorOptions
from mds.engines import VinaEngine
from mds.gui.main_window import MainWindow
from mds.gui.panels.appearance_panel import (MOLECULE_DEFAULT_STYLE, MOLECULE_STYLES,
                                             PROTEIN_DEFAULT_STYLE, PROTEIN_STYLES)
from mds.project import Project

DATA = Path(__file__).parent / "data"
needs_vina = pytest.mark.skipif(not VinaEngine().available(), reason="vina.exe not bundled")


@pytest.fixture
def window(qapp, settings, tmp_path):
    w = MainWindow(settings, use_viewer=False)
    w.controller.set_project(Project.create(tmp_path / "proj"))
    yield w
    w.job_runner.wait_for_done(10000)
    w.close()


def test_panel_starts_empty_with_no_objects(window):
    panel = window.appearance
    assert panel.object_combo.count() == 0
    assert not panel.visible_btn.isEnabled()
    for btn in panel.protein_group.buttons.values():
        assert not btn.isEnabled()


def test_adding_an_object_populates_panel(window):
    window.project_explorer.upsert_object("Proteins", "test_prot",
                                          {"Type": "x"}, visible=True, kind="protein")
    panel = window.appearance
    assert panel.object_combo.count() == 1
    assert panel.object_combo.currentData() == "test_prot"
    assert panel.visible_btn.isEnabled() and panel.visible_btn.isChecked()
    assert panel.protein_group.buttons[PROTEIN_DEFAULT_STYLE].isChecked()
    assert panel.protein_group.buttons[PROTEIN_DEFAULT_STYLE].isEnabled()
    # Molecule styles disabled for a protein
    for btn in panel.molecule_group.buttons.values():
        assert not btn.isEnabled()


def test_switching_to_a_ligand_enables_molecule_styles(window):
    ex = window.project_explorer
    ex.upsert_object("Proteins", "prot", {"Type": "x"}, kind="protein")
    ex.upsert_object("Ligands", "lig", {"Type": "x"}, kind="ligand")
    panel = window.appearance
    panel.object_combo.setCurrentIndex(panel.object_combo.findData("lig"))
    assert panel.molecule_group.buttons[MOLECULE_DEFAULT_STYLE].isChecked()
    assert panel.molecule_group.buttons[MOLECULE_DEFAULT_STYLE].isEnabled()
    for btn in panel.protein_group.buttons.values():
        assert not btn.isEnabled()


def test_clicking_a_style_updates_state_and_viewer(window, monkeypatch):
    calls = []
    monkeypatch.setattr(window, "_viewer_call", lambda *a: calls.append(a))
    # The AppearancePanel captured the OLD _viewer_call during init; replace there too.
    window.appearance._viewer_call = window._viewer_call
    window.project_explorer.upsert_object("Proteins", "prot", {"Type": "x"}, kind="protein")
    panel = window.appearance
    # Click "Ribbon"
    panel.protein_group.buttons["ribbon"].click()
    assert ("set_representation", "prot", "ribbon") in calls
    assert panel._state["prot"]["rep"] == "ribbon"
    # Switch away and back: the ribbon selection persists
    window.project_explorer.upsert_object("Ligands", "lig", {"Type": "x"}, kind="ligand")
    panel.object_combo.setCurrentIndex(panel.object_combo.findData("lig"))
    panel.object_combo.setCurrentIndex(panel.object_combo.findData("prot"))
    assert panel.protein_group.buttons["ribbon"].isChecked()


def test_visibility_button_round_trip(window, monkeypatch):
    calls = []
    monkeypatch.setattr(window, "_viewer_call", lambda *a: calls.append(a))
    window.appearance._viewer_call = window._viewer_call
    window.project_explorer.upsert_object("Proteins", "prot", {"Type": "x"}, kind="protein")
    panel = window.appearance
    panel.visible_btn.setChecked(False)
    assert ("set_visible", "prot", False) in calls
    assert not window.project_explorer.is_visible("prot")
    # Toggle visibility via the Explorer; the panel must reflect it
    window.project_explorer.set_visible("prot", True)
    window.project_explorer.visibility_changed.emit("prot", True)
    assert panel.visible_btn.isChecked()


def test_removing_an_object_updates_the_panel(window):
    ex = window.project_explorer
    ex.upsert_object("Proteins", "a", {"Type": "x"}, kind="protein")
    ex.upsert_object("Proteins", "b", {"Type": "x"}, kind="protein")
    panel = window.appearance
    assert panel.object_combo.count() == 2
    ex.remove_object("a")
    assert panel.object_combo.count() == 1
    assert panel.object_combo.currentData() == "b"


def test_colour_scheme_change(window, monkeypatch):
    calls = []
    monkeypatch.setattr(window, "_viewer_call", lambda *a: calls.append(a))
    window.appearance._viewer_call = window._viewer_call
    window.project_explorer.upsert_object("Proteins", "prot", {"Type": "x"}, kind="protein")
    panel = window.appearance
    panel.colour_combo.setCurrentIndex(panel.colour_combo.findData("spectrum"))
    assert ("set_color_scheme", "prot", "spectrum") in calls
    assert panel._state["prot"]["colour"] == "spectrum"


def test_reset_restores_default(window, monkeypatch):
    calls = []
    monkeypatch.setattr(window, "_viewer_call", lambda *a: calls.append(a))
    window.appearance._viewer_call = window._viewer_call
    window.project_explorer.upsert_object("Proteins", "prot", {"Type": "x"}, kind="protein")
    panel = window.appearance
    panel.protein_group.buttons["surface"].click()
    panel.colour_combo.setCurrentIndex(panel.colour_combo.findData("spectrum"))
    # Reset
    calls.clear()
    for i in range(panel.layout().count()):
        pass  # (layout walk not needed, we test via the method)
    panel._reset()
    assert panel._state["prot"]["rep"] == PROTEIN_DEFAULT_STYLE
    assert panel.protein_group.buttons[PROTEIN_DEFAULT_STYLE].isChecked()
    assert ("set_representation", "prot", PROTEIN_DEFAULT_STYLE) in calls
