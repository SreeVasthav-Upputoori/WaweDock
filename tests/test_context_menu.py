"""Native context menu actions target the supplied molecular state."""
from PySide6.QtCore import Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QWidget
import pytest

from mds.render.context_menu import build_context_menu


class Viewer(QWidget):
    edit_requested = Signal(str)

    def __init__(self):
        super().__init__()
        self.calls = []

    def _call(self, *args):
        self.calls.append(args)

    def set_representation(self, *args): self._call("representation", *args)
    def set_color_scheme(self, *args): self._call("color", *args)
    def set_visible(self, *args): self._call("visible", *args)
    def zoom_to(self, *args): self._call("zoom", *args)
    def clear_selection(self): self._call("clear")
    def set_selection_mode(self, *args): self._call("mode", *args)


def find(menu, text):
    for action in menu.actions():
        if action.text() == text:
            return action
        if action.menu():
            found = find(action.menu(), text)
            if found is not None:
                return found
    return None


@pytest.fixture
def viewer(qapp):
    widget = Viewer()
    yield widget
    widget.close()


def state(**changes):
    data = dict(name="protein", kind="protein", target="ALA25:A CA (C)",
                representation="cartoon", color="chain", editable=True, atomIndex=23,
                components=dict(hydrogens=True, waters=True, ligands=True),
                selectionCount=5, selectionMode="residue", background="#ffffff",
                objects=[dict(name="protein", visible=True), dict(name="pose", visible=False)])
    data.update(changes)
    return data


def test_representation_color_visibility_and_components(viewer):
    menu = build_context_menu(viewer, state())
    assert find(menu, "Cartoon").isChecked()
    find(menu, "Ball and Stick").trigger()
    find(menu, "Color by Residue").trigger()
    find(menu, "Show Waters").trigger()
    find(menu, "Hide Representation").trigger()
    find(menu, "pose").trigger()
    assert viewer.calls == [("representation", "protein", "ball-and-stick"),
                            ("color", "protein", "residue"),
                            ("setComponentVisible", "protein", "waters", False),
                            ("visible", "protein", False), ("visible", "pose", True)]


def test_selection_operations_and_read_only_target(viewer):
    edits = []
    viewer.edit_requested.connect(edits.append)
    menu = build_context_menu(viewer, state())
    find(menu, "Atom").trigger()
    find(menu, "Focus on Selection").trigger()
    find(menu, "Delete Selected Atoms / Residues").trigger()
    assert viewer.calls == [("selectContextTarget", "protein", 23, "atom"), ("focusSelection",)]
    assert edits == ["delete"]
    read_only = build_context_menu(viewer, state(editable=False, atomIndex=None, selectionCount=0))
    assert find(read_only, "Select Clicked…") is None
    assert not find(read_only, "Delete Selected Atoms / Residues").isEnabled()


def test_empty_viewport_has_background_without_molecular_actions(viewer):
    menu = build_context_menu(viewer, dict(background="#ffffff", selectionCount=0))
    assert find(menu, "Protein Representation") is None
    assert find(menu, "White").isChecked()
    find(menu, "Dark Grey").trigger()
    assert viewer.calls == [("setBackgroundColor", "#1b1d21")]
    assert not find(menu, "Focus on Selection").isEnabled()


def test_ligand_menu_uses_molecular_representations(viewer):
    menu = build_context_menu(viewer, state(kind="ligand", representation="sticks"))
    assert find(menu, "Molecular Representation") is not None
    assert find(menu, "Cartoon") is None
    assert find(menu, "Sticks").isChecked()


def test_read_only_visual_selection_can_be_focused_and_cleared(viewer):
    menu = build_context_menu(viewer, state(editable=False, atomIndex=None,
                                          selectionCount=0, focusCount=3))
    assert find(menu, "Focus on Selection").isEnabled()
    assert find(menu, "Clear Selection").isEnabled()
    assert not find(menu, "Delete Selected Atoms / Residues").isEnabled()


def test_color_picker_applies_only_valid_colors(viewer, monkeypatch):
    monkeypatch.setattr("mds.render.context_menu.QColorDialog.getColor", lambda *args: QColor("#123456"))
    menu = build_context_menu(viewer, state())
    find(menu, "Custom Color…").trigger()
    find(menu, "Custom Background Color…").trigger()
    assert viewer.calls == [("setCustomColor", "protein", "#123456"),
                            ("setBackgroundColor", "#123456")]
    monkeypatch.setattr("mds.render.context_menu.QColorDialog.getColor", lambda *args: QColor())
    find(menu, "Custom Color…").trigger()
    assert len(viewer.calls) == 2
