"""Native molecular context menus driven by a snapshot of the viewer state."""
from __future__ import annotations

from PySide6.QtGui import QActionGroup, QColor
from PySide6.QtWidgets import QColorDialog, QMenu

REPRESENTATIONS = (
    ("Cartoon", "cartoon"), ("Ribbon", "ribbon"), ("Backbone", "backbone"),
    ("Sticks", "sticks"), ("Ball and Stick", "ball-and-stick"),
    ("Spheres", "spheres"), ("Lines", "lines"), ("Molecular Surface", "surface"),
)
COLORS = (("Chain", "chain"), ("Element", "element"), ("Residue", "residue"),
          ("Secondary Structure", "secondary structure"), ("Spectrum", "spectrum"))


def build_context_menu(viewer, state: dict) -> QMenu:
    """Build without executing, allowing action wiring to be tested without WebGL."""
    menu = QMenu(viewer)
    name = state.get("name")

    def choice(parent, text, checked, callback):
        action = parent.addAction(text)
        action.setCheckable(True)
        action.setChecked(checked)
        action.triggered.connect(callback)
        return action

    def pick_color(background=False):
        initial = state.get("background") if background else state.get("customColor")
        color = QColorDialog.getColor(QColor(initial or "#9aa5b1"), viewer,
                                      "Background Color" if background else "Molecular Color")
        if color.isValid():
            if background:
                viewer._call("setBackgroundColor", color.name())
            else:
                viewer._call("setCustomColor", name, color.name())

    if name:
        menu.addSection(str(state.get("target") or name))
        rep = menu.addMenu("Protein Representation" if state.get("kind") == "protein" else "Molecular Representation")
        group = QActionGroup(rep)
        for label, value in REPRESENTATIONS:
            if state.get("kind") != "protein" and value in {"cartoon", "ribbon", "backbone"}:
                continue
            action = choice(rep, label, state.get("representation") == value,
                            lambda _on=False, v=value: viewer.set_representation(name, v))
            group.addAction(action)
        rep.addSeparator()
        rep.addAction("Hide Representation", lambda: viewer.set_visible(name, False))
        colors = menu.addMenu("Color Options")
        group = QActionGroup(colors)
        for label, value in COLORS:
            action = choice(colors, "Color by " + label,
                            not state.get("customColor") and state.get("color") == value,
                            lambda _on=False, v=value: viewer.set_color_scheme(name, v))
            group.addAction(action)
        colors.addSeparator()
        colors.addAction("Custom Color…", pick_color)
        operations = menu.addMenu("Molecular Operations")
        for label, value in (("Hydrogens", "hydrogens"), ("Waters", "waters"), ("Ligands", "ligands")):
            choice(operations, "Show " + label, state.get("components", {}).get(value, True),
                   lambda on, v=value: viewer._call("setComponentVisible", name, v, on))
        choice(operations, "Show Labels", state.get("labels", False),
               lambda on: viewer._call("setObjectLabels", name, on))
        operations.addSeparator()
        operations.addAction("Zoom to Selected Structure", lambda: viewer.zoom_to(name))
        if state.get("editable") and state.get("atomIndex") is not None:
            select = menu.addMenu("Select Clicked…")
            for label, mode in (("Atom", "atom"), ("Residue", "residue"),
                                ("Chain", "chain"), ("Component", "component")):
                select.addAction(label, lambda _on=False, m=mode: viewer._call(
                    "selectContextTarget", name, state["atomIndex"], m))
        else:
            menu.addSection("Read-only overlay — edit the imported structure")

    background = menu.addMenu("Background Settings")
    group = QActionGroup(background)
    for label, value in (("Black", "#000000"), ("White", "#ffffff"), ("Dark Grey", "#1b1d21")):
        group.addAction(choice(background, label, state.get("background") == value,
                               lambda _on=False, v=value: viewer._call("setBackgroundColor", v)))
    background.addAction("Custom Background Color…", lambda: pick_color(True))
    objects = state.get("objects", [])
    if objects:
        visibility = menu.addMenu("Show / Hide Structures")
        for obj in objects:
            choice(visibility, obj["name"], obj["visible"],
                   lambda on, n=obj["name"]: viewer.set_visible(n, on))
    menu.addSeparator()
    selected = bool(state.get("selectionCount"))
    highlighted = bool(state.get("focusCount", state.get("selectionCount")))
    menu.addAction("Focus on Selection", lambda: viewer._call("focusSelection")).setEnabled(highlighted)
    menu.addAction("Clear Selection", viewer.clear_selection).setEnabled(highlighted)
    menu.addAction("Delete Selected Atoms / Residues",
                   lambda: viewer.edit_requested.emit("delete")).setEnabled(selected)
    menu.addAction("Reset Camera", lambda: viewer.zoom_to(None))
    menu.addSeparator()
    menu.addAction("Undo Structure Edit", lambda: viewer.edit_requested.emit("undo"))
    menu.addAction("Redo Structure Edit", lambda: viewer.edit_requested.emit("redo"))
    return menu
