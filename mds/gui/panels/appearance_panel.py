"""Appearance panel: pick representation, colour and visibility for a workspace object.

Design goals (keep it simple):
* The active object is picked from a dropdown; the current style and colour are shown
  with a checked push-button and a dropdown respectively — no hidden menus.
* Four protein-oriented styles + four small-molecule styles, labelled with everyday names
  (``Ball & Stick`` instead of ``ball-and-stick``). Styles not meaningful for an object's
  kind are disabled.
* Visibility is a prominent toggle; zoom and remove are big buttons at the bottom.
* The panel stays in sync with the Project Explorer: when objects appear or disappear
  elsewhere, or the user changes visibility from the tree, this panel updates.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QButtonGroup, QComboBox, QFrame, QGridLayout, QHBoxLayout,
                               QLabel, QPushButton, QSizePolicy, QVBoxLayout, QWidget)

from mds.gui.panels.project_explorer import ProjectExplorer

# (label shown to the user, backend name)
PROTEIN_STYLES = [("Cartoon", "cartoon"), ("Ribbon", "ribbon"),
                  ("Backbone", "backbone"), ("Surface", "surface")]
# "&&" escapes Qt's keyboard-mnemonic parsing so the ampersand is rendered as a character.
MOLECULE_STYLES = [("Sticks", "sticks"), ("Ball && Stick", "ball-and-stick"),
                   ("Spheres", "spheres"), ("Lines", "lines")]
COLOUR_SCHEMES = [("Chain", "chain"), ("Element", "element"), ("Residue", "residue"),
                  ("Secondary structure", "secondary structure"),
                  ("Spectrum", "spectrum"), ("Object colour", "object")]

PROTEIN_DEFAULT_STYLE = "cartoon"
MOLECULE_DEFAULT_STYLE = "sticks"
PROTEIN_DEFAULT_COLOUR = "chain"
MOLECULE_DEFAULT_COLOUR = "element"

# Only objects whose kind is in this set are listed in the Appearance panel.  The docking
# box has its own show/hide; categories don't need appearance controls.
VISIBLE_KINDS = ("protein", "ligand", "result")


class AppearancePanel(QWidget):
    def __init__(self, explorer: ProjectExplorer, viewer_call: Callable[..., None],
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._explorer = explorer
        self._viewer_call = viewer_call
        # Per-object state: {name: {"rep": str, "colour": str}}.  Updated whenever the user
        # changes a control so the dropdown switch restores exactly what they chose last.
        self._state: dict[str, dict[str, str]] = {}
        self._applying = False       # set while we drive controls programmatically

        self._build()
        explorer.objects_changed.connect(self._refresh_objects)
        explorer.visibility_changed.connect(self._on_external_visibility)
        self._refresh_objects()

    # ------------------------------------------------------------------ layout
    def _build(self) -> None:
        self.object_combo = QComboBox()
        self.object_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        # Elide (not truncate silently) if the dock is narrow, and expand to fill the row.
        self.object_combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        self.object_combo.setMinimumContentsLength(16)
        self.object_combo.currentIndexChanged.connect(self._on_object_changed)
        self.object_combo.setToolTip("Which object in the workspace to style")

        self.visible_btn = QPushButton("Visible")
        self.visible_btn.setCheckable(True)
        self.visible_btn.setChecked(True)
        self.visible_btn.toggled.connect(self._on_visibility_toggled)
        self.visible_btn.setToolTip("Show or hide the selected object in the 3D viewer")

        self.protein_group = self._style_group("Protein style", PROTEIN_STYLES, self._set_rep)
        self.molecule_group = self._style_group("Molecule style", MOLECULE_STYLES, self._set_rep)

        self.colour_combo = QComboBox()
        for label, name in COLOUR_SCHEMES:
            self.colour_combo.addItem(label, name)
        self.colour_combo.currentIndexChanged.connect(self._on_colour_changed)

        reset_btn = QPushButton("Reset to default")
        reset_btn.setToolTip("Return this object to its default style and colour.")
        reset_btn.clicked.connect(self._reset)
        zoom_btn = QPushButton("Zoom to object")
        zoom_btn.clicked.connect(self._zoom)
        remove_btn = QPushButton("Remove from workspace")
        remove_btn.setToolTip("Remove this object from the 3D view (does not delete the file)")
        remove_btn.clicked.connect(self._remove)

        head = QHBoxLayout()
        head.addWidget(QLabel("Object:"))
        head.addWidget(self.object_combo, 1)
        head.addWidget(self.visible_btn)

        colour_row = QHBoxLayout()
        colour_row.addWidget(QLabel("Colour by:"))
        colour_row.addWidget(self.colour_combo, 1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(10)
        layout.addLayout(head)
        layout.addWidget(_hline())
        layout.addWidget(self.protein_group)
        layout.addWidget(self.molecule_group)
        layout.addWidget(_hline())
        layout.addLayout(colour_row)
        layout.addWidget(_hline())
        layout.addWidget(reset_btn)
        layout.addWidget(zoom_btn)
        layout.addWidget(remove_btn)
        layout.addStretch(1)

    def _style_group(self, title: str, items: list[tuple[str, str]],
                     slot: Callable[[str], None]) -> QFrame:
        """A labelled frame with an exclusive 2x2 grid of checkable push-buttons."""
        frame = QFrame()
        frame.setFrameShape(QFrame.Shape.NoFrame)
        v = QVBoxLayout(frame)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        v.addWidget(QLabel(f"<b>{title}</b>"))
        grid = QGridLayout()
        grid.setSpacing(4)
        group = QButtonGroup(frame)
        group.setExclusive(True)
        buttons: dict[str, QPushButton] = {}
        for i, (label, name) in enumerate(items):
            b = QPushButton(label)
            b.setCheckable(True)
            b.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            b.clicked.connect(lambda _c=False, n=name: slot(n))
            group.addButton(b)
            buttons[name] = b
            grid.addWidget(b, i // 2, i % 2)
        v.addLayout(grid)
        frame.buttons = buttons      # type: ignore[attr-defined] — for easy lookup later
        return frame

    # -------------------------------------------------------------- object sync
    def _refresh_objects(self) -> None:
        """Rebuild the object dropdown from the Project Explorer, keeping the selection."""
        self._applying = True
        current = self.object_combo.currentData()
        self.object_combo.clear()
        for name in self._explorer.all_object_names(VISIBLE_KINDS):
            self.object_combo.addItem(name, name)
            self._state.setdefault(name, self._default_state_for(name))
        # Drop state for objects that no longer exist.
        existing = set(self._explorer.all_object_names(VISIBLE_KINDS))
        for stale in [n for n in self._state if n not in existing]:
            del self._state[stale]
        idx = self.object_combo.findData(current) if current else -1
        if idx < 0 and self.object_combo.count():
            idx = 0
        self.object_combo.setCurrentIndex(idx)
        self._applying = False
        self._apply_controls()

    def _default_state_for(self, name: str) -> dict[str, str]:
        kind = self._explorer.kind_of(name)
        if kind == "protein":
            return {"rep": PROTEIN_DEFAULT_STYLE, "colour": PROTEIN_DEFAULT_COLOUR}
        return {"rep": MOLECULE_DEFAULT_STYLE, "colour": MOLECULE_DEFAULT_COLOUR}

    def _active_name(self) -> str:
        return self.object_combo.currentData() or ""

    def _apply_controls(self) -> None:
        """Reflect the active object's current state in all controls (no viewer calls)."""
        name = self._active_name()
        has_object = bool(name)
        self.visible_btn.setEnabled(has_object)
        self.colour_combo.setEnabled(has_object)
        for frame in (self.protein_group, self.molecule_group):
            for btn in frame.buttons.values():       # type: ignore[attr-defined]
                btn.setEnabled(False)
                btn.setChecked(False)
        if not has_object:
            return
        self._applying = True
        state = self._state.setdefault(name, self._default_state_for(name))
        kind = self._explorer.kind_of(name)
        # Protein controls are only meaningful for proteins; small molecules the other way.
        for frame, active_for in ((self.protein_group, "protein"),
                                  (self.molecule_group, "ligand")):
            enable = (kind == active_for) or (kind == "result" and active_for == "ligand")
            for name_, btn in frame.buttons.items():   # type: ignore[attr-defined]
                btn.setEnabled(enable)
                btn.setChecked(enable and state["rep"] == name_)
        idx = self.colour_combo.findData(state["colour"])
        self.colour_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.visible_btn.setChecked(self._explorer.is_visible(name))
        self._applying = False

    def _on_object_changed(self) -> None:
        if not self._applying:
            self._apply_controls()

    def _on_external_visibility(self, name: str, visible: bool) -> None:
        if name == self._active_name():
            was = self._applying
            self._applying = True
            self.visible_btn.setChecked(visible)
            self._applying = was

    # ------------------------------------------------------------ user actions
    def _on_visibility_toggled(self, visible: bool) -> None:
        if self._applying:
            return
        name = self._active_name()
        if not name:
            return
        self._explorer.set_visible(name, visible)
        self._viewer_call("set_visible", name, visible)

    def _set_rep(self, rep: str) -> None:
        if self._applying:
            return
        name = self._active_name()
        if not name:
            return
        self._state.setdefault(name, {})["rep"] = rep
        self._viewer_call("set_representation", name, rep)

    def _on_colour_changed(self) -> None:
        if self._applying:
            return
        name = self._active_name()
        colour = self.colour_combo.currentData()
        if not name or not colour:
            return
        self._state.setdefault(name, {})["colour"] = colour
        self._viewer_call("set_color_scheme", name, colour)

    def _zoom(self) -> None:
        name = self._active_name()
        if name:
            self._viewer_call("zoom_to", name)

    def _reset(self) -> None:
        name = self._active_name()
        if not name:
            return
        self._state[name] = self._default_state_for(name)
        self._apply_controls()
        self._viewer_call("set_representation", name, self._state[name]["rep"])
        self._viewer_call("set_color_scheme", name, self._state[name]["colour"])

    def _remove(self) -> None:
        name = self._active_name()
        if name:
            # The Explorer's remove_requested signal is wired to the main window's remove flow,
            # which also removes the object from the viewer and the explorer.
            self._explorer.remove_requested.emit(name)


def _hline() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line
