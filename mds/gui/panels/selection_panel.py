"""Viewport selection, reversible structural edits and opt-in component cleanup."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFrame, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QPushButton, QScrollArea,
                               QSizePolicy, QToolButton, QVBoxLayout, QWidget)

if TYPE_CHECKING:
    from mds.gui.controller import Controller


class SelectionPanel(QWidget):
    """Small side panel; buttons are enabled only while something is selected."""

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller

        self.summary = QLabel("Nothing selected.")
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet("font-weight: 600;")
        self.summary.setMinimumWidth(0)

        tip = QLabel(
            "<small>Click selects a residue; double-click selects the whole protein or all waters. "
            "Selections appear yellow. Shift/Ctrl-click toggles. Alt-drag draws a freehand "
            "lasso through the structure's depth; Shift+Alt adds. Delete or right-click "
            "to delete. Esc clears. Edits can be undone.</small>")
        tip.setWordWrap(True)
        tip.setVisible(False)
        help_btn = QToolButton()
        help_btn.setText("Selection shortcuts")
        help_btn.setCheckable(True)
        help_btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        help_btn.setArrowType(Qt.ArrowType.RightArrow)
        help_btn.toggled.connect(tip.setVisible)
        help_btn.toggled.connect(lambda opened: help_btn.setArrowType(
            Qt.ArrowType.DownArrow if opened else Qt.ArrowType.RightArrow))

        self.expand_btn = self._button("Expand to chain", "Grow the selection to every "
                                       "residue in the same chain(s).")
        self.expand_btn.clicked.connect(self._expand)
        self.invert_btn = self._button("Invert selection", "Select every residue NOT currently selected.")
        self.invert_btn.clicked.connect(self._invert)
        self.clear_btn = self._button("Clear selection",
                                      "Deselect everything.")
        self.clear_btn.clicked.connect(self._clear)

        self.delete_btn = self._button("Delete selected", "Remove the selected atoms "
                                       "from the active protein and reload the viewer.")
        self.delete_btn.clicked.connect(lambda: self._edit(keep=False))

        self.keep_btn = self._button("Keep only selected", "Remove everything EXCEPT the "
                                     "selected atoms from the active protein.")
        self.keep_btn.clicked.connect(lambda: self._edit(keep=True))

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll = scroll
        content = QWidget()
        scroll.setWidget(content)
        outer.addWidget(scroll)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        layout.addWidget(self.summary)
        layout.addWidget(help_btn)
        layout.addWidget(tip)
        self.edit_source_btn = self._button("Show editable structure", "Show the imported structure for selection and editing.")
        self.edit_source_btn.clicked.connect(self.c.editing_view_requested.emit)
        layout.addWidget(self.edit_source_btn)
        self.mode = QComboBox()
        for label, mode in (("Residues", "residue"), ("Atoms", "atom"),
                            ("Chains", "chain"), ("Components", "component")):
            self.mode.addItem(label, mode)
        self.mode.setMinimumContentsLength(8)
        self.mode.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.mode.setToolTip(
            "Atoms: select atom centres inside the lasso. Residues / Chains: expand "
            "each hit to its whole residue / chain. The lasso selects through depth.")
        self.mode.currentIndexChanged.connect(
            lambda _: self.c.selection_mode_requested.emit(self.mode.currentData()))
        layout.addWidget(QLabel("Select by"))
        layout.addWidget(self.mode)
        selection_group = QGroupBox("Selection")
        selection_layout = QVBoxLayout(selection_group)
        selection_layout.addWidget(self.expand_btn)
        selection_layout.addWidget(self.invert_btn)
        selection_layout.addWidget(self.clear_btn)
        layout.addWidget(selection_group)
        edit_group = QGroupBox("Edit structure")
        edit_layout = QVBoxLayout(edit_group)
        edit_layout.addWidget(self.delete_btn)
        edit_layout.addWidget(self.keep_btn)
        history = QHBoxLayout()
        self.undo_btn = QPushButton("Undo")
        self.redo_btn = QPushButton("Redo")
        self.undo_btn.setEnabled(False)
        self.redo_btn.setEnabled(False)
        self.undo_btn.clicked.connect(self.c.undo_structure_edit)
        self.redo_btn.clicked.connect(self.c.redo_structure_edit)
        history.addWidget(self.undo_btn)
        history.addWidget(self.redo_btn)
        edit_layout.addLayout(history)
        layout.addWidget(edit_group)

        spatial = QGroupBox("Spatial selection")
        spatial_layout = QVBoxLayout(spatial)
        self.radius = QDoubleSpinBox()
        self.radius.setRange(0.1, 100.0)
        self.radius.setValue(5.0)
        self.radius.setSuffix(" Å")
        self.radius.setToolTip("Union of spheres around the selected atoms, in 3D coordinates.")
        spatial_layout.addWidget(self.radius)
        self.radius_btn = self._button("Select nearby atoms", "Expand the current selection to atoms within the specified radius.")
        self.radius_btn.setEnabled(False)
        self.radius_btn.clicked.connect(lambda: self.c.selection_radius_requested.emit(self.radius.value()))
        spatial_layout.addWidget(self.radius_btn)
        layout.addWidget(spatial)

        cleanup = QGroupBox("Automatic cleanup")
        cleanup_layout = QVBoxLayout(cleanup)
        self.remove_water = QCheckBox("Remove waters")
        self.remove_hetero = QCheckBox("Remove all HETATM")
        self.remove_ligands = QCheckBox("Remove ligands")
        self.retain_selected = QCheckBox("Protect selected atoms")
        self.remove_hetero.setToolTip("Remove all HETATM components, including waters, ligands and ions.")
        self.remove_ligands.setToolTip("Remove non-protein ligands.")
        self.retain_selected.setToolTip("Keep selected atoms and components during automatic cleanup.")
        self.retain_selected.setChecked(True)
        for widget in (self.remove_water, self.remove_hetero, self.remove_ligands):
            cleanup_layout.addWidget(widget)
        cleanup_layout.addWidget(QLabel("Remove chains"))
        self.remove_chains = QLineEdit()
        self.remove_chains.setPlaceholderText("Chain IDs, e.g. B, C")
        cleanup_layout.addWidget(self.remove_chains)
        cleanup_layout.addWidget(self.retain_selected)
        self.cleanup_btn = QPushButton("Apply cleanup")
        self.cleanup_btn.clicked.connect(self._cleanup)
        cleanup_layout.addWidget(self.cleanup_btn)
        layout.addWidget(cleanup)
        layout.addStretch(1)
        # Keep the content within a narrow dock; only vertical scrolling is needed.
        for widget in content.findChildren(QWidget):
            if isinstance(widget, (QPushButton, QComboBox, QLineEdit, QDoubleSpinBox)):
                widget.setMinimumWidth(0)
                widget.setMinimumHeight(28)
                widget.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            elif isinstance(widget, QLabel):
                widget.setWordWrap(True)
                widget.setMinimumWidth(0)
                widget.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        content.setMinimumWidth(0)
        content.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)

        self._update_buttons(empty=True)
        controller.selection_changed.connect(self._on_selection)
        controller.structure_loaded.connect(lambda *_: self._on_selection({}, ""))
        controller.edit_history_changed.connect(self._history)

    def _button(self, label: str, tooltip: str) -> QPushButton:
        b = QPushButton(label)
        b.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        b.setToolTip(tooltip)
        return b

    def _on_selection(self, detail: dict, summary: str) -> None:
        total = sum(len(v) for v in detail.values())
        self.summary.setText(summary or "Nothing selected.")
        self._update_buttons(empty=total == 0)
        self.radius_btn.setEnabled(total > 0)

    def _history(self, undo: bool, redo: bool) -> None:
        self.undo_btn.setEnabled(undo)
        self.redo_btn.setEnabled(redo)

    def _cleanup(self) -> None:
        self.c.cleanup_structure(
            waters=self.remove_water.isChecked(), hetero=self.remove_hetero.isChecked(),
            ligands=self.remove_ligands.isChecked(),
            remove_chains=[s.strip() for s in self.remove_chains.text().split(",") if s.strip()],
            retain_selected=self.retain_selected.isChecked())

    def _update_buttons(self, *, empty: bool) -> None:
        for b in (self.expand_btn, self.invert_btn, self.clear_btn,
                  self.delete_btn, self.keep_btn):
            b.setEnabled(not empty or b is self.invert_btn)

    def _expand(self) -> None:
        self.c.expand_selection_to_chain()

    def _invert(self) -> None:
        self.c.invert_selection()

    def _clear(self) -> None:
        self.c.clear_selection()

    def _edit(self, *, keep: bool) -> None:
        total = sum(len(v) for v in self.c.selection.values())
        if total == 0:
            return
        if keep:
            self.c.keep_only_selected_residues()
        else:
            self.c.delete_selected_residues()
