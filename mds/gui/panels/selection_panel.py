"""Selection panel: shows the current viewer selection and the actions on it.

The researcher's mental model from Discovery Studio:
  * Click in the viewer picks residue(s)
  * The side panel shows a plain-English count ("3 residues in chain A")
  * Big buttons: Expand to chain / Invert / Clear / Delete / Keep only this

Delete rewrites the original protein structure for the active project (so the next
docking preparation uses the edited version) and reloads it in the viewer.  The
original file in ``originals/`` is NEVER modified.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (QFrame, QLabel, QMessageBox, QPushButton, QSizePolicy,
                               QVBoxLayout, QWidget)

from mds.gui.controller import Controller


class SelectionPanel(QWidget):
    """Small side panel; buttons are enabled only while something is selected."""

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller

        self.summary = QLabel("Nothing selected.")
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet("font-size: 11pt;")

        tip = QLabel(
            "<small>Click a residue to select it, Shift-click to add, "
            "Alt-drag to lasso-box, Esc to clear.</small>")
        tip.setWordWrap(True)
        tip.setEnabled(False)

        self.expand_btn = self._button("Expand to chain", "Grow the selection to every "
                                       "residue in the same chain(s).")
        self.expand_btn.clicked.connect(self._expand)
        self.invert_btn = self._button("Invert selection", "Select every residue NOT currently selected.")
        self.invert_btn.clicked.connect(self._invert)
        self.clear_btn = self._button("Clear selection (Esc)",
                                      "Deselect everything.")
        self.clear_btn.clicked.connect(self._clear)

        self.delete_btn = self._button("Delete selected residues", "Remove the selected residues "
                                       "from the active protein and reload the viewer.")
        self.delete_btn.setStyleSheet("background: #fecaca; font-weight: 600;")
        self.delete_btn.clicked.connect(lambda: self._edit(keep=False))

        self.keep_btn = self._button("Keep only selected", "Remove everything EXCEPT the "
                                     "selected residues from the active protein.")
        self.keep_btn.setStyleSheet("background: #fde68a; font-weight: 600;")
        self.keep_btn.clicked.connect(lambda: self._edit(keep=True))

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)
        layout.addWidget(self.summary)
        layout.addWidget(tip)
        layout.addWidget(_hline())
        layout.addWidget(self.expand_btn)
        layout.addWidget(self.invert_btn)
        layout.addWidget(self.clear_btn)
        layout.addWidget(_hline())
        layout.addWidget(self.delete_btn)
        layout.addWidget(self.keep_btn)
        layout.addStretch(1)

        self._update_buttons(empty=True)
        controller.selection_changed.connect(self._on_selection)
        controller.structure_loaded.connect(lambda *_: self._on_selection({}, ""))

    def _button(self, label: str, tooltip: str) -> QPushButton:
        b = QPushButton(label)
        b.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        b.setToolTip(tooltip)
        return b

    def _on_selection(self, detail: dict, summary: str) -> None:
        total = sum(len(v) for v in detail.values())
        self.summary.setText(summary or "Nothing selected.")
        self._update_buttons(empty=total == 0)

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
        verb, action = ("Keep only", "keep ONLY") if keep else ("Delete", "DELETE")
        msg = (f"{action} the {total} selected residue(s)?\n\n"
               "The active protein will be re-saved into the project's prepared/ folder "
               "and reloaded. The original file in originals/ is not modified.\n\n"
               "You will need to run Prepare Receptor again before docking.")
        answer = QMessageBox.question(self, f"{verb} residues", msg)
        if answer != QMessageBox.StandardButton.Yes:
            return
        if keep:
            self.c.keep_only_selected_residues()
        else:
            self.c.delete_selected_residues()


def _hline() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line
