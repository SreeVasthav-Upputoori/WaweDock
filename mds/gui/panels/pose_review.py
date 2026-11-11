"""Compact controls for fixed-reference, pose-by-pose visual comparison.

Everything sits on one row. The two display toggles live behind a *Display* button rather than
as permanent check boxes: they are set once and then left alone, so they do not earn standing
space above a table the user reads constantly. The no-reference hint occupies a second row that
collapses to nothing as soon as a reference exists, which is the only state where it helps.
"""
from PySide6.QtCore import QSignalBlocker, Qt, Signal
from PySide6.QtGui import QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import (QCheckBox, QComboBox, QGridLayout, QHBoxLayout, QLabel, QMenu,
                               QToolButton, QVBoxLayout, QWidget)

NO_REFERENCE_HINT = ("No reference ligand yet — prepare a receptor with a co-crystal ligand, or "
                     "right-click any ligand in Project Explorer ▸ Use as pose reference.")
LOCK_TIP = ("Keep the chosen reference visible and hide the other references while browsing poses. "
            "Reference coordinates are never changed.")
CAMERA_TIP = "Keep the current camera when the pose changes. Forced on while a reference is locked."


class PoseReviewBar(QWidget):
    settings_changed = Signal()
    step_requested = Signal(int)
    fit_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.previous = QToolButton()
        self.previous.setText("◀")
        self.next = QToolButton()
        self.next.setText("▶")
        self.previous.setToolTip("Previous pose in the current table order (Alt+Left)")
        self.next.setToolTip("Next pose in the current table order (Alt+Right)")
        self.position = QLabel("Pose 0 / 0")
        self.position.setTextFormat(Qt.TextFormat.PlainText)
        self.position.setToolTip("Position within the poses currently listed by the table's sort and filter.")
        self.reference = QComboBox()
        self.reference.setMinimumContentsLength(10)
        self.reference.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.reference.setToolTip("Choose the reference ligand to keep visible while browsing.")
        self.reference_label = QLabel("Reference compound:")
        self.reference_label.setBuddy(self.reference)
        self.lock_reference = QCheckBox("Lock reference")
        self.lock_reference.setObjectName("lockReference")
        self.lock_reference.setToolTip(LOCK_TIP)

        self.display = QToolButton()
        self.display.setText("Display")
        self.display.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.display.setToolTip("Viewer options that apply while you browse poses.")
        menu = QMenu(self.display)
        self.keep_camera = QAction("Keep camera", menu, checkable=True, checked=True)
        self.keep_camera.setToolTip(CAMERA_TIP)
        self.hide_protein = QAction("Hide protein", menu, checkable=True)
        self.hide_protein.setToolTip("Hide protein objects without deleting them; uncheck to "
                                     "restore their visibility.")
        menu.addAction(self.keep_camera)
        menu.addAction(self.hide_protein)
        self.display.setMenu(menu)

        self.fit = QToolButton()
        self.fit.setText("Fit")
        self.fit.setToolTip("Frame the locked reference and current pose together. Never realigns atoms.")
        self.hint = QLabel(NO_REFERENCE_HINT)
        self.hint.setWordWrap(True)
        self.hint.setTextFormat(Qt.TextFormat.PlainText)
        self.hint.setEnabled(False)

        row = QHBoxLayout()
        row.setSpacing(4)
        for widget in (self.previous, self.next, self.position):
            row.addWidget(widget)
        row.addStretch(1)
        for widget in (self.reference_label, self.reference, self.lock_reference,
                       self.display, self.fit):
            row.addWidget(widget)
        # A grid, so the hidden hint's row collapses to nothing instead of leaving a gap.
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setVerticalSpacing(2)
        grid.addLayout(row, 0, 0)
        grid.addWidget(self.hint, 1, 0)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(grid)

        self.previous.clicked.connect(lambda: self.step_requested.emit(-1))
        self.next.clicked.connect(lambda: self.step_requested.emit(1))
        self.fit.clicked.connect(self.fit_requested)
        self.lock_reference.toggled.connect(self._lock_changed)
        self.reference.currentIndexChanged.connect(lambda _i: self.settings_changed.emit())
        self.hide_protein.toggled.connect(lambda _on: self.settings_changed.emit())
        for keys, button in (("Alt+Left", self.previous), ("Alt+Right", self.next)):
            shortcut = QShortcut(QKeySequence(keys), self)
            shortcut.setContext(Qt.ShortcutContext.WindowShortcut)
            shortcut.activated.connect(lambda b=button: self._step_shortcut(b))
        self.set_references([])
        self.set_position(-1, 0)

    def _step_shortcut(self, button) -> None:
        """Only step while the bar is on screen, so the keys stay inert in other workspaces."""
        if self.isVisible() and button.isEnabled():
            button.click()

    def _lock_changed(self, on):
        self.lock_reference.setText("Reference locked" if on else "Lock reference")
        if on:
            self.keep_camera.setChecked(True)
        self.keep_camera.setEnabled(not on)
        self.keep_camera.setToolTip(CAMERA_TIP if not on else
                                    "Enforced while the reference is locked, so pose changes cannot move the view.")
        self.settings_changed.emit()

    def set_references(self, names):
        previous = self.reference.currentText()
        blocker = QSignalBlocker(self.reference)
        self.reference.clear()
        self.reference.addItems(names)
        if previous in names:
            self.reference.setCurrentText(previous)
        del blocker
        self.reference.setEnabled(bool(names))
        self.lock_reference.setEnabled(bool(names))
        self.reference_label.setEnabled(bool(names))
        if not names:
            self.lock_reference.setChecked(False)
        else:
            self.lock_reference.setText("Reference locked"
                                        if self.lock_reference.isChecked() else "Lock reference")
        self.hint.setVisible(not names)
        self.settings_changed.emit()

    def select_reference(self, name: str) -> bool:
        """Show ``name`` as the locked reference. False when it is not an offered choice."""
        if self.reference.findText(name) < 0:
            return False
        self.reference.setCurrentText(name)
        self.lock_reference.setChecked(True)
        return True

    def set_position(self, index, total, detail=""):
        shown = f"{index+1 if index >= 0 else 0} / {total}"
        self.position.setText(f"Pose {shown}" + (f" — {detail}" if detail and index >= 0 else ""))
        self.previous.setEnabled(index > 0)
        self.next.setEnabled(0 <= index < total-1)
        self.fit.setEnabled(total > 0)
