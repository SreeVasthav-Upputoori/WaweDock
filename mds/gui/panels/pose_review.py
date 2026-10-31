"""Compact controls for fixed-reference, pose-by-pose visual comparison."""
from PySide6.QtCore import QSignalBlocker, Qt, Signal
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (QCheckBox, QComboBox, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                               QPushButton, QVBoxLayout, QWidget)

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
        self.previous = QPushButton("◀ Previous")
        self.next = QPushButton("Next ▶")
        self.previous.setToolTip("Previous pose in the current table order (Alt+Left)")
        self.next.setToolTip("Next pose in the current table order (Alt+Right)")
        self.position = QLabel("0 / 0")
        self.position.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.position.setTextFormat(Qt.TextFormat.PlainText)
        self.position.setToolTip("Position within the poses currently listed by the table's sort and filter.")
        self.reference = QComboBox()
        self.reference.setMinimumContentsLength(8)
        self.reference.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.reference.setToolTip("Choose the reference ligand to keep visible while browsing.")
        self.reference_label = QLabel("Reference:")
        self.reference_label.setBuddy(self.reference)
        self.lock_reference = QCheckBox("Lock reference")
        self.lock_reference.setToolTip(LOCK_TIP)
        self.keep_camera = QCheckBox("Keep camera")
        self.keep_camera.setChecked(True)
        self.keep_camera.setToolTip(CAMERA_TIP)
        self.hide_protein = QCheckBox("Hide protein")
        self.hide_protein.setToolTip("Hide protein objects without deleting them; uncheck to restore their visibility.")
        self.fit = QPushButton("Fit view")
        self.fit.setToolTip("Frame the locked reference and current pose together. Never realigns atoms.")
        self.hint = QLabel(NO_REFERENCE_HINT)
        self.hint.setWordWrap(True)
        self.hint.setTextFormat(Qt.TextFormat.PlainText)
        self.hint.setEnabled(False)

        box = QGroupBox("Pose review")
        grid = QGridLayout(box)
        grid.setHorizontalSpacing(6)
        grid.setVerticalSpacing(4)
        grid.addWidget(self.previous, 0, 0)
        grid.addWidget(self.position, 0, 1)
        grid.addWidget(self.next, 0, 2)
        grid.setColumnStretch(1, 1)
        reference_row = QHBoxLayout()
        reference_row.addWidget(self.reference_label)
        reference_row.addWidget(self.reference, 1)
        reference_row.addWidget(self.lock_reference)
        grid.addLayout(reference_row, 1, 0, 1, 3)
        toggles = QHBoxLayout()
        toggles.addWidget(self.keep_camera)
        toggles.addWidget(self.hide_protein)
        toggles.addStretch(1)
        toggles.addWidget(self.fit)
        grid.addLayout(toggles, 2, 0, 1, 3)
        grid.addWidget(self.hint, 3, 0, 1, 3)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(box)

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
