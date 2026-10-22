"""Compact controls for fixed-reference, pose-by-pose visual comparison."""
from PySide6.QtCore import QSignalBlocker, Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget


class PoseReviewBar(QWidget):
    settings_changed = Signal()
    step_requested = Signal(int)
    fit_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.previous = QPushButton("Previous")
        self.next = QPushButton("Next")
        self.position = QLabel("0 / 0")
        self.reference = QComboBox()
        self.reference.setMinimumContentsLength(8)
        self.reference.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.reference.setToolTip("Choose the extracted crystal ligand to keep visible while browsing.")
        self.lock_reference = QCheckBox("Lock reference")
        self.keep_camera = QCheckBox("Keep camera")
        self.keep_camera.setChecked(True)
        self.hide_protein = QCheckBox("Hide protein")
        self.hide_protein.setToolTip("Hide protein objects without deleting them; uncheck to restore their visibility.")
        self.fit = QPushButton("Fit view")
        self.fit.setToolTip("Frame the locked reference and current pose together.")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        for widgets in ((self.previous,self.position,self.next),
                        (self.reference,self.lock_reference),
                        (self.keep_camera,self.hide_protein,self.fit)):
            row = QHBoxLayout()
            for widget in widgets:
                row.addWidget(widget)
            layout.addLayout(row)
        self.previous.clicked.connect(lambda:self.step_requested.emit(-1))
        self.next.clicked.connect(lambda:self.step_requested.emit(1))
        self.fit.clicked.connect(self.fit_requested)
        self.lock_reference.toggled.connect(self._lock_changed)
        self.reference.currentIndexChanged.connect(lambda _i:self.settings_changed.emit())
        self.hide_protein.toggled.connect(lambda _on:self.settings_changed.emit())
        self.set_references([])
        self.set_position(-1,0)

    def _lock_changed(self, on):
        if on:
            self.keep_camera.setChecked(True)
        self.keep_camera.setEnabled(not on)
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
        if not names:
            self.lock_reference.setChecked(False)
        self.settings_changed.emit()

    def set_position(self, index, total):
        self.position.setText(f"{index+1 if index>=0 else 0} / {total}")
        self.previous.setEnabled(index>0)
        self.next.setEnabled(0<=index<total-1)
        self.fit.setEnabled(total>0)
