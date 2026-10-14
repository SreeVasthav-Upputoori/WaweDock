"""Clearly-labelled panel for workspaces whose functionality arrives in a later phase."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget


class PlannedFeaturePanel(QWidget):
    def __init__(self, title: str, phase: int, description: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        label = QLabel(f"<b>{title}</b><br><br>Planned — Phase {phase}<br><br>{description}")
        label.setWordWrap(True)
        label.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        label.setEnabled(False)
        layout = QVBoxLayout(self)
        layout.addWidget(label)
