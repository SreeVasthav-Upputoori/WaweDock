"""Light and dark application themes.

Uses Qt's *Fusion* style with an explicit ``QPalette`` so both themes look
identical on every Windows version, independent of the OS colour scheme.
"""

from __future__ import annotations

from enum import Enum

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication


class Theme(str, Enum):
    LIGHT = "light"
    DARK = "dark"


_COMMON_QSS = """
QDockWidget::title { padding: 4px 6px; font-weight: 600; }
QToolBar { spacing: 2px; }
QTreeView, QTableView, QPlainTextEdit { border: none; }
QStatusBar QLabel { padding: 0 6px; }
"""


def _dark_palette() -> QPalette:
    window = QColor(45, 47, 51)
    base = QColor(30, 31, 34)
    text = QColor(220, 222, 226)
    disabled = QColor(120, 122, 128)
    accent = QColor(64, 140, 230)

    p = QPalette()
    p.setColor(QPalette.ColorRole.Window, window)
    p.setColor(QPalette.ColorRole.WindowText, text)
    p.setColor(QPalette.ColorRole.Base, base)
    p.setColor(QPalette.ColorRole.AlternateBase, QColor(38, 40, 44))
    p.setColor(QPalette.ColorRole.ToolTipBase, QColor(60, 62, 66))
    p.setColor(QPalette.ColorRole.ToolTipText, text)
    p.setColor(QPalette.ColorRole.PlaceholderText, disabled)
    p.setColor(QPalette.ColorRole.Text, text)
    p.setColor(QPalette.ColorRole.Button, window)
    p.setColor(QPalette.ColorRole.ButtonText, text)
    p.setColor(QPalette.ColorRole.BrightText, QColor(255, 90, 90))
    p.setColor(QPalette.ColorRole.Link, accent)
    p.setColor(QPalette.ColorRole.Highlight, accent)
    p.setColor(QPalette.ColorRole.HighlightedText, QColor(255, 255, 255))
    for role in (
        QPalette.ColorRole.WindowText,
        QPalette.ColorRole.Text,
        QPalette.ColorRole.ButtonText,
    ):
        p.setColor(QPalette.ColorGroup.Disabled, role, disabled)
    return p


def apply_theme(app: QApplication, theme: Theme) -> None:
    app.setStyle("Fusion")
    palette = _dark_palette() if theme is Theme.DARK else app.style().standardPalette()
    app.setPalette(palette)
    app.setStyleSheet(_COMMON_QSS)
