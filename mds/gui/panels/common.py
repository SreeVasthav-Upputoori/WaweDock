"""Small shared widgets for the analysis panels."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor, QFontMetrics
from PySide6.QtWidgets import QHBoxLayout, QLabel, QMessageBox, QSizePolicy, QWidget


class ElidedLabel(QLabel):
    """A single-line label that shortens its text to fit, keeping the whole of it as a tooltip.

    Wrapping is what makes a panel grow taller as its message gets longer, which is exactly
    backwards for a status or summary line sitting above a table.
    """

    def __init__(self, text: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._full = ""
        self.setTextFormat(Qt.TextFormat.PlainText)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.setText(text)

    def setText(self, text: str) -> None:  # noqa: N802 — overrides QLabel.setText
        self._full = text or ""
        self.setToolTip(self._full)
        self._elide()

    def text(self) -> str:
        """The text as set, not as shortened to fit.

        Elision is a rendering detail: a caller that sets a message and reads it back must get
        the message, not whatever happened to survive the current widget width.
        """
        return self._full

    def displayedText(self) -> str:  # noqa: N802 — Qt spelling
        """What is actually painted, shortened to the current width."""
        return super().text()

    def _elide(self) -> None:
        metrics = QFontMetrics(self.font())
        super().setText(metrics.elidedText(self._full, Qt.TextElideMode.ElideRight,
                                           max(self.width(), 40)))

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._elide()


class CaveatBar(QWidget):
    """A caveat kept to one line, with the full wording a hover or a click away.

    These notes are what stops a predicted number being read as a measured one, so they stay
    on screen rather than moving into a menu. Stacked as wrapped paragraphs under every tab,
    though, they cost more height than the data they qualify: three of them plus the score
    disclaimer took about a fifth of the results panel. Elided to a single line the warning is
    still always visible, and the complete text is one hover or click away.
    """

    def __init__(self, text: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._text = ""
        self._title = "About these numbers"
        self.icon = QLabel("ⓘ")          # circled i
        self.icon.setToolTip("Click for the full note")
        self.label = ElidedLabel()
        for widget in (self, self.icon, self.label):
            widget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        layout = QHBoxLayout(self)
        layout.setContentsMargins(2, 0, 2, 0)
        layout.setSpacing(4)
        layout.addWidget(self.icon)
        layout.addWidget(self.label, 1)
        self.setText(text)

    def setText(self, text: str) -> None:  # noqa: N802 — matches QLabel's spelling
        self._text = text
        self.setToolTip(text)
        self.label.setText(text)

    def text(self) -> str:
        return self._text

    def mousePressEvent(self, event) -> None:  # noqa: N802
        super().mousePressEvent(event)
        self.show_full_text()

    def show_full_text(self) -> None:
        box = QMessageBox(QMessageBox.Icon.Information, self._title, self._text,
                          QMessageBox.StandardButton.Close, self)
        box.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        box.exec()
