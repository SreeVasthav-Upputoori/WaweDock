"""Central workspace host for the 3D viewer.

Phase 2 will insert a ``ViewerBackend`` widget via :meth:`ViewportHost.set_viewer`.
Until then the host shows an explicit notice — it does not imitate a viewer.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFrame, QLabel, QVBoxLayout, QWidget


class ViewportHost(QFrame):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("viewportHost")
        self.setFrameShape(QFrame.Shape.NoFrame)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

        self._notice = QLabel(
            "<h2>3D Workspace</h2>"
            "<p>The molecular viewer is unavailable in this session.</p>"
            "<p>You can still prepare structures and run docking using the menus.</p>"
        )
        self._notice.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._notice.setEnabled(False)  # renders with the theme's muted text colour
        self._layout.addWidget(self._notice)
        self._viewer: QWidget | None = None

    def set_viewer(self, viewer: QWidget) -> None:
        """Replace the notice with a real viewer widget."""
        if self._viewer is not None:
            self._layout.removeWidget(self._viewer)
            self._viewer.deleteLater()
        self._notice.hide()
        self._viewer = viewer
        self._layout.addWidget(viewer)
