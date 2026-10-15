"""Persistence of window geometry and dock arrangements.

``QMainWindow.saveState()`` serialises the position, size, visibility and
tabbing of every dock and toolbar (identified by ``objectName``) into a
``QByteArray``.  We store it in ``QSettings``:

* ``session/*``          — restored automatically at the next start
* ``layouts/<name>/*``   — user-named layouts (View ▸ Layouts)

``LAYOUT_VERSION`` must be incremented whenever docks are added/renamed;
``restoreState`` then rejects stale data instead of producing a broken layout.
"""

from __future__ import annotations

import logging

from PySide6.QtCore import QByteArray, QSettings
from PySide6.QtWidgets import QMainWindow

log = logging.getLogger("mds.gui.layouts")

LAYOUT_VERSION = 3


class LayoutManager:
    def __init__(self, window: QMainWindow, settings: QSettings) -> None:
        self._window = window
        self._settings = settings
        self._default_state: QByteArray | None = None

    def capture_default(self) -> None:
        """Remember the as-built layout so *Reset Layout* can return to it."""
        self._default_state = self._window.saveState(LAYOUT_VERSION)

    def reset(self) -> None:
        if self._default_state is not None:
            self._window.restoreState(self._default_state, LAYOUT_VERSION)
            log.info("Layout reset to default")

    # ---------------------------------------------------------- session
    def save_session(self) -> None:
        self._settings.setValue("session/geometry", self._window.saveGeometry())
        self._settings.setValue("session/state", self._window.saveState(LAYOUT_VERSION))

    def restore_session(self) -> bool:
        geometry = self._settings.value("session/geometry")
        state = self._settings.value("session/state")
        if geometry is not None:
            self._window.restoreGeometry(geometry)
        if state is None:
            return False
        ok = self._window.restoreState(state, LAYOUT_VERSION)
        if not ok:
            log.warning("Stored session layout is incompatible; using default layout")
        return ok

    # ---------------------------------------------------- named layouts
    def names(self) -> list[str]:
        self._settings.beginGroup("layouts")
        try:
            return sorted(self._settings.childGroups())
        finally:
            self._settings.endGroup()

    def save_named(self, name: str) -> None:
        self._settings.setValue(f"layouts/{name}/state", self._window.saveState(LAYOUT_VERSION))
        log.info("Layout saved: %s", name)

    def restore_named(self, name: str) -> bool:
        state = self._settings.value(f"layouts/{name}/state")
        ok = state is not None and self._window.restoreState(state, LAYOUT_VERSION)
        if ok:
            log.info("Layout restored: %s", name)
        else:
            log.warning("Layout '%s' could not be restored", name)
        return ok

    def delete_named(self, name: str) -> None:
        self._settings.remove(f"layouts/{name}")
        log.info("Layout deleted: %s", name)
