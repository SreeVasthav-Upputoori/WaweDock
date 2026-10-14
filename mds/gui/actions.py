"""Central registry of ``QAction`` objects.

A single ``QAction`` can appear in a menu, a toolbar and a context menu at
the same time; enabling/disabling it updates every place at once.  Keeping
them in one registry, addressed by a stable string key, lets later phases
look actions up (``registry["docking.run"]``) instead of passing widgets
around.

Actions for features that do not exist yet are created *disabled* with a
"Planned — Phase N" tooltip, so the menu structure is visible without
pretending the feature works.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtGui import QAction, QIcon, QKeySequence
from PySide6.QtWidgets import QWidget


class ActionRegistry:
    def __init__(self, parent: QWidget) -> None:
        self._parent = parent
        self._actions: dict[str, QAction] = {}

    def add(
        self,
        key: str,
        text: str,
        *,
        slot: Callable[..., object] | None = None,
        shortcut: QKeySequence | QKeySequence.StandardKey | str | None = None,
        icon: QIcon | None = None,
        tip: str = "",
        planned_phase: int | None = None,
        checkable: bool = False,
    ) -> QAction:
        if key in self._actions:
            raise KeyError(f"Duplicate action key: {key}")
        action = QAction(text, self._parent)
        action.setObjectName(key)
        if icon is not None:
            action.setIcon(icon)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
        action.setCheckable(checkable)

        if planned_phase is not None:
            action.setEnabled(False)
            tip = f"{tip + ' — ' if tip else ''}Planned — Phase {planned_phase}"
        elif slot is not None:
            signal = action.toggled if checkable else action.triggered
            signal.connect(slot)

        if tip:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        self._actions[key] = action
        return action

    def __getitem__(self, key: str) -> QAction:
        return self._actions[key]

    def __contains__(self, key: str) -> bool:
        return key in self._actions

    def keys(self) -> list[str]:
        return list(self._actions)
