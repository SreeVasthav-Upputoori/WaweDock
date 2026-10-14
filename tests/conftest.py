"""Shared test fixtures. Tests run headless via Qt's 'offscreen' platform."""

from __future__ import annotations

import os
import time
from collections.abc import Callable

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PySide6.QtCore import QCoreApplication, QSettings  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="session")
def qapp() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture
def settings(tmp_path) -> QSettings:
    """Isolated INI-file settings so tests never touch the user's registry."""
    return QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)


def wait_until(predicate: Callable[[], bool], timeout: float = 5.0) -> bool:
    """Process Qt events until ``predicate`` is true or the timeout expires."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        QCoreApplication.processEvents()
        if predicate():
            return True
        time.sleep(0.01)
    return False
