"""Execution log panel fed by Python ``logging``.

``QtLogHandler`` is a ``logging.Handler`` that forwards each record through a
Qt signal.  Records may be produced on worker threads; because the emitter
``QObject`` lives in the GUI thread, Qt queues the signal and the panel is
always updated on the GUI thread.
"""

from __future__ import annotations

import html
import logging
from collections import deque
from pathlib import Path

from PySide6.QtCore import QObject, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from mds.core.logging_setup import LOG_FORMAT

_LEVELS = {"Debug": logging.DEBUG, "Info": logging.INFO, "Warning": logging.WARNING, "Error": logging.ERROR}
_COLOURS = {logging.WARNING: "#d08a00", logging.ERROR: "#e05050", logging.CRITICAL: "#ff3030"}
_MAX_RECORDS = 5000


class _Emitter(QObject):
    record = Signal(int, str)  # level, formatted text


class QtLogHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.setFormatter(logging.Formatter(LOG_FORMAT, datefmt="%H:%M:%S"))
        self.emitter = _Emitter()

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.emitter.record.emit(record.levelno, self.format(record))
        except RuntimeError:
            pass  # emitter already destroyed during shutdown


class LogPanel(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._records: deque[tuple[int, str]] = deque(maxlen=_MAX_RECORDS)

        self._level = QComboBox()
        self._level.addItems(list(_LEVELS))
        self._level.setCurrentText("Info")
        self._level.currentTextChanged.connect(self._refilter)

        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear)
        save_btn = QPushButton("Save…")
        save_btn.clicked.connect(self._save)

        self._view = QPlainTextEdit()
        self._view.setReadOnly(True)
        self._view.setMaximumBlockCount(_MAX_RECORDS)
        self._view.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)

        bar = QHBoxLayout()
        bar.addWidget(QLabel("Level:"))
        bar.addWidget(self._level)
        bar.addStretch()
        bar.addWidget(clear_btn)
        bar.addWidget(save_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(bar)
        layout.addWidget(self._view)

        self.handler = QtLogHandler()
        self.handler.emitter.record.connect(self.append_record)

    def _min_level(self) -> int:
        return _LEVELS[self._level.currentText()]

    def append_record(self, level: int, text: str) -> None:
        self._records.append((level, text))
        if level >= self._min_level():
            self._write(level, text)

    def _write(self, level: int, text: str) -> None:
        escaped = html.escape(text).replace("\n", "<br>")
        colour = _COLOURS.get(level)
        self._view.appendHtml(
            f'<span style="color:{colour}">{escaped}</span>' if colour else escaped
        )

    def _refilter(self) -> None:
        self._view.clear()
        threshold = self._min_level()
        for level, text in self._records:
            if level >= threshold:
                self._write(level, text)

    def clear(self) -> None:
        self._records.clear()
        self._view.clear()

    def _save(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save Log", "mds-log.txt", "Text files (*.txt)")
        if path:
            Path(path).write_text("\n".join(t for _, t in self._records), encoding="utf-8")
