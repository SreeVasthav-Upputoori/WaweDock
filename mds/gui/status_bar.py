"""Status bar: transient messages, current selection and background-job progress."""

from __future__ import annotations

from PySide6.QtWidgets import QLabel, QProgressBar, QStatusBar, QWidget


class StudioStatusBar(QStatusBar):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._selection = QLabel("No selection")
        self._jobs = QLabel("Idle")
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setFixedWidth(160)
        self._progress.setTextVisible(True)
        self._progress.hide()

        self.addPermanentWidget(self._selection)
        self.addPermanentWidget(self._jobs)
        self.addPermanentWidget(self._progress)

    def set_selection(self, text: str) -> None:
        self._selection.setText(text or "No selection")

    def set_active_jobs(self, count: int) -> None:
        self._jobs.setText("Idle" if count == 0 else f"{count} job(s) running")
        self._progress.setVisible(count > 0)
        if count == 0:
            self._progress.setValue(0)

    def set_progress(self, percent: int) -> None:
        self._progress.setValue(percent)
