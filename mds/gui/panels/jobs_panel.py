"""Table of background jobs with live progress and cancellation."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QHeaderView,
    QProgressBar,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from mds.jobs import JobRunner, JobState

_COLUMNS = ["ID", "Job", "State", "Progress", "Message"]


class JobsPanel(QWidget):
    def __init__(self, runner: JobRunner, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._runner = runner
        self._rows: dict[int, int] = {}

        self._table = QTableWidget(0, len(_COLUMNS))
        self._table.setHorizontalHeaderLabels(_COLUMNS)
        self._table.verticalHeader().hide()
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setStretchLastSection(True)

        cancel_btn = QPushButton("Cancel Selected")
        cancel_btn.clicked.connect(self._cancel_selected)
        clear_btn = QPushButton("Remove Completed")
        clear_btn.clicked.connect(self._remove_completed)

        bar = QHBoxLayout()
        bar.addStretch()
        bar.addWidget(cancel_btn)
        bar.addWidget(clear_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(bar)
        layout.addWidget(self._table)

        runner.job_submitted.connect(self._on_submitted)
        runner.job_started.connect(lambda i: self._set_state(i, JobState.RUNNING))
        runner.job_progress.connect(self._on_progress)
        runner.job_finished.connect(lambda i, _r: self._set_state(i, JobState.FINISHED, 100))
        runner.job_failed.connect(lambda i, _t: self._set_state(i, JobState.FAILED))
        runner.job_cancelled.connect(lambda i: self._set_state(i, JobState.CANCELLED))

    def row_count(self) -> int:
        return self._table.rowCount()

    def _on_submitted(self, job_id: int, name: str) -> None:
        row = self._table.rowCount()
        self._table.insertRow(row)
        self._rows[job_id] = row
        self._table.setItem(row, 0, QTableWidgetItem(str(job_id)))
        self._table.setItem(row, 1, QTableWidgetItem(name))
        self._table.setItem(row, 2, QTableWidgetItem(JobState.PENDING.value))
        bar = QProgressBar()
        bar.setRange(0, 100)
        self._table.setCellWidget(row, 3, bar)
        self._table.setItem(row, 4, QTableWidgetItem(""))

    def _on_progress(self, job_id: int, percent: int, message: str) -> None:
        row = self._rows.get(job_id)
        if row is None:
            return
        self._table.cellWidget(row, 3).setValue(percent)
        self._table.item(row, 4).setText(message)

    def _set_state(self, job_id: int, state: JobState, percent: int | None = None) -> None:
        row = self._rows.get(job_id)
        if row is None:
            return
        self._table.item(row, 2).setText(state.value)
        if percent is not None:
            self._table.cellWidget(row, 3).setValue(percent)

    def _cancel_selected(self) -> None:
        for index in self._table.selectionModel().selectedRows():
            self._runner.cancel(int(self._table.item(index.row(), 0).text()))

    def _remove_completed(self) -> None:
        for job_id, row in sorted(self._rows.items(), key=lambda kv: kv[1], reverse=True):
            if self._runner.state(job_id).is_terminal:
                self._table.removeRow(row)
                del self._rows[job_id]
        # Re-index the remaining rows.
        self._rows = {int(self._table.item(r, 0).text()): r for r in range(self._table.rowCount())}
