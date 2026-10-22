"""Run at most one job at a time, retaining only the latest queued UI request."""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from PySide6.QtCore import QObject
from mds.jobs.runner import Job, JobRunner


@dataclass
class _Request:
    job: Job
    done: Callable
    failed: Callable
    revision: int


class LatestJobQueue(QObject):
    def __init__(self, runner: JobRunner, parent=None):
        super().__init__(parent)
        self.runner = runner
        self.revision = 0
        self.active = None
        self.pending = None
        runner.job_finished.connect(self._finished)
        runner.job_failed.connect(self._failed)
        runner.job_cancelled.connect(lambda job_id: self._failed(job_id, "Cancelled"))

    def submit(self, job: Job, done: Callable, failed: Callable):
        self.revision += 1
        self.pending = _Request(job, done, failed, self.revision)
        self._start_pending()

    def clear(self):
        self.revision += 1
        self.pending = None
        if self.active is not None:
            self.runner.cancel(self.active[0])

    def _start_pending(self):
        if self.active is not None or self.pending is None:
            return
        request = self.pending
        self.pending = None
        # job_submitted observers can enqueue another request during submit().
        self.active = (-1, request)
        try:
            self.active = (self.runner.submit(request.job), request)
        except Exception:
            self.active = None
            raise

    def _finished(self, job_id, result):
        self._complete(job_id, result, None)

    def _failed(self, job_id, message):
        self._complete(job_id, None, message)

    def _complete(self, job_id, result, error):
        if self.active is None or self.active[0] != job_id:
            return
        _, request = self.active
        self.active = None
        try:
            if request.revision == self.revision:
                if error is None:
                    request.done(result)
                else:
                    request.failed(error)
        finally:
            self._start_pending()
