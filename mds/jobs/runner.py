"""Thread-pool job runner with progress reporting and cooperative cancellation.

How it fits together
--------------------
* A :class:`Job` is a plain Python object whose :meth:`Job.run` does the work.
  It never touches widgets.
* :class:`JobRunner` executes jobs on a ``QThreadPool``.  Worker threads emit
  signals on a private ``_Signals`` object that lives in the GUI thread, so Qt
  delivers them as *queued* calls — every ``JobRunner`` signal is therefore
  received on the GUI thread and slots may update widgets safely.
* Cancellation is cooperative: :meth:`JobRunner.cancel` sets an event and the
  job calls :meth:`JobContext.check_cancelled` at safe points.  Jobs that wrap
  external programs (e.g. Vina) poll the event and terminate the process.
"""

from __future__ import annotations

import itertools
import logging
import threading
import traceback
from abc import ABC, abstractmethod
from enum import Enum
from typing import Any

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal

log = logging.getLogger("mds.jobs")


class JobCancelled(Exception):
    """Raised inside a job when the user requested cancellation."""


class JobState(Enum):
    PENDING = "Pending"
    RUNNING = "Running"
    FINISHED = "Finished"
    FAILED = "Failed"
    CANCELLED = "Cancelled"

    @property
    def is_terminal(self) -> bool:
        return self in (JobState.FINISHED, JobState.FAILED, JobState.CANCELLED)


class _Signals(QObject):
    started = Signal(int)
    progress = Signal(int, int, str)
    finished = Signal(int, object)
    failed = Signal(int, str)
    cancelled = Signal(int)


class JobContext:
    """Handed to :meth:`Job.run`; the job's only channel back to the application."""

    def __init__(self, job_id: int, signals: _Signals, cancel_event: threading.Event) -> None:
        self.job_id = job_id
        self._signals = signals
        self._cancel_event = cancel_event

    def progress(self, percent: float, message: str = "") -> None:
        self._signals.progress.emit(self.job_id, max(0, min(100, int(percent))), message)

    @property
    def cancelled(self) -> bool:
        return self._cancel_event.is_set()

    def check_cancelled(self) -> None:
        if self._cancel_event.is_set():
            raise JobCancelled()


class Job(ABC):
    """Base class for all background work."""

    name: str = "Job"

    @abstractmethod
    def run(self, ctx: JobContext) -> Any:
        """Do the work in a worker thread and return a result object."""


class _JobRunnable(QRunnable):
    def __init__(
        self, job_id: int, job: Job, signals: _Signals, cancel_event: threading.Event
    ) -> None:
        super().__init__()
        self._job_id = job_id
        self._job = job
        self._signals = signals
        self._cancel_event = cancel_event

    def run(self) -> None:
        if self._cancel_event.is_set():
            self._signals.cancelled.emit(self._job_id)
            return
        self._signals.started.emit(self._job_id)
        ctx = JobContext(self._job_id, self._signals, self._cancel_event)
        try:
            result = self._job.run(ctx)
        except JobCancelled:
            self._signals.cancelled.emit(self._job_id)
        except Exception:  # noqa: BLE001 — any job failure must reach the GUI
            self._signals.failed.emit(self._job_id, traceback.format_exc())
        else:
            if self._cancel_event.is_set():
                self._signals.cancelled.emit(self._job_id)
            else:
                self._signals.finished.emit(self._job_id, result)


class JobRunner(QObject):
    """Submits jobs to a thread pool and publishes their lifecycle as signals."""

    job_submitted = Signal(int, str)  # id, name
    job_started = Signal(int)
    job_progress = Signal(int, int, str)  # id, percent, message
    job_finished = Signal(int, object)  # id, result
    job_failed = Signal(int, str)  # id, traceback text
    job_cancelled = Signal(int)
    active_count_changed = Signal(int)

    def __init__(self, max_threads: int | None = None, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._pool = QThreadPool(self)
        if max_threads:
            self._pool.setMaxThreadCount(max_threads)
        self._ids = itertools.count(1)
        self._names: dict[int, str] = {}
        self._states: dict[int, JobState] = {}
        self._cancel_events: dict[int, threading.Event] = {}

        self._signals = _Signals(self)
        self._signals.started.connect(self._on_started)
        self._signals.progress.connect(self.job_progress)
        self._signals.finished.connect(self._on_finished)
        self._signals.failed.connect(self._on_failed)
        self._signals.cancelled.connect(self._on_cancelled)

    # ----------------------------------------------------------------- API
    def submit(self, job: Job) -> int:
        job_id = next(self._ids)
        event = threading.Event()
        self._names[job_id] = job.name
        self._states[job_id] = JobState.PENDING
        self._cancel_events[job_id] = event
        self._pool.start(_JobRunnable(job_id, job, self._signals, event))
        log.info("Job %d submitted: %s", job_id, job.name)
        self.job_submitted.emit(job_id, job.name)
        self.active_count_changed.emit(self.active_count())
        return job_id

    def cancel(self, job_id: int) -> None:
        event = self._cancel_events.get(job_id)
        if event is not None and not self._states[job_id].is_terminal:
            log.info("Cancellation requested for job %d (%s)", job_id, self._names[job_id])
            event.set()

    def cancel_all(self) -> None:
        for job_id in list(self._cancel_events):
            self.cancel(job_id)

    def state(self, job_id: int) -> JobState:
        return self._states[job_id]

    def name(self, job_id: int) -> str:
        return self._names[job_id]

    def active_count(self) -> int:
        return sum(1 for s in self._states.values() if not s.is_terminal)

    def wait_for_done(self, msecs: int = -1) -> bool:
        return self._pool.waitForDone(msecs)

    # ------------------------------------------------------------ internals
    def _set_terminal(self, job_id: int, state: JobState) -> None:
        self._states[job_id] = state
        self._cancel_events.pop(job_id, None)
        self.active_count_changed.emit(self.active_count())

    def _on_started(self, job_id: int) -> None:
        self._states[job_id] = JobState.RUNNING
        log.debug("Job %d started", job_id)
        self.job_started.emit(job_id)

    def _on_finished(self, job_id: int, result: object) -> None:
        self._set_terminal(job_id, JobState.FINISHED)
        log.info("Job %d finished: %s", job_id, self._names[job_id])
        self.job_finished.emit(job_id, result)

    def _on_failed(self, job_id: int, tb: str) -> None:
        self._set_terminal(job_id, JobState.FAILED)
        log.error("Job %d failed: %s\n%s", job_id, self._names[job_id], tb)
        self.job_failed.emit(job_id, tb)

    def _on_cancelled(self, job_id: int) -> None:
        self._set_terminal(job_id, JobState.CANCELLED)
        log.warning("Job %d cancelled: %s", job_id, self._names[job_id])
        self.job_cancelled.emit(job_id)
