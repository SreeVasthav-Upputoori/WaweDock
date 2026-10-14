from __future__ import annotations

import threading

from conftest import wait_until

from mds.jobs import Job, JobContext, JobRunner, JobState
from mds.jobs.diagnostics import JobSystemSelfTest


class _Add(Job):
    name = "add"

    def run(self, ctx: JobContext) -> int:
        ctx.progress(50, "half")
        return 2 + 3


class _Boom(Job):
    name = "boom"

    def run(self, ctx: JobContext) -> None:
        raise ValueError("deliberate failure")


def test_job_finishes_with_result_on_gui_thread(qapp):
    runner = JobRunner()
    results, progress, threads = [], [], []
    runner.job_finished.connect(lambda i, r: (results.append(r), threads.append(threading.current_thread())))
    runner.job_progress.connect(lambda i, p, m: progress.append((p, m)))

    job_id = runner.submit(_Add())
    assert wait_until(lambda: runner.state(job_id) is JobState.FINISHED)
    assert results == [5]
    assert (50, "half") in progress
    assert threads == [threading.main_thread()]  # slots run on the GUI thread
    assert runner.active_count() == 0


def test_job_failure_is_reported_with_traceback(qapp):
    runner = JobRunner()
    errors = []
    runner.job_failed.connect(lambda i, tb: errors.append(tb))
    job_id = runner.submit(_Boom())
    assert wait_until(lambda: runner.state(job_id) is JobState.FAILED)
    assert "deliberate failure" in errors[0]


def test_job_cancellation(qapp):
    runner = JobRunner()
    job_id = runner.submit(JobSystemSelfTest(steps=200, interval_s=0.02))
    assert wait_until(lambda: runner.state(job_id) is JobState.RUNNING)
    runner.cancel(job_id)
    assert wait_until(lambda: runner.state(job_id) is JobState.CANCELLED)
    assert runner.wait_for_done(2000)
