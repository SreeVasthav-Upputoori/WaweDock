"""Diagnostic jobs used to verify the job system itself (not scientific work)."""

from __future__ import annotations

import time

from mds.jobs.runner import Job, JobContext


class JobSystemSelfTest(Job):
    """Sleeps in small steps, reporting progress and honouring cancellation.

    Used from *Tools ▸ Diagnostics* to confirm the GUI stays responsive while
    a background job runs.
    """

    name = "Job system self-test"

    def __init__(self, steps: int = 50, interval_s: float = 0.1) -> None:
        self.steps = steps
        self.interval_s = interval_s

    def run(self, ctx: JobContext) -> str:
        for i in range(1, self.steps + 1):
            ctx.check_cancelled()
            time.sleep(self.interval_s)
            ctx.progress(100 * i / self.steps, f"step {i}/{self.steps}")
        return f"Completed {self.steps} steps"
