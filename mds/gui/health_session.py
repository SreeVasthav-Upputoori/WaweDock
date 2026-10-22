"""Protein health workflow state: job lifecycle stays separate from widgets."""
from pathlib import Path
from PySide6.QtCore import QObject, Signal
from mds.jobs.latest import LatestJobQueue
from mds.jobs.protein_health import ProteinHealthJob, ProteinRepairPlanJob, ProteinRepairJob


class ProteinHealthSession(QObject):
    invalidated = Signal()
    report_ready = Signal(object)
    busy_changed = Signal(bool)
    progress_changed = Signal(int, str)
    failed = Signal(str)
    cancelled = Signal()
    plan_ready = Signal(object)
    candidate_ready = Signal(object)

    def __init__(self, runner, source, parent=None):
        super().__init__(parent)
        self.source = source
        self.queue = LatestJobQueue(runner,self)
        self.busy = False
        runner.job_progress.connect(self._progress)

    def invalidate(self, *_):
        self.queue.clear()
        self._busy(False)
        self.invalidated.emit()

    def _busy(self, value):
        self.busy = value
        self.busy_changed.emit(value)

    def check(self):
        if self.busy or self.source() is None: return
        source = Path(self.source())
        self.invalidated.emit()
        self._busy(True)
        try:
            self.queue.submit(ProteinHealthJob(source),self._done,self._failed)
        except Exception as exc:
            self._failed(str(exc))

    def cancel(self):
        self.queue.clear()
        self._busy(False)
        self.cancelled.emit()

    def _progress(self, job_id, percent, message):
        active = self.queue.active
        if self.busy and active and active[0]==job_id and active[1].revision==self.queue.revision:
            self.progress_changed.emit(percent,message)

    def _done(self, report):
        self._busy(False)
        if self.source() is not None and Path(report.source)==Path(self.source()):
            self.report_ready.emit(report)

    def _failed(self, message):
        self._busy(False)
        self.failed.emit(str(message))

    def propose(self, selected_residues=None, ph=7.4):
        if self.busy or self.source() is None: return
        self.invalidated.emit()
        self._busy(True)
        try: self.queue.submit(ProteinRepairPlanJob(self.source(),selected_residues,ph),self._plan_done,self._failed)
        except Exception as exc: self._failed(str(exc))

    def _plan_done(self,plan):
        self._busy(False)
        if self.source() is not None and Path(plan.source)==Path(self.source()): self.plan_ready.emit(plan)

    def build(self,plan,selected,output):
        if self.busy or self.source() is None: return
        if Path(plan.source)!=Path(self.source()):
            self.failed.emit('Repair proposal is stale; generate a new preview.')
            return
        self._busy(True)
        try: self.queue.submit(ProteinRepairJob(plan,selected,output),self._candidate_done,self._failed)
        except Exception as exc: self._failed(str(exc))

    def _candidate_done(self,result):
        self._busy(False)
        if self.source() is not None and Path(result.plan.source)==Path(self.source()): self.candidate_ready.emit(result)
