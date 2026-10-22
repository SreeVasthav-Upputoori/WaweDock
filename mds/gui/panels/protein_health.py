"""Collapsible health dock content with asynchronous, stale-safe diagnostics."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QCheckBox, QFileDialog, QLabel, QListWidget,
                               QListWidgetItem, QProgressBar, QPushButton, QTextBrowser,
                               QVBoxLayout, QWidget)


class ProteinHealthPanel(QWidget):
    issue_requested = Signal(object)
    highlights_changed = Signal(object)

    def __init__(self, controller, parent=None):
        super().__init__(parent)
        self.c = controller
        self.report = None
        self.session = controller.protein_health
        self._busy = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8,8,8,8)
        title = QLabel('Protein Health & Repair')
        title.setWordWrap(True)
        title.setStyleSheet('font-weight: 600;')
        layout.addWidget(title)
        self.summary = QLabel('Check the loaded structure. Screening does not change coordinates.')
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.check_btn = QPushButton('Check Protein')
        self.validate_btn = QPushButton('Validate Structure')
        self.validate_btn.setToolTip('Repeat the supported geometric screens on the current structure. Not a full stereochemical or docking-readiness validation.')
        for button in (self.check_btn,self.validate_btn):
            layout.addWidget(button)
            button.clicked.connect(self.check)
        self.progress = QProgressBar()
        self.progress.setRange(0,100)
        self.progress.hide()
        layout.addWidget(self.progress)
        self.cancel_btn = QPushButton('Cancel check')
        self.cancel_btn.clicked.connect(self.session.cancel)
        self.cancel_btn.hide()
        layout.addWidget(self.cancel_btn)
        self.highlights = QCheckBox('Show issue markers')
        self.highlights.setChecked(True)
        self.highlights.setToolTip('Small crosses: red errors, orange warnings, blue notes. Selection remains yellow. At most 250 markers are drawn; the issue list is authoritative.')
        self.highlights.toggled.connect(self._highlights)
        layout.addWidget(self.highlights)
        self.issues = QListWidget()
        self.issues.setWordWrap(True)
        self.issues.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.issues.currentRowChanged.connect(self._inspect)
        layout.addWidget(self.issues,2)
        self.details = QTextBrowser()
        self.details.setMinimumWidth(0)
        self.details.setMinimumHeight(100)
        self.details.setOpenExternalLinks(False)
        layout.addWidget(self.details,1)
        self.export_btn = QPushButton('Export validation report')
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self._export)
        layout.addWidget(self.export_btn)
        self.session.progress_changed.connect(self._progress)
        self.session.invalidated.connect(self.invalidate)
        self.session.busy_changed.connect(self._set_busy)
        self.session.report_ready.connect(self._done)
        self.session.failed.connect(self._failed)
        self.session.cancelled.connect(lambda: self.summary.setText('Check cancelled. No coordinates changed.'))
        self._availability()

    def _availability(self):
        ready = self.c.source_pdb is not None and not self._busy
        self.check_btn.setEnabled(ready)
        self.validate_btn.setEnabled(ready)

    def invalidate(self, *_):
        self.report = None
        self.issues.clear()
        self.details.clear()
        self.export_btn.setEnabled(False)
        self.summary.setText('Structure changed. Run Check Protein on the current coordinates.')
        self._set_busy(False)
        self._highlights()

    def _set_busy(self, value):
        self._busy = value
        if value:
            self.summary.setText("Checking protein...")
            self.progress.setValue(0)
        self.progress.setVisible(value)
        self.cancel_btn.setVisible(value)
        self._availability()

    def check(self):
        self.session.check()

    def _progress(self, percent, message):
        self.progress.setValue(percent)
        self.summary.setText(message)

    def _done(self, report):
        self._set_busy(False)
        if self.c.source_pdb is None or Path(report.source) != Path(self.c.source_pdb): return
        self.report = report
        self.summary.setText(report.summary)
        self.issues.setUpdatesEnabled(False)
        try:
            for issue in report.issues:
                item = QListWidgetItem(f'{issue.severity} · {issue.location}\n{issue.code.replace("_", " ")}')
                item.setToolTip(issue.message)
                self.issues.addItem(item)
        finally:
            self.issues.setUpdatesEnabled(True)
        self.details.setPlainText('SUPPORTED CHECKS\n' + '\n'.join(report.coverage) + '\n\nLIMITATIONS\n' + '\n'.join(report.limitations))
        self.export_btn.setEnabled(True)
        self._highlights()

    def _failed(self, message):
        self._set_busy(False)
        self.summary.setText('Protein check failed. No coordinates changed.')
        self.details.setPlainText(str(message))

    def _inspect(self, row):
        if self.report is None or row < 0 or row >= len(self.report.issues): return
        issue = self.report.issues[row]
        self.details.setPlainText(f'{issue.severity} · {issue.location}\n{issue.message}\n\nRecommended action\n{issue.recommendation}\n\nStatus: {issue.status}\n\nRead-only screening; no repair has been performed.')
        self.issue_requested.emit(issue)

    def inspect_selection(self, selections):
        """Find a diagnosed issue for an exact existing viewer atom selection."""
        if self.report is None:
            self.summary.setText('Run Check Protein to inspect issues in this selection.')
            return
        indices = {item['atom_index'] for item in selections if 'atom_index' in item}
        for row,issue in enumerate(self.report.issues):
            if any(a['atom_index'] in indices for a in issue.atoms):
                if self.issues.currentRow() == row:
                    self._inspect(row)
                else:
                    self.issues.setCurrentRow(row)
                return
        self.details.setPlainText('No reported issue overlaps this selection. This does not establish that the region is free of abnormalities; review screening limits.')

    def _highlights(self, *_):
        self.highlights_changed.emit(self.report if self.highlights.isChecked() else None)

    def _export(self):
        if self.report is None: return
        path,_ = QFileDialog.getSaveFileName(self,'Export validation report','protein-health.json','JSON report (*.json)')
        if path:
            try: self.report.export(Path(path))
            except (OSError,ValueError) as exc: self.details.setPlainText(f'Cannot export report: {exc}')
