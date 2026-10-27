"""Collapsible health dock content with asynchronous, stale-safe diagnostics."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, Signal
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QListWidget,
                               QListWidgetItem, QProgressBar, QPushButton, QScrollArea, QFrame, QSizePolicy, QTextBrowser,
                               QVBoxLayout, QWidget)


class ProteinHealthPanel(QWidget):
    issue_requested = Signal(object)
    highlights_changed = Signal(object)
    rebuilt_highlights_changed = Signal(bool)
    comparison_requested = Signal(str)

    def __init__(self, controller, parent=None):
        super().__init__(parent)
        self.c = controller
        self.report = None
        self.repair_plan = None
        self.repair_result = None
        self.accepted_result = None
        self._filling = False
        self._can_undo = False
        self.session = controller.protein_health
        self._busy = False
        outer=QVBoxLayout(self)
        outer.setContentsMargins(0,0,0,0)
        self.scroll=QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content=QWidget()
        self.scroll.setWidget(content)
        outer.addWidget(self.scroll)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(8,8,8,8)
        title = QLabel('Protein Health & Repair')
        title.setWordWrap(True)
        title.setStyleSheet('font-weight: 600;')
        layout.addWidget(title)
        self.summary = QLabel('Check the loaded structure. Screening is read-only.')
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.check_btn = QPushButton('Check Protein')
        self.repair_btn = QPushButton('Repair Protein')
        self.repair_btn.setToolTip('Preview selectable corrections before building a separate repair candidate.')
        self.repair_btn.clicked.connect(self.propose)
        self.validate_btn = QPushButton('Validate Structure')
        self.validate_btn.setToolTip('Repeat the supported geometric screens on the current structure. Not a full stereochemical or docking-readiness validation.')
        for button in (self.check_btn,self.repair_btn,self.validate_btn):
            layout.addWidget(button)
            if button is not self.repair_btn: button.clicked.connect(self.check)
        self.repair_group = QGroupBox('Proposed corrections')
        repair_layout = QVBoxLayout(self.repair_group)
        self.repair_scope = QComboBox()
        self.repair_scope.addItems(['All protein residues','Viewer selection'])
        self.repair_limits = QLabel('Atom completion and optional short-loop models (up to 8 residues). Loops require sequence evidence and remain unchecked until you approve them. Backbone torsion quality needs independent review.')
        self.repair_limits.setWordWrap(True)
        repair_layout.addWidget(self.repair_limits)
        repair_layout.addWidget(self.repair_scope)
        row = QHBoxLayout()
        row.addWidget(QLabel('Hydrogen pH'))
        self.repair_ph = QDoubleSpinBox()
        self.repair_ph.setRange(0,14)
        self.repair_ph.setValue(7.4)
        self.repair_ph.setDecimals(1)
        self.repair_ph.setSingleStep(0.1)
        row.addWidget(self.repair_ph)
        repair_layout.addLayout(row)
        self.proposed_changes = QListWidget()
        self.proposed_changes.setWordWrap(True)
        self.proposed_changes.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.proposed_changes.setMaximumHeight(150)
        repair_layout.addWidget(self.proposed_changes)
        self.build_btn = QPushButton('Build selected repairs')
        self.build_btn.clicked.connect(self._build)
        repair_layout.addWidget(self.build_btn)
        self.comparison = QComboBox()
        self.comparison.addItems(['Before: original','After: reconstructed'])
        self.comparison.currentIndexChanged.connect(lambda i: self.comparison_requested.emit('after' if i else 'before'))
        repair_layout.addWidget(self.comparison)
        self.candidate_summary=QLabel()
        self.candidate_summary.setWordWrap(True)
        self.candidate_summary.hide()
        repair_layout.addWidget(self.candidate_summary)
        self.accept_btn = QPushButton('Use repaired structure')
        self.accept_btn.clicked.connect(self._accept)
        repair_layout.addWidget(self.accept_btn)
        self.undo_btn = QPushButton('Undo structure edit')
        self.undo_btn.clicked.connect(self.c.undo_structure_edit)
        repair_layout.addWidget(self.undo_btn)
        self.structure_export_btn = QPushButton('Export repaired structure')
        self.structure_export_btn.clicked.connect(self._export_structure)
        repair_layout.addWidget(self.structure_export_btn)
        self.repair_scope.currentIndexChanged.connect(self._settings_changed)
        self.repair_ph.valueChanged.connect(self._settings_changed)
        self.proposed_changes.itemChanged.connect(self._choices_changed)
        self.repair_group.hide()
        layout.addWidget(self.repair_group)
        self.progress = QProgressBar()
        self.progress.setRange(0,100)
        self.progress.hide()
        layout.addWidget(self.progress)
        self.cancel_btn = QPushButton('Cancel operation')
        self.cancel_btn.clicked.connect(self.session.cancel)
        self.cancel_btn.hide()
        layout.addWidget(self.cancel_btn)
        self.highlights = QCheckBox('Show issue markers')
        self.highlights.setChecked(True)
        self.highlights.setToolTip('Small crosses: red errors, orange warnings, blue notes. Selection remains yellow. At most 250 markers are drawn; the issue list is authoritative.')
        self.highlights.toggled.connect(self._highlights)
        layout.addWidget(self.highlights)
        self.rebuilt_highlights = QCheckBox('Show rebuilt atoms')
        self.rebuilt_highlights.setChecked(True)
        self.rebuilt_highlights.setToolTip('Mark computationally rebuilt atoms (completed side chains and modelled loops) in cyan. Turn off to view the structure without the overlay.')
        self.rebuilt_highlights.toggled.connect(self.rebuilt_highlights_changed)
        layout.addWidget(self.rebuilt_highlights)
        self.issue_heading=QLabel('Current structure issues')
        self.issue_heading.setWordWrap(True)
        layout.addWidget(self.issue_heading)
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
        self.report_summary_btn = QPushButton('Protein report')
        self.report_summary_btn.setEnabled(False)
        self.report_summary_btn.setToolTip('Consolidated quality summary with a heuristic score. A convenience roll-up of the screens, not a validated metric.')
        self.report_summary_btn.clicked.connect(self._show_report_summary)
        layout.addWidget(self.report_summary_btn)
        self.rama_btn = QPushButton('Ramachandran plot')
        self.rama_btn.setEnabled(False)
        self.rama_btn.setToolTip('Backbone phi/psi plot with a favored/allowed/outlier region screen. Coarse region boxes, not MolProbity contours.')
        self.rama_btn.clicked.connect(self._show_ramachandran)
        layout.addWidget(self.rama_btn)
        self.session.progress_changed.connect(self._progress)
        self.session.invalidated.connect(self.invalidate)
        self.session.busy_changed.connect(self._set_busy)
        self.session.report_ready.connect(self._done)
        self.session.failed.connect(self._failed)
        self.session.cancelled.connect(lambda: self.summary.setText('Operation cancelled; the active structure is unchanged.'))
        self.session.plan_ready.connect(self._plan_ready)
        self.session.candidate_ready.connect(self._candidate_ready)
        self.c.edit_history_changed.connect(self._history)
        for widget in content.findChildren(QWidget):
            if isinstance(widget,(QPushButton,QComboBox,QDoubleSpinBox,QListWidget,QTextBrowser)):
                policy=widget.sizePolicy()
                policy.setHorizontalPolicy(QSizePolicy.Policy.Ignored)
                widget.setSizePolicy(policy)
                widget.setMinimumWidth(0)
            if isinstance(widget,QComboBox):
                widget.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
                widget.setMinimumContentsLength(8)
        content.setMinimumWidth(0)
        content.setSizePolicy(QSizePolicy.Policy.Ignored,QSizePolicy.Policy.Preferred)
        self._availability()

    def _availability(self):
        ready = self.c.source_pdb is not None and not self._busy
        self.check_btn.setEnabled(ready)
        self.validate_btn.setEnabled(ready)
        self.repair_btn.setEnabled(ready)
        checked=any(self.proposed_changes.item(i).checkState()==Qt.CheckState.Checked for i in range(self.proposed_changes.count()))
        self.build_btn.setEnabled(ready and self.repair_plan is not None and checked)
        self.accept_btn.setText('Use loop model' if self.repair_result is not None and self.repair_result.loops else 'Use repaired structure')
        self.accept_btn.setEnabled(ready and self.repair_result is not None and self.repair_result.can_accept)
        self.comparison.setEnabled(not self._busy and self.repair_result is not None)
        self.undo_btn.setEnabled(not self._busy and self._can_undo)
        result=self._exportable_result()
        self.structure_export_btn.setEnabled(not self._busy and result is not None and result.can_accept)
        for widget in (self.repair_scope,self.repair_ph,self.proposed_changes): widget.setEnabled(not self._busy)

    def invalidate(self, *_):
        self.comparison.setCurrentIndex(0)
        self.comparison_requested.emit('clear')
        self.repair_plan = self.repair_result = None
        self.candidate_summary.hide()
        self.proposed_changes.clear()
        self.repair_group.hide()
        self.report = None
        self.issues.clear()
        self.details.clear()
        self.export_btn.setEnabled(False)
        self.report_summary_btn.setEnabled(False)
        self.rama_btn.setEnabled(False)
        self.summary.setText('Structure changed. Run Check Protein on the current coordinates.')
        self._set_busy(False)
        self._highlights()

    def _set_busy(self, value):
        self._busy = value
        if value:
            self.summary.setText("Working on protein...")
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
        self.issue_heading.setText('Current structure issues')
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
        self.report_summary_btn.setEnabled(bool(report.assessment))
        self.rama_btn.setEnabled(bool(report.ramachandran))
        self._highlights()

    def _failed(self, message):
        self._set_busy(False)
        self.summary.setText('Protein operation failed. The active structure is unchanged.')
        self.details.setPlainText(str(message))

    def _inspect(self, row):
        if self.report is None or row < 0 or row >= len(self.report.issues): return
        issue = self.report.issues[row]
        self.comparison.setCurrentIndex(0)
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

    def _show_report_summary(self):
        if self.report is None or not self.report.assessment: return
        a = self.report.assessment
        m = a['metrics']
        colour = {'Excellent': '#2e7d32', 'Good': '#558b2f', 'Fair': '#e8820c', 'Poor': '#d32f2f'}.get(a['grade'], '#555')
        rows = [
            ('Protein residues', m['residues']),
            ('Ramachandran favored', f"{m['ramachandran_favored_percent']}%  ({m['ramachandran_outliers']} outlier)"),
            ('Rotamer favored', f"{m['rotamer_favored_percent']}%  ({m['rotamer_outliers']} outlier)"),
            ('Steric clashes', m['steric_clashes']),
            ('Residues missing atoms', m['residues_missing_atoms']),
            ('Unmodelled regions', m['unmodelled_regions']),
            ('Alternate conformations', m['alternate_conformations']),
            ('Disulfide bonds', m['disulfides']),
            ('Free cysteines', m['free_cysteines']),
        ]
        table = ''.join(f'<tr><td>{name}</td><td align="right"><b>{value}</b></td></tr>' for name, value in rows)
        html = (f'<div style="font-size:30px;font-weight:700;color:{colour}">{a["score"]}'
                f'<span style="font-size:14px;color:#777">/100</span>'
                f'&nbsp;<span style="font-size:16px;color:{colour}">{a["grade"]}</span></div>'
                f'<table cellpadding="4" width="100%">{table}</table>'
                f'<p style="color:#777;font-size:11px">Score = {a["score_formula"]}.</p>'
                f'<p style="color:#b00;font-size:11px">{a["note"]}</p>')
        dialog = QDialog(self)
        dialog.setWindowTitle('Protein report')
        dlayout = QVBoxLayout(dialog)
        view = QTextBrowser()
        view.setHtml(html)
        dlayout.addWidget(view)
        dialog.resize(420, 480)
        dialog.exec()

    def _show_ramachandran(self):
        if self.report is None or not self.report.ramachandran: return
        from mds.chem.ramachandran import render_svg
        s = self.report.ramachandran_summary
        dialog = QDialog(self)
        dialog.setWindowTitle('Ramachandran plot')
        dlayout = QVBoxLayout(dialog)
        caption = QLabel(f"{s['total']} residues  ·  favored {s['percent']['favored']}%  ·  "
                         f"allowed {s['percent']['allowed']}%  ·  outlier {s['percent']['outlier']}%\n"
                         "Region-based screen (general/glycine/proline); hover a point for its residue. "
                         "Coarse region boxes, not MolProbity contours.")
        caption.setWordWrap(True)
        dlayout.addWidget(caption)
        plot = QSvgWidget()
        plot.load(QByteArray(render_svg(self.report.ramachandran).encode('utf-8')))
        plot.setMinimumSize(440, 466)
        dlayout.addWidget(plot, 1)
        dialog.resize(480, 560)
        dialog.exec()

    def _export(self):
        if self.report is None: return
        path,_ = QFileDialog.getSaveFileName(self,'Export validation report','protein-health.json','JSON report (*.json)')
        if path:
            try:
                result=self._exportable_result()
                if result: result.export_report(Path(path))
                else: self.report.export(Path(path))
            except (OSError,ValueError) as exc: self.details.setPlainText(f'Cannot export report: {exc}')

    def _history(self,undo,redo):
        self._can_undo=undo
        self._availability()

    def _exportable_result(self):
        if self.repair_result is not None: return self.repair_result
        if self.accepted_result and self.c.source_pdb is not None and Path(self.accepted_result.candidate)==Path(self.c.source_pdb):
            return self.accepted_result
        return None

    def _settings_changed(self,*_):
        self.repair_plan=self.repair_result=None
        self.candidate_summary.hide()
        self.comparison.setCurrentIndex(0)
        self.comparison_requested.emit('clear')
        self.proposed_changes.clear()
        self._availability()
        self.details.setPlainText('Repair settings changed. Click Repair Protein to generate a new proposal.')

    def propose(self):
        scope=None
        if self.repair_scope.currentIndex()==1:
            source=self.c.source_pdb
            selected=self.c.selection.get(source.stem,[]) if source else []
            scope={(a['model'],a['chain_index'],a['residue_index']) for a in selected
                   if all(k in a for k in ('model','chain_index','residue_index'))}
            if not scope:
                self.repair_group.show()
                self.details.setPlainText('Select residues on the imported protein in the viewer first.')
                return
        self.session.propose(scope,self.repair_ph.value())
        self.repair_group.show()

    def _plan_ready(self,plan):
        self._done(plan.before)
        self.repair_plan=plan
        self.issue_heading.setText('Original structure issues (candidate not active)')
        self.repair_group.show()
        self._filling=True
        self.proposed_changes.clear()
        for change in plan.changes:
            item=QListWidgetItem(change.label)
            item.setData(Qt.ItemDataRole.UserRole,change.id)
            item.setFlags(item.flags()|Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if change.recommended else Qt.CheckState.Unchecked)
            self.proposed_changes.addItem(item)
        self._filling=False
        self.details.setPlainText('PROPOSED REPAIRS - no coordinates changed\nCheck the corrections to build. Hydrogens and short loops are opt-in. Build loops separately from atom completion. Choose only one conformer per residue.\n\n'+'\n'.join(plan.notes))
        self.summary.setText(f'{len(plan.changes)} proposed correction(s). Review choices before building.')
        self._availability()

    def _choices_changed(self,item):
        if self._filling or self.repair_plan is None: return
        choices={c.id:c for c in self.repair_plan.changes}
        change=choices.get(item.data(Qt.ItemDataRole.UserRole))
        if change and item.checkState()==Qt.CheckState.Checked:
            self._filling=True
            for i in range(self.proposed_changes.count()):
                other=self.proposed_changes.item(i)
                c=choices[other.data(Qt.ItemDataRole.UserRole)]
                if (change.kind=='loop') != (c.kind=='loop'):
                    other.setCheckState(Qt.CheckState.Unchecked)
            self._filling=False
        if change and change.kind=='altloc' and item.checkState()==Qt.CheckState.Checked:
            self._filling=True
            for i in range(self.proposed_changes.count()):
                other=self.proposed_changes.item(i)
                c=choices[other.data(Qt.ItemDataRole.UserRole)]
                if other is not item and c.kind=='altloc' and (c.chain_index,c.residue_index)==(change.chain_index,change.residue_index):
                    other.setCheckState(Qt.CheckState.Unchecked)
            self._filling=False
        self.repair_result=None
        self.candidate_summary.hide()
        self.comparison.setCurrentIndex(0)
        self.comparison_requested.emit('clear')
        self._availability()

    def _build(self):
        if self.repair_plan is None or self.c.project is None: return
        selected=[self.proposed_changes.item(i).data(Qt.ItemDataRole.UserRole)
                  for i in range(self.proposed_changes.count())
                  if self.proposed_changes.item(i).checkState()==Qt.CheckState.Checked]
        self.repair_result=None
        self.candidate_summary.hide()
        self.comparison.setCurrentIndex(0)
        self.comparison_requested.emit('clear')
        self.session.build(self.repair_plan,selected,self.c.project.prepared_dir('protein_repair'))

    def _candidate_ready(self,result):
        self.repair_result=result
        before_warnings=sum(i.severity=='Warning' for i in result.plan.before.issues)
        after_warnings=sum(i.severity=='Warning' for i in result.after.issues)
        self.candidate_summary.setText(
            f'Before: {before_warnings} warning(s). After: {after_warnings} warning(s).\n'
            f'{len(result.added)} new atom(s); {len(result.loops)} loop model(s) built. '+(('Geometry gates passed; backbone torsion review required.' if result.loops else 'Review required.') if result.can_accept else
            f'Blocked: {len(result.blockers)} validation finding(s).\n'+'\n'.join(result.blockers[:3])+
            ('\nSee validation details for remaining findings.' if len(result.blockers)>3 else '')))
        self.candidate_summary.show()
        self.summary.setText('Candidate ready - new atoms cyan; review required.' if result.can_accept else 'Candidate blocked by geometry checks.')
        self.details.setPlainText('BEFORE\n'+result.plan.before.summary+'\n\nAFTER\n'+result.after.summary+
            f'\n\nComputational atoms added: {len(result.added)}\nAlternate atoms removed: {len(result.removed)}\nAlternate labels resolved: {len(result.normalised)}'+
            '\n\n'+('BLOCKERS\n'+'\n'.join(result.blockers) if result.blockers else 'Supported geometry gates passed; this is not a full stereochemical validation.')+
            '\n\n'+'\n'.join(result.warnings))
        if result.loops:
            text=self.details.toPlainText()+f'\n\nLOOP MODELLING\n{result.modelling.get("backend", "")}\nSelected seed: {result.modelling.get("selected_seed")}\n'
            text+='Phi/psi angles below require independent conformational review; no statistical Ramachandran classification was performed.\n'
            for metric in result.modelling.get('geometry',[]):
                if 'residue' in metric:
                    text+=f'{metric["residue"]}: phi {metric["phi_degrees"]:.1f}, psi {metric["psi_degrees"]:.1f} degrees\n'
            self.details.setPlainText(text)
        self._availability()
        if self.isVisible(): self.comparison.setCurrentIndex(1)

    def _accept(self):
        result=self.repair_result
        if result is None: return
        self.comparison.setCurrentIndex(0)
        previous=self.accepted_result
        self.accepted_result=result
        if self.c.accept_protein_repair(result):
            self.session.check()
            self.repair_group.show()
        else: self.accepted_result=previous
        self._availability()

    def _export_structure(self):
        from mds.chem.repair import export_repaired_structure
        result=self._exportable_result()
        if result is None: return
        path,_=QFileDialog.getSaveFileName(self,'Export repaired structure','protein_repaired'+Path(result.candidate).suffix,
                                        'PDB (*.pdb);;mmCIF (*.cif *.mmcif)')
        if path:
            try: export_repaired_structure(result,Path(path))
            except (OSError,ValueError,RuntimeError) as exc: self.details.setPlainText(f'Cannot export structure: {exc}')
