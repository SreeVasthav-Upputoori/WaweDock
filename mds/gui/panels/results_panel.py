"""Docking results: ranked pose table, contacts/interactions, exports, and preparation reports."""

from __future__ import annotations

import html
from pathlib import Path

from PySide6.QtCore import QItemSelectionModel, QSignalBlocker, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QFileDialog, QHBoxLayout, QHeaderView, QLabel,
                               QLineEdit, QMenu, QPushButton, QSplitter, QStackedWidget, QTableWidget,
                               QTableWidgetItem, QTabWidget, QTextEdit, QToolButton, QVBoxLayout, QWidget)

from mds.chem.contacts import CAVEAT
from mds.chem.interactions import CAVEAT as INTERACTION_CAVEAT
from mds.chem.interactions import InteractionResult
from mds.chem.workflow import LigandOutcome, SCORE_DISCLAIMER
from mds.core.paths import assets_dir
from mds.gui.controller import Controller
from mds.gui.panels.score_chart import ScoreChart
from mds.gui.panels.pose_review import PoseReviewBar
from mds.gui.panels.pharmacophore_panel import PharmacophorePanel

COLUMNS = ["Ligand", "Pose", "Score (kcal/mol)", "RMSD lb", "RMSD ub", "RMSD to crystal (Å)", "Conformer"]


class _NumItem(QTableWidgetItem):
    """Sorts numerically instead of alphabetically."""

    def __init__(self, value: float | None, text: str) -> None:
        super().__init__(text)
        self._v = float("inf") if value is None else value

    def __lt__(self, other: QTableWidgetItem) -> bool:  # noqa: D105
        return self._v < getattr(other, "_v", float("inf"))


class InteractionsPanel(QWidget):
    """ProLIF-classified interactions for the selected pose: table + 2D network diagram."""

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        self._current: tuple[LigandOutcome, int] | None = None

        self.analyse_btn = QPushButton("Analyse Interactions (ProLIF)")
        self.analyse_btn.setToolTip("Classifies contacts as hydrogen bonds, hydrophobic contacts, "
                                    "pi-stacking, salt bridges, etc. using distance AND angle rules.")
        self.analyse_btn.clicked.connect(self._analyse)
        self.status = QLabel("Select a pose, then click Analyse.")
        self.status.setWordWrap(True)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Residue", "Interaction type(s)"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)

        self._diagram_stack = QStackedWidget()
        self._placeholder = QLabel("<i>No diagram yet.</i>")
        self._placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._diagram_stack.addWidget(self._placeholder)
        self._diagram = None   # created lazily: QtWebEngine import is somewhat heavy

        caveat = QLabel(f"<small>{INTERACTION_CAVEAT}</small>")
        caveat.setWordWrap(True)

        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.table)
        split.addWidget(self._diagram_stack)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 2)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        bar = QHBoxLayout()
        bar.addWidget(self.analyse_btn)
        bar.addWidget(self.status, 1)
        layout.addLayout(bar)
        layout.addWidget(split, 1)
        layout.addWidget(caveat)

        controller.interactions_started.connect(self._started)
        controller.interactions_ready.connect(self._ready)
        controller.runner.job_failed.connect(lambda job_id, _tb: self._job_stopped(job_id, "failed"))
        controller.runner.job_cancelled.connect(lambda job_id: self._job_stopped(job_id, "cancelled"))

    def clear_pose(self) -> None:
        self._current = None
        self.table.setRowCount(0)
        self.status.setText("Select a pose, then click Analyse.")
        self._diagram_stack.setCurrentWidget(self._placeholder)
        self._placeholder.setText("<i>No diagram yet.</i>")

    def _job_stopped(self, job_id: int, status: str) -> None:
        try:
            is_interaction = self.c.runner.name(job_id).startswith("Interaction analysis")
        except KeyError:
            return
        if is_interaction:
            self.analyse_btn.setEnabled(True)
            self.status.setText(f"Interaction analysis {status}. See Log for details; you can retry.")

    def set_pose(self, outcome: LigandOutcome, pose_index: int) -> None:
        self._current = (outcome, pose_index)
        self.table.setRowCount(0)
        self._diagram_stack.setCurrentWidget(self._placeholder)
        self._placeholder.setText("<i>No diagram yet — click Analyse Interactions.</i>")
        self.status.setText(f"Ready to analyse {outcome.name} pose #{outcome.poses[pose_index].rank}.")

    def _analyse(self) -> None:
        if self._current is None:
            self.status.setText("Select a pose in the results table first.")
            return
        self.c.analyze_interactions(*self._current)

    def _started(self) -> None:
        self.analyse_btn.setEnabled(False)
        self.status.setText("Analysing… (building protein model, running ProLIF)")

    def _ready(self, outcome: LigandOutcome, pose_index: int, result: InteractionResult) -> None:
        self.analyse_btn.setEnabled(True)
        if self._current != (outcome, pose_index):
            return  # selection moved on while this was running; discard
        self.status.setText(f"{result.n_interactions} interaction(s) across {len(result.residues)} "
                            f"residue(s) for {outcome.name} pose #{outcome.poses[pose_index].rank}.")
        self.table.setRowCount(len(result.residues))
        for r, ri in enumerate(result.residues):
            self.table.setItem(r, 0, QTableWidgetItem(ri.residue))
            self.table.setItem(r, 1, QTableWidgetItem(", ".join(ri.interactions)))
        if result.network_html:
            self._show_diagram(result.network_html)
        else:
            self._placeholder.setText("<i>No 2D diagram available for this pose "
                                      "(no interactions detected, or the diagram could not be built).</i>")
            self._diagram_stack.setCurrentWidget(self._placeholder)

    def _show_diagram(self, html_text: str) -> None:
        if self._diagram is None:
            from PySide6.QtWebEngineWidgets import QWebEngineView
            self._diagram = QWebEngineView()
            self._diagram_stack.addWidget(self._diagram)
        from PySide6.QtCore import QUrl
        base = QUrl.fromLocalFile(str(assets_dir() / "viewer") + "/")
        self._diagram.setHtml(html_text, base)
        self._diagram_stack.setCurrentWidget(self._diagram)


class ResultsPanel(QWidget):
    pose_cleared = Signal()
    pose_selected = Signal(object, int)         # LigandOutcome, pose index

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        controller.contacts_ready.connect(self._contacts_ready)
        controller.contacts_failed.connect(self._contacts_failed)
        self._rows: list[tuple[LigandOutcome, int]] = []
        self._current_pose = None
        self.review_bar = PoseReviewBar(self)
        self.review_bar.step_requested.connect(self._step_pose)

        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().hide()
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setToolTip("Select a pose to view it. Ctrl/Shift-select rows to compare them in Compare Poses.")
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._selection_changed)
        self.table.currentCellChanged.connect(lambda *_args: self._selection_changed())

        self.contact_status=QLabel("Select a pose to preview contacts.")
        self.contact_status.setTextFormat(Qt.TextFormat.PlainText)
        self.contact_status.setWordWrap(True)
        self.contacts = QTableWidget(0, 3)
        self.contacts.setHorizontalHeaderLabels(["Residue", "Min. distance (Å)", "Distance-based candidates"])
        self.contacts.verticalHeader().hide()
        self.contacts.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.contacts.horizontalHeader().setStretchLastSection(True)
        self.contacts.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.contacts.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        caveat = QLabel(f"<small>{CAVEAT}</small>")
        caveat.setWordWrap(True)
        contact_box = QWidget()
        cl = QVBoxLayout(contact_box)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.addWidget(self.contact_status)
        cl.addWidget(self.contacts)
        cl.addWidget(caveat)

        self.interactions_panel = InteractionsPanel(controller)
        self.pharmacophore_panel = PharmacophorePanel(controller)
        self.pharmacophore_panel.pose_requested.connect(self.select_pose)
        self.analysis_tabs = QTabWidget()
        self.analysis_tabs.setUsesScrollButtons(True)
        self.analysis_tabs.setElideMode(Qt.TextElideMode.ElideRight)
        self.analysis_tabs.addTab(contact_box, "Contacts")
        self.analysis_tabs.addTab(self.interactions_panel, "Interactions")
        self.analysis_tabs.addTab(self.pharmacophore_panel, "Pharmacophore")
        self.comparison = QTableWidget(0, len(COLUMNS))
        self.comparison.setHorizontalHeaderLabels(COLUMNS)
        self.comparison.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.comparison.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.analysis_tabs.addTab(self.comparison, "Compare Poses")
        self.score_chart = ScoreChart()
        self.score_chart.set_outcomes([])
        self.score_chart.pose_selected.connect(self.select_pose)
        self.analysis_tabs.addTab(self.score_chart, "Score Comparison")

        self.splitter = QSplitter(Qt.Orientation.Vertical)
        self.splitter.addWidget(self.table)
        self.splitter.addWidget(self.analysis_tabs)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)

        self.summary = QLabel("No docking results yet.")
        self.summary.setWordWrap(True)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search ligands, poses or scores...")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._filter_rows)
        self.match_count = QLabel("0 poses")
        self.export_button = QToolButton()
        self.export_button.setText("Export")
        self.export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        export_menu = QMenu(self.export_button)
        for text, fmt in (("CSV", "csv"), ("Excel", "xlsx"), ("SDF", "sdf"), ("HTML Report", "html")):
            export_menu.addAction(text, lambda _on=False, f=fmt: self._export(f))
        self.export_button.setMenu(export_menu)
        self.export_button.setToolTip("Export all results, including poses hidden by search.")
        open_btn = QToolButton()
        open_btn.setText("Run Folder")
        open_btn.clicked.connect(self._open_folder)
        buttons = QHBoxLayout()
        buttons.addWidget(self.search, 1)
        buttons.addWidget(self.export_button)
        buttons.addWidget(open_btn)

        disclaimer = QLabel(f"<small><i>{SCORE_DISCLAIMER}</i></small>")
        disclaimer.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.summary)
        layout.addLayout(buttons)
        layout.addWidget(self.review_bar)
        layout.addWidget(self.match_count)
        layout.addWidget(self.splitter, 1)
        layout.addWidget(disclaimer)
        controller.docking_finished.connect(self.show_outcomes)

    def show_outcomes(self, outcomes: list[LigandOutcome]) -> None:
        blocker = QSignalBlocker(self.table)
        self.table.setSortingEnabled(False)
        self.contacts.setRowCount(0)
        self.interactions_panel.clear_pose()
        self.pharmacophore_panel.clear()
        self.comparison.setRowCount(0)
        self.table.setRowCount(0)
        self._rows = []
        self._current_pose = None
        failed = []
        for o in outcomes:
            if o.status != "ok":
                failed.append(f"{o.name}: {o.status} {o.message}")
                continue
            for i, p in enumerate(o.poses):
                r = o.reference_rmsd[i] if i < len(o.reference_rmsd) else None
                row = self.table.rowCount()
                self.table.insertRow(row)
                cells = [QTableWidgetItem(o.name), _NumItem(p.rank, str(p.rank)),
                         _NumItem(p.score, f"{p.score:.3f}"), _NumItem(p.rmsd_lb, f"{p.rmsd_lb:.3f}"),
                         _NumItem(p.rmsd_ub, f"{p.rmsd_ub:.3f}"),
                         _NumItem(r, "" if r is None else f"{r:.2f}"), _NumItem(p.conformer + 1, str(p.conformer + 1))]
                cells[0].setData(Qt.ItemDataRole.UserRole, len(self._rows))
                for col, item in enumerate(cells):
                    self.table.setItem(row, col, item)
                if r is not None and r < 2.0 and p.rank == 1:
                    for item in cells:
                        item.setBackground(Qt.GlobalColor.darkGreen)
                        item.setForeground(Qt.GlobalColor.white)
                self._rows.append((o, i))
        self.table.setSortingEnabled(True)
        self.table.sortItems(2, Qt.SortOrder.AscendingOrder)
        ok = [o for o in outcomes if o.status == "ok"]
        text = f"{len(ok)} ligand(s) docked" if outcomes else "No docking results yet."
        refs = [o.reference_rmsd[0] for o in ok if o.reference_rmsd and o.reference_rmsd[0] is not None]
        if refs:
            verdict = "SUCCESS (<2 Å)" if refs[0] < 2.0 else "NOT reproduced (≥2 Å)"
            text += f" · redocking top-pose RMSD {refs[0]:.2f} Å — {verdict}"
        if failed:
            text += f" · {len(failed)} failed: " + "; ".join(failed)
        self.summary.setText(text)
        self.score_chart.set_outcomes(outcomes)
        del blocker
        self._filter_rows()

    def _filter_rows(self, *_args) -> None:
        query = self.search.text().strip().casefold()
        blocker = QSignalBlocker(self.table)
        visible = []
        for row in range(self.table.rowCount()):
            matches = not query or any(query in self.table.item(row, col).text().casefold()
                                       for col in range(self.table.columnCount()))
            self.table.setRowHidden(row, not matches)
            if matches:
                visible.append(row)
            else:
                self.table.selectionModel().select(self.table.model().index(row, 0),
                    QItemSelectionModel.SelectionFlag.Deselect | QItemSelectionModel.SelectionFlag.Rows)
        selected = self.table.selectionModel().selectedRows()
        if visible and not selected:
            self.table.selectRow(visible[0])
        elif not visible:
            self.table.clearSelection()
        del blocker
        self.match_count.setText(f"{len(visible)} / {self.table.rowCount()} poses")
        self._selection_changed()

    def select_pose(self, outcome: LigandOutcome, pose_index: int) -> None:
        self.search.clear()
        for row in range(self.table.rowCount()):
            idx = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
            candidate, index = self._rows[idx]
            if candidate is outcome and index == pose_index:
                self.table.selectRow(row)
                self.table.scrollToItem(self.table.item(row, 0))
                return

    def _step_pose(self, direction: int) -> None:
        rows = [row for row in range(self.table.rowCount()) if not self.table.isRowHidden(row)]
        if not rows:
            return
        current = self.table.currentRow()
        index = rows.index(current) if current in rows else 0
        target = max(0, min(len(rows)-1, index+direction))
        self.table.selectRow(rows[target])
        self.table.scrollToItem(self.table.item(rows[target],0))

    def _selection_changed(self) -> None:
        rows = [index.row() for index in self.table.selectionModel().selectedRows()
                if not self.table.isRowHidden(index.row())]
        self.comparison.setRowCount(len(rows))
        for target, source in enumerate(rows):
            for col in range(len(COLUMNS)):
                self.comparison.setItem(target, col, QTableWidgetItem(self.table.item(source, col).text()))
        if not rows:
            self.review_bar.set_position(-1, sum(not self.table.isRowHidden(r) for r in range(self.table.rowCount())))
            self._current_pose = None
            self.c.clear_contact_requests()
            self.contact_status.setText("Select a pose to preview contacts.")
            self.contacts.setRowCount(0)
            self.interactions_panel.clear_pose()
            self.pharmacophore_panel.clear_pose()
            self.pose_cleared.emit()
            return
        row = self.table.currentRow()
        row = row if row in rows else rows[0]
        visible_rows = [r for r in range(self.table.rowCount()) if not self.table.isRowHidden(r)]
        idx = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        outcome, pose_i = self._rows[idx]
        pose = outcome.poses[pose_i]
        self.review_bar.set_position(visible_rows.index(row), len(visible_rows),
                                     f"{outcome.name} #{pose.rank} · {pose.score:.2f} kcal/mol")
        if self._current_pose is not None and self._current_pose[0] is outcome and self._current_pose[1] == pose_i:
            return
        self._current_pose = (outcome, pose_i)
        self._show_contacts(outcome, pose_i)
        self.interactions_panel.set_pose(outcome, pose_i)
        self.pharmacophore_panel.set_pose(outcome, pose_i)
        self.pose_selected.emit(outcome, pose_i)

    def _show_contacts(self, outcome:LigandOutcome,pose_i:int) -> None:
        self.contacts.setRowCount(0)
        self.contact_status.setText("Loading contacts...")
        self.contacts.setToolTip("Loading contact preview...")
        self.c.request_contacts(outcome,pose_i)

    def _contacts_ready(self,outcome,pose_i,contacts) -> None:
        if self._current_pose is None or self._current_pose[0] is not outcome or self._current_pose[1]!=pose_i:
            return
        self.contact_status.setText(f"{len(contacts)} contacting residue(s)." if self.c.receptor is not None
                                    else "Prepare a receptor to preview contacts.")
        self.contacts.setToolTip(CAVEAT)
        self.contacts.setUpdatesEnabled(False)
        try:
            self.contacts.setRowCount(len(contacts))
            for row,contact in enumerate(contacts):
                self.contacts.setItem(row,0,QTableWidgetItem(contact.residue))
                self.contacts.setItem(row,1,QTableWidgetItem(f"{contact.min_distance:.2f}"))
                self.contacts.setItem(row,2,QTableWidgetItem("; ".join(contact.candidates)))
        finally:
            self.contacts.setUpdatesEnabled(True)

    def _contacts_failed(self,outcome,pose_i,message) -> None:
        if self._current_pose is not None and self._current_pose[0] is outcome and self._current_pose[1]==pose_i:
            self.contact_status.setText(f"Contact preview unavailable: {message}")
            self.contacts.setToolTip(f"Contact preview unavailable: {message}")

    def _export(self, fmt: str) -> None:
        if not self.c.outcomes:
            self.c.error.emit("Nothing to export", "Run a docking calculation first.")
            return
        project = self.c.ensure_project()
        default = str(project.root / "exports" / f"docking_results.{fmt}")
        path, _ = QFileDialog.getSaveFileName(self, "Export results", default, f"*.{fmt}")
        if not path:
            return
        if fmt == "html":
            self.export_report(Path(path))
        else:
            {"csv": self.c.export_csv, "xlsx": self.c.export_excel, "sdf": self.c.export_sdf}[fmt](Path(path))
        self.c.message.emit(f"Exported {path}")

    def export_report(self, path: Path) -> None:
        """Export all pose metrics and run failures without fabricating analysis."""
        rows = []
        for outcome, i in self._rows:
            pose = outcome.poses[i]
            rmsd = outcome.reference_rmsd[i] if i < len(outcome.reference_rmsd) else None
            values = [outcome.name, pose.rank, pose.score, pose.rmsd_lb, pose.rmsd_ub,
                      "" if rmsd is None else rmsd, pose.conformer + 1]
            rows.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in values) + "</tr>")
        failures = [f"{o.name}: {o.status} {o.message}" for o in self.c.outcomes if o.status != "ok"]
        headers = "".join(f"<th>{html.escape(column)}</th>" for column in COLUMNS)
        document = ('<!doctype html><meta charset="utf-8"><title>WaweDock Results</title>'
                    '<style>body{font-family:system-ui;margin:2rem}table{border-collapse:collapse}'
                    'td,th{padding:.5rem;border:1px solid #ccc;text-align:left}</style>'
                    '<h1>WaweDock Results</h1><p>' + html.escape(SCORE_DISCLAIMER) + '</p>'
                    '<table><thead><tr>' + headers + '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table>'
                    '<h2>Run failures</h2><p>' + html.escape('; '.join(failures) or 'None reported.') + '</p>')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(document, encoding="utf-8")

    def _open_folder(self) -> None:
        if self.c.run_dir:
            from PySide6.QtCore import QUrl
            from PySide6.QtGui import QDesktopServices
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.c.run_dir)))


class ReportPanel(QWidget):
    """Shows preparation reports; warnings and errors are colour-coded."""

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(False)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.tabs)
        controller.report.connect(self.show_report)
        self._views: dict[str, QTextEdit] = {}

    def show_report(self, title: str, text: str) -> None:
        view = self._views.get(title)
        if view is None:
            view = QTextEdit()
            view.setReadOnly(True)
            self._views[title] = view
            self.tabs.addTab(view, title)
        colours = {"[WARNING]": "#d08a00", "[ERROR  ]": "#e05050"}
        lines = []
        for line in text.splitlines():
            esc = html.escape(line)
            colour = next((c for k, c in colours.items() if line.startswith(k)), None)
            lines.append(f'<span style="color:{colour}">{esc}</span>' if colour else esc)
        view.setHtml("<pre style='font-family:Consolas,monospace'>" + "<br>".join(lines) + "</pre>")
        self.tabs.setCurrentWidget(view)
