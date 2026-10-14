"""Docking results: ranked pose table, contacts/interactions, exports, and preparation reports."""

from __future__ import annotations

import html
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QFileDialog, QHBoxLayout, QHeaderView, QLabel,
                               QPushButton, QSplitter, QStackedWidget, QTableWidget,
                               QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget)

from mds.chem.contacts import CAVEAT
from mds.chem.interactions import CAVEAT as INTERACTION_CAVEAT
from mds.chem.interactions import InteractionResult
from mds.chem.workflow import LigandOutcome, SCORE_DISCLAIMER
from mds.core.paths import assets_dir
from mds.gui.controller import Controller

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
    pose_selected = Signal(object, int)         # LigandOutcome, pose index

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        self._rows: list[tuple[LigandOutcome, int]] = []

        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().hide()
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._selection_changed)

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
        cl.addWidget(self.contacts)
        cl.addWidget(caveat)

        self.interactions_panel = InteractionsPanel(controller)
        side_tabs = QTabWidget()
        side_tabs.addTab(contact_box, "Distance-based contacts")
        side_tabs.addTab(self.interactions_panel, "Interactions (ProLIF)")

        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(self.table)
        split.addWidget(side_tabs)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)

        self.summary = QLabel("No docking results yet.")
        buttons = QHBoxLayout()
        buttons.addWidget(self.summary, 1)
        for text, fmt in (("Export CSV…", "csv"), ("Export Excel…", "xlsx"), ("Export SDF…", "sdf")):
            b = QPushButton(text)
            b.clicked.connect(lambda _=False, f=fmt: self._export(f))
            buttons.addWidget(b)
        open_btn = QPushButton("Open run folder")
        open_btn.clicked.connect(self._open_folder)
        buttons.addWidget(open_btn)

        disclaimer = QLabel(f"<small><i>{SCORE_DISCLAIMER}</i></small>")
        disclaimer.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(buttons)
        layout.addWidget(split, 1)
        layout.addWidget(disclaimer)
        controller.docking_finished.connect(self.show_outcomes)

    def show_outcomes(self, outcomes: list[LigandOutcome]) -> None:
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self._rows = []
        failed = []
        for o in outcomes:
            if o.status != "ok":
                failed.append(f"{o.name}: {o.status} {o.message}")
                continue
            for i, (p, r) in enumerate(zip(o.poses, o.reference_rmsd)):
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
        text = f"{len(ok)} ligand(s) docked"
        refs = [o.reference_rmsd[0] for o in ok if o.reference_rmsd and o.reference_rmsd[0] is not None]
        if refs:
            verdict = "SUCCESS (<2 Å)" if refs[0] < 2.0 else "NOT reproduced (≥2 Å)"
            text += f" · redocking top-pose RMSD {refs[0]:.2f} Å — {verdict}"
        if failed:
            text += f" · {len(failed)} failed: " + "; ".join(failed)
        self.summary.setText(text)
        if self._rows:
            self.table.selectRow(0)

    def _selection_changed(self) -> None:
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return
        idx = self.table.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        outcome, pose_i = self._rows[idx]
        self._show_contacts(outcome, pose_i)
        self.interactions_panel.set_pose(outcome, pose_i)
        self.pose_selected.emit(outcome, pose_i)

    def _show_contacts(self, outcome: LigandOutcome, pose_i: int) -> None:
        try:
            contacts = self.c.contacts_for(outcome, pose_i)
        except Exception as exc:  # noqa: BLE001
            self.contacts.setRowCount(0)
            self.summary.setText(f"Contact analysis unavailable: {exc}")
            return
        self.contacts.setRowCount(len(contacts))
        for r, c in enumerate(contacts):
            self.contacts.setItem(r, 0, QTableWidgetItem(c.residue))
            self.contacts.setItem(r, 1, QTableWidgetItem(f"{c.min_distance:.2f}"))
            self.contacts.setItem(r, 2, QTableWidgetItem("; ".join(c.candidates)))

    def _export(self, fmt: str) -> None:
        if not self.c.outcomes:
            self.c.error.emit("Nothing to export", "Run a docking calculation first.")
            return
        project = self.c.ensure_project()
        default = str(project.root / "exports" / f"docking_results.{fmt}")
        path, _ = QFileDialog.getSaveFileName(self, "Export results", default, f"*.{fmt}")
        if not path:
            return
        {"csv": self.c.export_csv, "xlsx": self.c.export_excel, "sdf": self.c.export_sdf}[fmt](Path(path))
        self.c.message.emit(f"Exported {path}")

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
