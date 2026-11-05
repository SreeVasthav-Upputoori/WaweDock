"""Pose consistency tab: compare every docked pose with one reference binding mode."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QFileDialog, QFrame, QHBoxLayout,
                               QHeaderView, QLabel, QMenu, QPushButton, QSplitter, QTableWidget,
                               QTableWidgetItem, QTextEdit, QToolButton, QVBoxLayout, QWidget)

from mds.chem.pose_consistency import (CAVEAT, CONSERVED, DIFFERENT, FLIPPED, METHOD_WEIGHTS,
                                       OUTLIER, PARTIAL, VERDICT_THRESHOLDS, ConsistencyOptions,
                                       rows)
from mds.gui.panels.common import CaveatBar

COLUMNS = [
    ("Pose", "The docked pose being compared. Select a row to view it."),
    ("Verdict", "Conserved / Partially conserved / Different / Flipped / Outlier, from the "
                "consensus score and the flip test."),
    ("Consensus", "Weighted mean of the axes that could be computed: "
                  + ", ".join(f"{name} {weight:.0%}" for name, weight in METHOD_WEIGHTS.items())
                  + f". Conserved at ≥ {VERDICT_THRESHOLDS['conserved']:.2f}, partial at "
                    f"≥ {VERDICT_THRESHOLDS['partial']:.2f}. These weights and cut-offs are "
                    "conventions for ordering the list, not a fitted model."),
    ("RMSD (Å)", "Heavy-atom RMSD to the reference, symmetry-aware and measured where the "
                      "poses sit — nothing is re-superposed. Blank when the ligand differs "
                      "from the reference."),
    ("MCS RMSD (Å)", "RMSD over the largest substructure the pose and the reference share, "
                          "so a different ligand can still be compared on the common scaffold."),
    ("MCS", "How much of the reference the shared substructure covers."),
    ("Interaction Tc", "Tanimoto over the (residue, interaction type) pairs ProLIF classifies, "
                       "van der Waals contacts excluded. Needs the optional Open Babel backend."),
    ("Features", "Reference pharmacophore features the pose reproduces within tolerance, "
                 "matched only within their own kind."),
    ("Centroid Δ (Å)", "Distance between the two poses' heavy-atom centroids."),
    ("Axis °", "Angle between the two long axes, folded into 0–90°. Blank when a "
                    "ligand is too isotropic for its long axis to mean anything."),
    ("Volume (%)", "Overlap of the two poses' van der Waals volumes where they actually sit."),
    ("USR", "Shape similarity of the conformers. Position-independent by construction, so it is "
            "reported but deliberately kept out of the consensus."),
]

_VERDICT_COLORS = {
    CONSERVED: "#1e7a3c", PARTIAL: "#8a6d00", FLIPPED: "#8a3fb0",
    DIFFERENT: "#8a4b00", OUTLIER: "#7a2222",
}


def _separator() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.VLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line


class _NumItem(QTableWidgetItem):
    """Sorts numerically while displaying formatted text; blanks sort last."""

    def __init__(self, value: float | None, text: str) -> None:
        super().__init__(text)
        self._v = float("-inf") if value is None else float(value)

    def __lt__(self, other: QTableWidgetItem) -> bool:  # noqa: D105
        return self._v < getattr(other, "_v", float("-inf"))


class ConsistencyPanel(QWidget):
    """Ranks docked poses by how much of a reference binding mode they reproduce.

    The reference is chosen explicitly — the co-crystal ligand, or whichever docked pose the
    user decides to trust — because the comparison is only as meaningful as that choice, and
    guessing it would hide the most important assumption in the analysis.
    """

    pose_requested = Signal(object, int)         # LigandOutcome, pose index

    def __init__(self, controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        self._current: tuple[object, int] | None = None
        self._pose_reference: tuple[object, int] | None = None
        self._rows: list[tuple[object, int]] = []
        self._results: list = []
        self._awaiting = False

        self.reference = QComboBox()
        self.reference.setToolTip("The trusted pose every other pose is measured against.")
        self.reference.currentIndexChanged.connect(lambda _i: self._reference_changed())
        self.use_pose_btn = QPushButton("Use selected pose")
        self.use_pose_btn.setToolTip("Make the pose selected in the results table the reference.")
        self.use_pose_btn.clicked.connect(self._use_selected_pose)
        self.use_pose_btn.setEnabled(False)

        # The three axis switches are set once for a run and then left alone, so they sit behind
        # a button rather than taking a standing row above the results they produce.
        self.axes_btn = QToolButton()
        self.axes_btn.setText("Axes")
        self.axes_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        axes_menu = QMenu(self.axes_btn)
        self.interactions = QAction("Interaction fingerprint (ProLIF)", axes_menu,
                                    checkable=True, checked=True)
        self.interactions.setToolTip("Classify each pose's contacts and compare the sets. The "
                                     "slowest axis, and the only one needing Open Babel.")
        self.shape = QAction("Shape and volume overlap", axes_menu, checkable=True, checked=True)
        self.features = QAction("Pharmacophore features", axes_menu, checkable=True, checked=True)
        for action in (self.interactions, self.shape, self.features):
            axes_menu.addAction(action)
            action.toggled.connect(lambda _on: self._axes_changed())
        self.axes_btn.setMenu(axes_menu)
        self._axes_changed()

        self.compare_btn = QPushButton("Compare all poses")
        self.compare_btn.clicked.connect(self._compare)
        self.compare_btn.setEnabled(False)
        self.export_btn = QPushButton("Export CSV")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self._export)

        self.status = QLabel("Choose a reference pose, then compare.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)

        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels([name for name, _tip in COLUMNS])
        for column, (_name, tip) in enumerate(COLUMNS):
            self.table.horizontalHeaderItem(column).setToolTip(tip)
        self.table.verticalHeader().hide()
        self.table.setSortingEnabled(True)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setToolTip("Select a row to view that pose and read its per-pose notes.")
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._row_selected)

        self.detail = QTextEdit()
        self.detail.setReadOnly(True)
        self.detail.setPlaceholderText("Select a pose to see which reference contacts it keeps, "
                                       "loses and gains.")
        # Hidden until there is something in it: an empty pane is 75px taken from the table.
        self.detail.setVisible(False)

        label = QLabel("Reference:")
        label.setBuddy(self.reference)
        controls = QHBoxLayout()
        controls.setSpacing(4)
        for widget in (label, self.reference, self.use_pose_btn, _separator(), self.axes_btn,
                       self.compare_btn, self.export_btn):
            controls.addWidget(widget)
        controls.setStretchFactor(self.reference, 1)
        controls.addStretch(0)

        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.table)
        split.addWidget(self.detail)
        split.setStretchFactor(0, 4)
        split.setStretchFactor(1, 1)
        self.detail.setMaximumHeight(140)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 2, 0, 0)
        layout.setSpacing(3)
        layout.addLayout(controls)
        layout.addWidget(self.status)
        layout.addWidget(split, 1)
        layout.addWidget(CaveatBar(CAVEAT))

        controller.consistency_busy.connect(self._busy)
        controller.consistency_ready.connect(self._ready)
        controller.reference_changed.connect(lambda _name: self.refresh_references())
        controller.docking_finished.connect(lambda _outcomes: self.refresh_references())
        self.refresh_references()

    def _axes_changed(self) -> None:
        """Name the optional axes on the button, so a switched-off axis is never a hidden one."""
        off = [action.text().split(" (")[0] for action in
               (self.interactions, self.shape, self.features) if not action.isChecked()]
        self.axes_btn.setText("All axes" if not off else f"Axes ({3 - len(off)}/3)")
        self.axes_btn.setToolTip("Every axis is included."
                                 if not off else "Switched off: " + ", ".join(off))

    # ------------------------------------------------------------ reference choice
    def refresh_references(self) -> None:
        """Offer the co-crystal ligand and, once chosen, the pose being used as reference."""
        previous = self.reference.currentData()
        self.reference.blockSignals(True)
        self.reference.clear()
        if self.c.reference is not None:
            self.reference.addItem(f"Co-crystal ligand — {self.c.reference_name}", "crystal")
        if self._pose_reference is not None:
            outcome, index = self._pose_reference
            self.reference.addItem(
                f"Docked pose — {outcome.name} #{outcome.poses[index].rank}", "pose")
        if previous is not None:
            position = self.reference.findData(previous)
            if position >= 0:
                self.reference.setCurrentIndex(position)
        self.reference.blockSignals(False)
        self._reference_changed()

    def _reference_changed(self) -> None:
        has_reference = self.reference.count() > 0
        self.reference.setEnabled(has_reference)
        self.compare_btn.setEnabled(has_reference and bool(self.c.outcomes))
        if not has_reference:
            self.status.setText(
                "No reference pose yet. Prepare a receptor that has a co-crystal ligand, or "
                "select a docked pose you trust and press “Use selected pose”.")

    def _use_selected_pose(self) -> None:
        if self._current is None:
            return
        self._pose_reference = self._current
        self.refresh_references()
        position = self.reference.findData("pose")
        if position >= 0:
            self.reference.setCurrentIndex(position)

    def _reference_molecule(self) -> tuple[object, str] | tuple[None, str]:
        kind = self.reference.currentData()
        if kind == "crystal":
            return self.c.reference, self.c.reference_name or "co-crystal ligand"
        if kind == "pose" and self._pose_reference is not None:
            outcome, index = self._pose_reference
            return outcome.mols[index], f"{outcome.name} #{outcome.poses[index].rank}"
        return None, ""

    # -------------------------------------------------------------------- pose state
    def set_pose(self, outcome, pose_index: int) -> None:
        self._current = (outcome, pose_index)
        self.use_pose_btn.setEnabled(True)

    def clear_pose(self) -> None:
        self._current = None
        self.use_pose_btn.setEnabled(False)

    def clear(self) -> None:
        """Forget the results as well as the reference pose — used when results are replaced."""
        self.clear_pose()
        self._pose_reference = None
        self._rows, self._results = [], []
        self._awaiting = False
        self.table.setRowCount(0)
        self.detail.clear()
        self.detail.setVisible(False)
        self.export_btn.setEnabled(False)
        self.status.setText("Choose a reference pose, then compare.")
        self.refresh_references()       # last, so it can replace that with the no-reference hint

    # --------------------------------------------------------------------- comparing
    def _compare(self) -> None:
        reference, label = self._reference_molecule()
        if reference is None:
            self.status.setText("Choose a reference pose first.")
            return
        poses, self._rows = [], []
        for outcome in self.c.outcomes:
            if outcome.status != "ok":
                continue
            for index, mol in enumerate(outcome.mols):
                if (self._pose_reference is not None and self._pose_reference[0] is outcome
                        and self._pose_reference[1] == index):
                    continue        # the reference would trivially score 100% against itself
                poses.append(((outcome.name, outcome.poses[index].rank),
                              f"{outcome.name} #{outcome.poses[index].rank}", mol))
                self._rows.append((outcome, index))
        if not poses:
            self.status.setText("No other poses to compare against this reference.")
            return
        options = ConsistencyOptions(interactions=self.interactions.isChecked(),
                                     shape=self.shape.isChecked(),
                                     features=self.features.isChecked())
        self.status.setText(f"Comparing {len(poses)} pose(s) with {label}…")
        self.c.analyze_pose_consistency(reference, label, poses, options)

    def _busy(self, running: bool) -> None:
        # The controller emits this exactly when a job is submitted, so it -- not the button
        # press -- is what marks a comparison as genuinely pending.
        self._awaiting = self._awaiting or running
        self.compare_btn.setEnabled(not running and self.reference.count() > 0)
        if not running and self._awaiting:
            # _ready clears the flag, so reaching here means the job failed or was cancelled.
            self._awaiting = False
            self.status.setText("Comparison stopped before any result arrived — see the Log, "
                                "then try again.")

    def _ready(self, reference_label: str, results: list) -> None:
        if len(results) != len(self._rows):
            return      # the results changed underneath a running comparison
        self._awaiting = False
        self._results = [result for _key, result in results]
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(results))
        for row, (_key, r) in enumerate(results):
            cells = [
                QTableWidgetItem(r.label),
                QTableWidgetItem(r.verdict),
                _NumItem(r.consensus, f"{r.consensus:.2f}"),
                _NumItem(r.rmsd, "" if r.rmsd is None else f"{r.rmsd:.2f}"),
                _NumItem(r.mcs_rmsd, "" if r.mcs_rmsd is None else f"{r.mcs_rmsd:.2f}"),
                _NumItem(r.mcs_fraction, f"{r.mcs_atoms} atoms ({r.mcs_fraction:.0%})"),
                _NumItem(r.ifp_tanimoto,
                         "—" if r.ifp_tanimoto is None else f"{r.ifp_tanimoto:.2f}"),
                _NumItem(r.feature_recovery,
                         "—" if r.feature_recovery is None
                         else f"{r.matched_features} / {r.reference_features}"),
                _NumItem(r.centroid_shift,
                         "" if r.centroid_shift is None else f"{r.centroid_shift:.2f}"),
                _NumItem(r.axis_angle, "—" if r.axis_angle is None else f"{r.axis_angle:.0f}"),
                _NumItem(r.volume_overlap,
                         "—" if r.volume_overlap is None else f"{r.volume_overlap:.0f}"),
                _NumItem(r.usr, "—" if r.usr is None else f"{r.usr:.2f}"),
            ]
            cells[0].setData(Qt.ItemDataRole.UserRole, row)
            colour = _VERDICT_COLORS.get(r.verdict)
            if colour:
                cells[1].setForeground(QColor(colour))
            if r.notes:
                for item in cells:
                    item.setToolTip(" ".join(r.notes))
            for column, item in enumerate(cells):
                self.table.setItem(row, column, item)
        self.table.setSortingEnabled(True)
        self.table.sortItems(2, Qt.SortOrder.DescendingOrder)
        counts: dict[str, int] = {}
        for r in self._results:
            counts[r.verdict] = counts.get(r.verdict, 0) + 1
        summary = ", ".join(f"{n} {verdict.lower()}" for verdict, n in sorted(counts.items()))
        self.status.setText(f"{len(results)} pose(s) compared with {reference_label}: {summary}. "
                            f"The consensus column orders the list; read the individual axes "
                            f"before drawing a conclusion.")
        self.export_btn.setEnabled(True)

    def _row_selected(self) -> None:
        selected = self.table.selectionModel().selectedRows()
        if not selected:
            return
        index = self.table.item(selected[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        if index is None or index >= len(self._results):
            return
        self._show_detail(self._results[index])
        outcome, pose_index = self._rows[index]
        self.pose_requested.emit(outcome, pose_index)

    def _show_detail(self, result) -> None:
        def block(title: str, items) -> str:
            body = "<br>".join(items) if items else "<i>none</i>"
            return f"<p><b>{title}</b><br>{body}</p>"

        parts = [f"<p><b>{result.label}</b> — {result.verdict} "
                 f"(consensus {result.consensus:.2f})</p>"]
        if result.ifp_tanimoto is not None:
            parts += [block("Reference contacts kept", result.conserved),
                      block("Reference contacts lost", result.missing),
                      block("Contacts not in the reference", result.gained)]
        if result.missed_features:
            parts.append(block("Reference feature kinds not fully reproduced",
                               result.missed_features))
        if result.notes:
            parts.append(block("Notes", result.notes))
        self.detail.setHtml("".join(parts))
        self.detail.setVisible(True)

    # ------------------------------------------------------------------------ export
    def _export(self) -> None:
        if not self._results:
            return
        project = self.c.ensure_project()
        default = str(project.root / "exports" / "pose_consistency.csv")
        path, _filter = QFileDialog.getSaveFileName(self, "Export pose consistency", default,
                                                    "CSV (*.csv)")
        if not path:
            return
        from mds.chem.workflow import write_csv
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        write_csv(target, rows(self._results))
        self.status.setText(f"Exported {len(self._results)} row(s) to {target.name}.")
