"""Interaction pharmacophore tab: derive from a pose, show it in 3D, score the other poses."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QFileDialog, QHBoxLayout, QHeaderView,
                               QLabel, QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from mds.chem.pharmacophore import CAVEAT, Pharmacophore, derive

FEATURE_COLUMNS = ["Feature", "Partner residue(s)", "Tolerance (Å)"]
MATCH_COLUMNS = ["Ligand", "Pose", "Recovery", "Matched", "Missed features"]


class _NumItem(QTableWidgetItem):
    """Sorts numerically while displaying formatted text."""

    def __init__(self, value: float | None, text: str) -> None:
        super().__init__(text)
        self._v = float("-inf") if value is None else value

    def __lt__(self, other: QTableWidgetItem) -> bool:  # noqa: D105
        return self._v < getattr(other, "_v", float("-inf"))


class PharmacophorePanel(QWidget):
    """Turns the selected pose's classified interactions into a reusable 3D query.

    Deriving needs the ProLIF classification, so it depends on the same optional Open Babel
    backend the Interactions tab does. Matching other poses afterwards does not.
    """

    display_requested = Signal(object, bool)     # Pharmacophore | None, show
    pose_requested = Signal(object, int)         # LigandOutcome, pose index

    def __init__(self, controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        self._current: tuple[object, int] | None = None
        self._pharmacophore: Pharmacophore | None = None
        self._derived_from: tuple[object, int] | None = None
        self._await_interactions = False
        self._match_rows: list[tuple[object, int]] = []

        self.derive_btn = QPushButton("Derive from selected pose")
        self.derive_btn.setToolTip("Build a pharmacophore from this pose's classified interactions. "
                                   "Runs the ProLIF analysis first if it has not been run yet.")
        self.derive_btn.clicked.connect(self._derive)
        self.score_btn = QPushButton("Score all poses")
        self.score_btn.setToolTip("Report how much of this pharmacophore each docked pose reproduces.")
        self.score_btn.clicked.connect(self._score)
        self.score_btn.setEnabled(False)
        self.show_in_viewer = QCheckBox("Show in 3D")
        self.show_in_viewer.setToolTip("Draw each feature as a translucent sphere the size of its "
                                       "matching tolerance.")
        self.show_in_viewer.setEnabled(False)
        self.show_in_viewer.toggled.connect(self._display_changed)
        self.export_btn = QPushButton("Export JSON")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self._export)

        self.status = QLabel("Select a pose, then derive a pharmacophore from it.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)

        self.features = QTableWidget(0, len(FEATURE_COLUMNS))
        self.features.setHorizontalHeaderLabels(FEATURE_COLUMNS)
        self.features.verticalHeader().hide()
        self.features.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.features.horizontalHeader().setStretchLastSection(True)
        self.features.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)

        self.matches = QTableWidget(0, len(MATCH_COLUMNS))
        self.matches.setHorizontalHeaderLabels(MATCH_COLUMNS)
        self.matches.verticalHeader().hide()
        self.matches.setSortingEnabled(True)
        self.matches.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.matches.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.matches.setToolTip("Select a row to view that pose.")
        self.matches.horizontalHeader().setStretchLastSection(True)
        self.matches.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.matches.itemSelectionChanged.connect(self._match_selected)

        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.features)
        split.addWidget(self.matches)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)

        caveat = QLabel(f"<small>{CAVEAT}</small>")
        caveat.setWordWrap(True)

        bar = QHBoxLayout()
        bar.addWidget(self.derive_btn)
        bar.addWidget(self.score_btn)
        bar.addWidget(self.show_in_viewer)
        bar.addWidget(self.export_btn)
        bar.addStretch(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(bar)
        layout.addWidget(self.status)
        layout.addWidget(split, 1)
        layout.addWidget(caveat)

        controller.interactions_ready.connect(self._interactions_ready)
        controller.pharmacophore_matched.connect(self._matched)

    # ------------------------------------------------------------- pose state
    def set_pose(self, outcome, pose_index: int) -> None:
        self._current = (outcome, pose_index)
        self.derive_btn.setEnabled(True)
        if self._pharmacophore is None:
            self.status.setText(f"Ready to derive from {outcome.name} "
                                f"pose #{outcome.poses[pose_index].rank}.")

    def clear_pose(self) -> None:
        self._current = None
        self._await_interactions = False
        self.derive_btn.setEnabled(False)
        if self._pharmacophore is None:
            self.status.setText("Select a pose, then derive a pharmacophore from it.")

    def clear(self) -> None:
        """Forget the pharmacophore as well as the pose — used when results are replaced."""
        self.clear_pose()
        self._pharmacophore = None
        self._derived_from = None
        self._match_rows = []
        self.features.setRowCount(0)
        self.matches.setRowCount(0)
        for widget in (self.score_btn, self.show_in_viewer, self.export_btn):
            widget.setEnabled(False)
        self.show_in_viewer.setChecked(False)
        self.display_requested.emit(None, False)
        self.status.setText("Select a pose, then derive a pharmacophore from it.")

    # -------------------------------------------------------------- deriving
    def _derive(self) -> None:
        if self._current is None:
            self.status.setText("Select a pose in the results table first.")
            return
        self._await_interactions = True
        self.status.setText("Classifying interactions for this pose (ProLIF)…")
        # interactions_ready arrives for every analysis; _await_interactions marks ours.
        self.c.analyze_interactions(*self._current)

    def _interactions_ready(self, outcome, pose_index: int, result) -> None:
        if not self._await_interactions or self._current != (outcome, pose_index):
            return
        self._await_interactions = False
        pose = outcome.poses[pose_index]
        source = f"{outcome.name} pose #{pose.rank} ({pose.score:.2f} kcal/mol, predicted)"
        try:
            pharmacophore = derive(result.contacts, outcome.mols[pose_index], source=source)
        except ValueError as exc:
            self.status.setText(f"Could not derive a pharmacophore: {exc}")
            return
        self._pharmacophore = pharmacophore
        self._derived_from = (outcome, pose_index)
        self._show_features(pharmacophore)

    def _show_features(self, pharmacophore: Pharmacophore) -> None:
        self.features.setRowCount(len(pharmacophore.features))
        for row, feature in enumerate(pharmacophore.features):
            self.features.setItem(row, 0, QTableWidgetItem(feature.label))
            self.features.setItem(row, 1, QTableWidgetItem(", ".join(feature.residues)))
            self.features.setItem(row, 2, _NumItem(feature.tolerance, f"{feature.tolerance:.1f}"))
        has = bool(pharmacophore.features)
        for widget in (self.score_btn, self.show_in_viewer, self.export_btn):
            widget.setEnabled(has)
        summary = ", ".join(f"{n}×{kind}" for kind, n in sorted(pharmacophore.counts().items()))
        self.status.setText(f"{len(pharmacophore)} feature(s) from {pharmacophore.source}"
                            + (f" — {summary}." if has else
                               " — none of its interactions define a ligand-side feature."))
        if self.show_in_viewer.isChecked():
            self.display_requested.emit(pharmacophore, True)

    def _display_changed(self, on: bool) -> None:
        self.display_requested.emit(self._pharmacophore if on else None, on)

    # -------------------------------------------------------------- matching
    def _score(self) -> None:
        if self._pharmacophore is None:
            return
        poses, self._match_rows = [], []
        for outcome in self.c.outcomes:
            if outcome.status != "ok":
                continue
            for index, mol in enumerate(outcome.mols):
                poses.append(((outcome.name, outcome.poses[index].rank), mol))
                self._match_rows.append((outcome, index))
        self.status.setText(f"Scoring {len(poses)} pose(s) against {len(self._pharmacophore)} feature(s)…")
        self.c.match_pharmacophore(self._pharmacophore, poses)

    def _matched(self, results: list) -> None:
        if len(results) != len(self._match_rows):
            return  # a stale job, or the results changed underneath it
        self.matches.setSortingEnabled(False)
        self.matches.setRowCount(len(results))
        for row, ((name, rank), result) in enumerate(results):
            missed = ", ".join(sorted({f.label for f in result.missed()})) or "—"
            cells = [QTableWidgetItem(str(name)), _NumItem(rank, str(rank)),
                     _NumItem(result.recovery, f"{result.recovery * 100:.0f}%"),
                     _NumItem(result.n_matched, f"{result.n_matched} / {len(result.matches)}"),
                     QTableWidgetItem(missed)]
            cells[0].setData(Qt.ItemDataRole.UserRole, row)
            for column, item in enumerate(cells):
                self.matches.setItem(row, column, item)
        self.matches.setSortingEnabled(True)
        self.matches.sortItems(2, Qt.SortOrder.DescendingOrder)
        best = max((r.recovery for _k, r in results), default=0.0)
        self.status.setText(f"Scored {len(results)} pose(s); best recovery {best * 100:.0f}% of "
                            f"{len(self._pharmacophore)} feature(s). Recovery is geometric, "
                            f"not a predicted affinity.")

    def _match_selected(self) -> None:
        rows = self.matches.selectionModel().selectedRows()
        if not rows:
            return
        index = self.matches.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        if index is None or index >= len(self._match_rows):
            return
        outcome, pose_index = self._match_rows[index]
        self.pose_requested.emit(outcome, pose_index)

    # ---------------------------------------------------------------- export
    def _export(self) -> None:
        if self._pharmacophore is None:
            return
        path, _filter = QFileDialog.getSaveFileName(self, "Export pharmacophore", "pharmacophore.json",
                                                    "JSON (*.json)")
        if not path:
            return
        from pathlib import Path
        Path(path).write_text(self._pharmacophore.to_json(), encoding="utf-8")
        self.status.setText(f"Exported {len(self._pharmacophore)} feature(s) to {Path(path).name}.")
