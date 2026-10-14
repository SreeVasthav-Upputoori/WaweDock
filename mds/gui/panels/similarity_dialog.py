"""Similarity & clustering dialog: Tanimoto/Morgan fingerprints over the prepared ligands."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QAbstractItemView, QDialog, QDoubleSpinBox, QFormLayout,
                               QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton,
                               QTableWidget, QTableWidgetItem, QVBoxLayout)

from mds.chem.similarity import SimilarityResult, most_similar_pairs
from mds.gui.controller import Controller

CAVEAT = ("Tanimoto similarity on Morgan (ECFP4-like) fingerprints is a measure of structural "
         "resemblance. It does not predict biological activity, binding, or potency.")

# Distinguishable, colour-blind-friendlier palette for cluster highlighting; cycles for >8 clusters.
_CLUSTER_COLORS = [
    "#cfe8ff", "#ffe0b3", "#d8f5d0", "#f5d0e8", "#e8e0ff", "#fff3b0", "#c9f0ea", "#f0d0c9",
]


class SimilarityDialog(QDialog):
    def __init__(self, controller: Controller, parent=None) -> None:
        super().__init__(parent)
        self.c = controller
        self.setWindowTitle("Similarity & Clustering")
        self.resize(760, 560)

        self.cutoff = QDoubleSpinBox()
        self.cutoff.setRange(0.05, 0.95)
        self.cutoff.setSingleStep(0.05)
        self.cutoff.setValue(0.4)
        self.cutoff.setToolTip("Butina clustering distance threshold (1 - Tanimoto similarity). "
                               "Smaller = tighter, more numerous clusters.")
        self.reference = QLineEdit()
        self.reference.setPlaceholderText("Optional: SMILES of a reference compound to rank by similarity")
        self.compute_btn = QPushButton("Compute")
        self.compute_btn.clicked.connect(self._compute)

        form = QFormLayout()
        form.addRow("Cluster cutoff (distance)", self.cutoff)
        form.addRow("Reference compound", self.reference)
        bar = QHBoxLayout()
        bar.addLayout(form, 1)
        bar.addWidget(self.compute_btn)

        self.status = QLabel("Prepare at least two ligands, then click Compute.")
        self.status.setWordWrap(True)

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Ligand", "Cluster", "Similarity to reference"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)

        self.pairs_label = QLabel("<b>Most similar pairs</b>")
        self.pairs = QTableWidget(0, 3)
        self.pairs.setHorizontalHeaderLabels(["Ligand A", "Ligand B", "Tanimoto similarity"])
        self.pairs.verticalHeader().hide()
        self.pairs.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.pairs.horizontalHeader().setStretchLastSection(True)
        self.pairs.setMaximumHeight(160)

        caveat = QLabel(f"<small>{CAVEAT}</small>")
        caveat.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.addLayout(bar)
        layout.addWidget(self.status)
        layout.addWidget(self.table, 2)
        layout.addWidget(self.pairs_label)
        layout.addWidget(self.pairs, 1)
        layout.addWidget(caveat)

        controller.similarity_ready.connect(self._show_result)

    def _compute(self) -> None:
        self.status.setText("Computing…")
        self.compute_btn.setEnabled(False)
        self.c.compute_similarity(self.cutoff.value(), self.reference.text())

    def _show_result(self, result: SimilarityResult) -> None:
        self.compute_btn.setEnabled(True)
        n_clusters = len(set(result.clusters))
        text = f"{len(result.names)} ligand(s) in {n_clusters} cluster(s) at cutoff {self.cutoff.value():.2f}."
        if result.reference_name:
            text += f" Ranked by similarity to '{result.reference_name}'."
        self.status.setText(text)

        order = range(len(result.names))
        if result.reference_similarity is not None:
            order = sorted(order, key=lambda i: result.reference_similarity[i], reverse=True)

        self.table.setRowCount(len(result.names))
        for row, i in enumerate(order):
            name_item = QTableWidgetItem(result.names[i])
            cluster_item = QTableWidgetItem(str(result.clusters[i]))
            sim_text = f"{result.reference_similarity[i]:.3f}" if result.reference_similarity else ""
            sim_item = QTableWidgetItem(sim_text)
            color = QColor(_CLUSTER_COLORS[result.clusters[i] % len(_CLUSTER_COLORS)])
            for item in (name_item, cluster_item, sim_item):
                item.setBackground(color)
            self.table.setItem(row, 0, name_item)
            self.table.setItem(row, 1, cluster_item)
            self.table.setItem(row, 2, sim_item)

        pairs = most_similar_pairs(result, top_n=10)
        self.pairs.setRowCount(len(pairs))
        for row, (a, b, s) in enumerate(pairs):
            self.pairs.setItem(row, 0, QTableWidgetItem(a))
            self.pairs.setItem(row, 1, QTableWidgetItem(b))
            self.pairs.setItem(row, 2, QTableWidgetItem(f"{s:.3f}"))
