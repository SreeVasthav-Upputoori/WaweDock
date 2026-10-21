"""Clickable docking-score comparison, using native Qt graphics only."""
from __future__ import annotations

import math
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsRectItem, QGraphicsScene, QGraphicsView


class _ScoreBar(QGraphicsRectItem):
    def __init__(self, rectangle, chart, outcome, pose_index):
        super().__init__(*rectangle)
        self.chart, self.outcome, self.pose_index = chart, outcome, pose_index
        self.setBrush(QBrush(QColor('#408ce6')))
        self.setPen(QPen(Qt.PenStyle.NoPen))
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event):
        self.chart.pose_selected.emit(self.outcome, self.pose_index)
        event.accept()


class ScoreChart(QGraphicsView):
    pose_selected = Signal(object, int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.scene().setBackgroundBrush(QBrush(QColor("#ffffff")))
        self.setMinimumHeight(180)
        self.setAccessibleName('Best docking score per ligand. Click a bar to inspect its pose.')
        self.bars = []

    def set_outcomes(self, outcomes):
        scene = self.scene()
        scene.clear()
        self.bars = []
        points = []
        for outcome in outcomes:
            valid = [(i, p.score) for i, p in enumerate(outcome.poses) if math.isfinite(p.score)]
            if outcome.status == 'ok' and valid:
                i, score = min(valid, key=lambda pair: pair[1])
                points.append((outcome, i, score))
        points.sort(key=lambda item: item[2])
        if not points:
            scene.addText('No scores available. Run docking to compare ligands.')
            return
        scene.addText('Best score per ligand (kcal/mol) - click a bar to inspect')
        low = min(0, min(score for _, _, score in points))
        high = max(0, max(score for _, _, score in points))
        if high == low:
            high = low + 1
        span = high - low
        top, height, left, width = 42, 180, 70, 76
        y = lambda score: top + (high - score) / span * height
        pen = QPen(QColor('#9aa5b1'))
        for tick in range(5):
            score = low + span * tick / 4
            label = scene.addText(f'{score:.1f}')
            label.setPos(4, y(score) - 10)
            scene.addLine(left, y(score), left + len(points) * width, y(score), pen)
        for n, (outcome, i, score) in enumerate(points):
            x = left + n * width + 10
            bar = _ScoreBar((x, min(y(0), y(score)), 48, max(2, abs(y(score) - y(0)))), self, outcome, i)
            bar.setToolTip(f'{outcome.name}\nPose {outcome.poses[i].rank}: {score:.3f} kcal/mol')
            scene.addItem(bar)
            self.bars.append(bar)
            label = scene.addText(outcome.name)
            label.setTextWidth(width)
            label.setPos(x - 8, top + height + 6)
        scene.setSceneRect(scene.itemsBoundingRect().adjusted(0, 0, 12, 12))
