"""Draws the WaweDock logo (PNG sizes + multi-resolution ICO). Run: python tools/make_logo.py"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PIL import Image
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import (QBrush, QColor, QGuiApplication, QImage, QLinearGradient, QPainter,
                           QPainterPath, QPen, QRadialGradient)

OUT = Path(__file__).resolve().parent.parent / "mds" / "assets"


def draw_icon(size: int = 1024) -> QImage:
    img = QImage(size, size, QImage.Format.Format_ARGB32)
    img.fill(Qt.GlobalColor.transparent)
    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    s = size / 512.0
    p.scale(s, s)

    # rounded-square background (deep ocean gradient)
    bg = QLinearGradient(0, 0, 0, 512)
    bg.setColorAt(0, QColor("#12305c"))
    bg.setColorAt(1, QColor("#06142b"))
    path = QPainterPath()
    path.addRoundedRect(QRectF(16, 16, 480, 480), 104, 104)
    p.fillPath(path, QBrush(bg))
    p.setClipPath(path)

    # soft glow
    glow = QRadialGradient(QPointF(256, 230), 260)
    glow.setColorAt(0, QColor(0, 200, 255, 70))
    glow.setColorAt(1, QColor(0, 200, 255, 0))
    p.fillRect(QRectF(0, 0, 512, 512), QBrush(glow))

    # layered waves along the bottom
    for i, (y, amp, alpha, col) in enumerate([(330, 34, 70, "#1fb5e6"), (360, 30, 110, "#19d3c5"),
                                              (392, 26, 170, "#6ff2e0")]):
        w = QPainterPath()
        w.moveTo(0, 512)
        w.lineTo(0, y)
        for x in range(0, 521, 8):
            w.lineTo(x, y + amp * math.sin((x / 512) * 2 * math.pi * 1.5 + i * 0.9))
        w.lineTo(512, 512)
        w.closeSubpath()
        c = QColor(col)
        c.setAlpha(alpha)
        p.fillPath(w, QBrush(c))

    # hexagonal molecule ring
    cx, cy, r = 256, 214, 96
    pts = [QPointF(cx + r * math.cos(math.radians(60 * k - 90)),
                   cy + r * math.sin(math.radians(60 * k - 90))) for k in range(6)]
    ring = QPainterPath(pts[0])
    for pt in pts[1:]:
        ring.lineTo(pt)
    ring.closeSubpath()
    pen = QPen(QColor("#f2fbff"), 17)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(pen)
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPath(ring)
    inner = QPainterPath()  # aromatic-style inner bonds
    for k in (0, 2, 4):
        a, b = pts[k], pts[(k + 1) % 6]
        va = QPointF(cx, cy)
        a2 = QPointF(a.x() + (va.x() - a.x()) * 0.22, a.y() + (va.y() - a.y()) * 0.22)
        b2 = QPointF(b.x() + (va.x() - b.x()) * 0.22, b.y() + (va.y() - b.y()) * 0.22)
        inner.moveTo(a2)
        inner.lineTo(b2)
    p.setPen(QPen(QColor(150, 225, 255, 200), 9, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
    p.drawPath(inner)

    # atoms at vertices
    for k, pt in enumerate(pts):
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor("#ffffff"))
        p.drawEllipse(pt, 15, 15)
        p.setBrush(QColor("#2fd3ff") if k != 3 else QColor("#ffcf4a"))
        p.drawEllipse(pt, 8, 8)

    # ligand "docking" arrow: a wave-shaped stroke entering the ring
    arrow = QPainterPath(QPointF(66, 300))
    arrow.cubicTo(QPointF(120, 250), QPointF(150, 345), QPointF(205, 292))
    arrow.cubicTo(QPointF(235, 262), QPointF(250, 280), QPointF(300, 300))
    p.setPen(QPen(QColor("#ffcf4a"), 15, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPath(arrow)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor("#ffcf4a"))
    p.drawEllipse(QPointF(66, 300), 17, 17)
    p.end()
    return img


def main() -> None:
    app = QGuiApplication.instance() or QGuiApplication(sys.argv)
    OUT.mkdir(parents=True, exist_ok=True)
    big = draw_icon(1024)
    big.save(str(OUT / "wavedock.png"))
    pil = Image.open(OUT / "wavedock.png").convert("RGBA")
    pil.save(OUT / "wavedock.ico", sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])
    print("wrote", OUT / "wavedock.png", OUT / "wavedock.ico")
    _ = app


if __name__ == "__main__":
    main()
