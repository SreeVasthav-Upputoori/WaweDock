"""Render the approved Binding Pocket identity. Run: python tools/make_logo.py."""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PIL import Image
from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontDatabase, QGuiApplication, QImage, QPainter, QPainterPath
from PySide6.QtSvg import QSvgRenderer

OUT = Path(__file__).resolve().parent.parent / "mds" / "assets"
BLUE = "#064DFF"


def draw_icon(size: int = 1024) -> QImage:
    image = QImage(size, size, QImage.Format.Format_ARGB32)
    image.fill(Qt.GlobalColor.transparent)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    renderer = QSvgRenderer(str(OUT / "wavedock.svg"))
    if not renderer.isValid():
        raise ValueError("Invalid logo SVG")
    renderer.render(painter)
    painter.end()
    return image


def draw_splash() -> QImage:
    image = QImage(1000, 720, QImage.Format.Format_ARGB32)
    image.fill(QColor("white"))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
    painter.translate(258, 80)
    painter.scale(1.1, 1.1)
    pocket = QPainterPath()
    pocket.moveTo(0, 0)
    for x, y in [(300, 0), (300, 68), (186, 165), (300, 262), (300, 330), (0, 330)]:
        pocket.lineTo(x, y)
    pocket.closeSubpath()
    diamond = QPainterPath()
    diamond.moveTo(350, 92)
    for x, y in [(423, 165), (350, 238), (277, 165)]:
        diamond.lineTo(x, y)
    diamond.closeSubpath()
    painter.fillPath(pocket, QColor(BLUE))
    painter.fillPath(diamond, QColor(BLUE))
    painter.resetTransform()
    font = QFont("Segoe UI", 62)
    font.setWeight(QFont.Weight.DemiBold)
    painter.setFont(font)
    painter.setPen(QColor("#171717"))
    painter.drawText(QRectF(40, 500, 920, 145), Qt.AlignmentFlag.AlignCenter, "WaweDock")
    painter.end()
    return image


def main() -> None:
    app = QGuiApplication.instance() or QGuiApplication(sys.argv)
    # Offscreen Qt may not discover Windows fonts automatically.
    font_file = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / "segoeui.ttf"
    if font_file.is_file():
        QFontDatabase.addApplicationFont(str(font_file))
    icon = draw_icon()
    if not icon.save(str(OUT / "wavedock.png")):
        raise OSError("Could not save icon")
    with Image.open(OUT / "wavedock.png") as image:
        image.save(OUT / "wavedock.ico", sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])
    if not draw_splash().save(str(OUT / "wavedock-splash.png")):
        raise OSError("Could not save splash")
    print("Wrote Binding Pocket PNG, ICO and splash assets")
    _ = app


if __name__ == "__main__":
    main()
