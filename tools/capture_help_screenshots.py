"""Capture consistent WaweDock screenshots for the bundled help guide."""

from __future__ import annotations

import tempfile
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from mds.chem.box import DockingBox
from mds.gui.main_window import MainWindow
from mds.project import Project


OUT = ROOT / "docs" / "images"
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT / "tests" / "data"


def capture(window: MainWindow, name: str) -> None:
    window.show()
    QApplication.processEvents()
    window.grab().save(str(OUT / name), "PNG")


def main() -> None:
    app = QApplication.instance() or QApplication([])

    with tempfile.TemporaryDirectory(prefix="wavedock-help-") as folder:
        folder_path = Path(folder)
        settings = QSettings(str(folder_path / "settings.ini"), QSettings.Format.IniFormat)
        window = MainWindow(settings, use_viewer=False)
        window.resize(1500, 920)
        window.controller.set_project(Project.create(folder_path / "project"))
        window.controller.load_structure(DATA / "3PTB.pdb")
        QApplication.processEvents()

        window.workspace.select(0)
        window.docking_panel.show_step(0)
        window.dock_docking.show()
        window.dock_project.show()
        window.dock_selection.hide()
        window.dock_properties.hide()
        window.dock_report.hide()
        window.dock_log.hide()
        window.dock_jobs.hide()
        window.dock_health.hide()
        capture(window, "help-main-window.png")

        window.docking_panel.show_step(0)
        window.docking_panel.chains.setText("A")
        window.docking_panel.ph.setValue(7.4)
        capture(window, "help-receptor.png")

        window.docking_panel.show_step(1)
        window.controller.add_smiles("N=C(N)c1ccccc1", "benzamidine")
        QApplication.processEvents()
        capture(window, "help-ligands.png")

        window.docking_panel.show_step(2)
        window.docking_panel.site_method.setCurrentIndex(3)
        window.docking_panel.cx.setValue(9.0)
        window.docking_panel.cy.setValue(20.0)
        window.docking_panel.cz.setValue(20.0)
        window.docking_panel.sx.setValue(20.0)
        window.docking_panel.sy.setValue(20.0)
        window.docking_panel.sz.setValue(20.0)
        window.controller.set_box(DockingBox((9.0, 20.0, 20.0), (20.0, 20.0, 20.0)))
        QApplication.processEvents()
        capture(window, "help-box.png")

        window.workspace.select(2)
        window.results_panel.results_summary if hasattr(window.results_panel, "results_summary") else None
        QApplication.processEvents()
        capture(window, "help-review.png")

        window.close()


if __name__ == "__main__":
    main()
