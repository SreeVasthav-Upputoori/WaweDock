"""WaweDock main window.

Layout
------
* Centre ............ 3D viewer (3Dmol.js via QtWebEngine)
* Left dock ......... Compact Project Explorer; Selection & Cleanup on demand
* Right docks ....... Docking Setup; Properties on demand
* Bottom docks ...... Results | Preparation Report | Log | Jobs (hidden until requested)
* Menus/toolbars .... built from a shared ``ActionRegistry``

Every dock and toolbar has a stable ``objectName``; Qt needs it to save and restore layouts.
"""

from __future__ import annotations

import logging
import os
import platform
from pathlib import Path

from PySide6 import __version__ as pyside_version
from PySide6.QtCore import QSettings, QStandardPaths, Qt, QUrl, qVersion
from PySide6.QtGui import (QActionGroup, QCloseEvent, QDesktopServices, QDragEnterEvent, QDropEvent,
                           QIcon, QKeySequence, QPixmap)
from PySide6.QtWidgets import (QApplication, QDockWidget, QFileDialog, QInputDialog, QLabel,
                               QMainWindow, QMenu, QMessageBox, QStyle, QTableWidget,
                               QTableWidgetItem, QTabWidget, QToolBar, QVBoxLayout, QDialog, QWidget)
from rdkit import Chem

from mds import APP_NAME, __version__
from mds.chem.workflow import LigandOutcome
from mds.core.paths import assets_dir
from mds.gui.actions import ActionRegistry
from mds.gui.controller import Controller
from mds.gui.layouts import LayoutManager
from mds.gui.panels.docking_panel import DockingPanel
from mds.gui.panels.jobs_panel import JobsPanel
from mds.gui.panels.log_panel import LogPanel
from mds.gui.panels.project_explorer import ProjectExplorer
from mds.gui.panels.properties_panel import PropertiesPanel
from mds.gui.panels.selection_panel import SelectionPanel
from mds.gui.panels.results_panel import ReportPanel, ResultsPanel
from mds.gui.panels.similarity_dialog import SimilarityDialog
from mds.gui.status_bar import StudioStatusBar
from mds.gui.theme import Theme, apply_theme
from mds.gui.viewport import ViewportHost
from mds.jobs import JobRunner
from mds.jobs.diagnostics import JobSystemSelfTest
from mds.project import Project

log = logging.getLogger("mds.gui")

STRUCTURE_SUFFIXES = frozenset({".pdb", ".ent", ".cif", ".mmcif"})
LIGAND_SUFFIXES = frozenset({".sdf", ".mol", ".mol2", ".smi", ".smiles", ".csv", ".txt"})

SCIENTIFIC_NOTICE = (
    "<p>Docking scores produced by this software are <b>predictions</b> from an empirical "
    "scoring function (AutoDock Vina, reported in kcal/mol). They are not experimentally "
    "measured binding affinities and are not evidence of biological activity.</p>"
    "<p>Protein–ligand contacts shown are <b>calculated, distance-based candidates</b> and "
    "must be distinguished from experimentally verified interactions.</p>"
    "<p>Protonation states, tautomers, side-chain rebuilding and the chosen binding site are "
    "assumptions that affect results; validate a target by redocking its co-crystal ligand.</p>"
)


def app_icon() -> QIcon:
    return QIcon(str(assets_dir() / "wavedock.png"))


class MainWindow(QMainWindow):
    def __init__(self, settings: QSettings, log_dir: Path | None = None, *, use_viewer: bool = True) -> None:
        super().__init__()
        self.settings = settings
        self.log_dir = log_dir
        self.setWindowTitle(APP_NAME)
        self.setWindowIcon(app_icon())
        self.setObjectName("mainWindow")
        self.resize(1500, 920)
        self.setAcceptDrops(True)
        self.setDockOptions(QMainWindow.DockOption.AllowTabbedDocks | QMainWindow.DockOption.AllowNestedDocks
                            | QMainWindow.DockOption.AnimatedDocks)

        self.job_runner = JobRunner(parent=self)
        self.controller = Controller(self.job_runner, self)
        self.actions_ = ActionRegistry(self)
        self.viewer = None
        self._pose_object = ""
        self._receptor_objects: list[str] = []
        self._source_object = ""
        self._focus_state = None
        self._similarity_dialog: SimilarityDialog | None = None

        self._build_central(use_viewer)
        self._build_status_bar()
        self._build_docks()
        self._build_actions()
        self._build_menus()
        self._build_toolbars()
        self._connect_signals()

        self.layouts = LayoutManager(self, settings)
        self.layouts.capture_default()
        self.layouts.restore_session()
        self.set_theme(Theme(settings.value("appearance/theme", Theme.LIGHT.value)))
        log.info("%s %s started", APP_NAME, __version__)

    # ================================================================ build
    def _build_central(self, use_viewer: bool) -> None:
        self.viewport = ViewportHost()
        self.setCentralWidget(self.viewport)
        if not use_viewer:
            return
        try:
            from mds.render.mol3d import Mol3DViewer
            self.viewer = Mol3DViewer()
            self.viewport.set_viewer(self.viewer)
        except Exception as exc:  # noqa: BLE001 — QtWebEngine missing/blocked: run without 3D
            log.error("3D viewer unavailable: %s", exc)

    def _build_status_bar(self) -> None:
        self.status = StudioStatusBar()
        self.setStatusBar(self.status)

    def _add_dock(self, widget: QWidget, title: str, name: str, area: Qt.DockWidgetArea) -> QDockWidget:
        dock = QDockWidget(title, self)
        dock.setObjectName(name)
        dock.setWidget(widget)
        self.addDockWidget(area, dock)
        return dock

    def _build_docks(self) -> None:
        left, right, bottom = (Qt.DockWidgetArea.LeftDockWidgetArea, Qt.DockWidgetArea.RightDockWidgetArea,
                               Qt.DockWidgetArea.BottomDockWidgetArea)
        self.project_explorer = ProjectExplorer()
        self.selection_panel = SelectionPanel(self.controller)
        self.properties = PropertiesPanel()
        self.docking_panel = DockingPanel(self.controller)
        self.results_panel = ResultsPanel(self.controller)
        self.report_panel = ReportPanel(self.controller)
        self.log_panel = LogPanel()
        self.jobs_panel = JobsPanel(self.job_runner)

        self.dock_project = self._add_dock(self.project_explorer, "Project Explorer", "dockProject", left)
        self.dock_selection = self._add_dock(self.selection_panel, "Selection & Cleanup", "dockSelection", left)
        self.dock_docking = self._add_dock(self.docking_panel, "Docking Setup", "dockDocking", right)
        self.dock_properties = self._add_dock(self.properties, "Properties", "dockProperties", right)
        self.dock_results = self._add_dock(self.results_panel, "Results", "dockResults", bottom)
        self.dock_report = self._add_dock(self.report_panel, "Preparation Report", "dockReport", bottom)
        self.dock_log = self._add_dock(self.log_panel, "Log", "dockLog", bottom)
        self.dock_jobs = self._add_dock(self.jobs_panel, "Jobs", "dockJobs", bottom)

        self.tabifyDockWidget(self.dock_docking, self.dock_properties)
        for a, b in ((self.dock_results, self.dock_report), (self.dock_report, self.dock_log),
                     (self.dock_log, self.dock_jobs)):
            self.tabifyDockWidget(a, b)
        self.dock_docking.raise_()
        self.dock_results.raise_()
        self.tabifyDockWidget(self.dock_project, self.dock_selection)
        self.dock_project.raise_()
        for dock in (self.dock_selection, self.dock_properties, self.dock_results,
                     self.dock_report, self.dock_log, self.dock_jobs):
            dock.hide()
        self.resizeDocks([self.dock_project, self.dock_docking], [220, 340], Qt.Orientation.Horizontal)
        self.resizeDocks([self.dock_results], [260], Qt.Orientation.Vertical)
        logging.getLogger("mds").addHandler(self.log_panel.handler)

    def _icon(self, pixmap: QStyle.StandardPixmap) -> QIcon:
        return self.style().standardIcon(pixmap)

    def _build_actions(self) -> None:
        a, sp, c = self.actions_, QStyle.StandardPixmap, self.controller
        sk = QKeySequence.StandardKey
        # File
        a.add("file.new", "&New Project…", slot=self._new_project, shortcut=sk.New, icon=self._icon(sp.SP_FileIcon))
        a.add("file.open", "&Open Project…", slot=self._open_project, shortcut=sk.Open,
              icon=self._icon(sp.SP_DialogOpenButton))
        a.add("file.save", "&Save Project", slot=self._save_project, shortcut=sk.Save,
              icon=self._icon(sp.SP_DialogSaveButton))
        a.add("file.save_as", "Save Project &As (copy)…", slot=self._save_project_as, shortcut=sk.SaveAs)
        a.add("file.archive", "Archive Project (zip)…", slot=self._archive_project)
        a.add("file.open_structure", "Open &Structure…", slot=self._open_structure, shortcut="Ctrl+Shift+O",
              icon=self._icon(sp.SP_DirOpenIcon), tip="Open a PDB or mmCIF protein structure")
        a.add("file.fetch_pdb", "Fetch from &PDB…", slot=self._fetch_pdb, shortcut="Ctrl+Shift+P",
              icon=self._icon(sp.SP_ArrowDown), tip="Download a structure from the RCSB Protein Data Bank")
        a.add("file.import_ligands", "&Import Ligands…", slot=lambda: self.docking_panel._import_file(), shortcut="Ctrl+I")
        a.add("file.export_image", "Export &Image…", slot=self._export_image)
        a.add("file.exit", "E&xit", slot=self.close, shortcut=sk.Quit)
        # Edit
        a.add("edit.undo", "Undo Structure Edit", slot=c.undo_structure_edit, shortcut=sk.Undo)
        a.add("edit.redo", "Redo Structure Edit", slot=c.redo_structure_edit, shortcut=sk.Redo)
        a.add("edit.delete", "Delete Selected Atoms", slot=c.delete_selected_residues)
        a["edit.undo"].setEnabled(False)
        a["edit.redo"].setEnabled(False)
        a.add("edit.select_around", "Show Residues Around Ligand…", slot=self._select_around)
        a.add("edit.clear_selection", "Clear Residue Display", slot=self._clear_around)
        # View
        a.add("view.fullscreen", "&Full Screen", slot=self._toggle_fullscreen, shortcut=sk.FullScreen, checkable=True)
        a.add("view.theme.light", "&Light", slot=lambda on: on and self.set_theme(Theme.LIGHT), checkable=True)
        a.add("view.theme.dark", "&Dark", slot=lambda on: on and self.set_theme(Theme.DARK), checkable=True)
        a.add("view.theme.toggle", "Toggle Theme", slot=self._toggle_theme, shortcut="Ctrl+Shift+T")
        a.add("view.reset_camera", "Reset &Camera", slot=lambda: self._viewer_call("zoom_to", None), shortcut="Home")
        a.add("view.label_residues", "Label Protein Residues", slot=lambda on: self._viewer_call("label_residues", on),
              checkable=True)
        a.add("view.measure.distance", "Measure &Distance", slot=lambda: self._measure("distance"))
        a.add("view.measure.angle", "Measure &Angle", slot=lambda: self._measure("angle"))
        a.add("view.measure.dihedral", "Measure D&ihedral", slot=lambda: self._measure("dihedral"))
        a.add("view.measure.off", "Stop Measuring", slot=lambda: self._measure("none"))
        a.add("view.measure.clear", "Clear Measurements", slot=lambda: self._viewer_call("clear_measurements"))
        a.add("view.exit_fullscreen", "Exit Molecular Full Screen",
              slot=lambda: a["view.fullscreen"].setChecked(False), shortcut="Escape")
        self.addAction(a["view.exit_fullscreen"])
        # Build
        a.add("build.from_smiles", "Ligand from &SMILES…", slot=lambda: self.docking_panel._add_smiles())
        # Prepare
        a.add("prepare.protein", "Prepare &Protein", slot=lambda: self.docking_panel._prepare_receptor(),
              tip="Receptor preparation workflow")
        a.add("prepare.ligands", "Prepare &Ligands", slot=lambda: self.docking_panel._prepare_ligands(),
              tip="Ligand preparation workflow")
        # Docking
        a.add("docking.box", "Define Docking &Box", slot=self._focus_site)
        a.add("docking.run", "&Run Docking (AutoDock Vina)", slot=lambda: self.docking_panel._run(), shortcut="F5")
        a.add("docking.cancel", "&Cancel Docking", slot=c.cancel_docking)
        a.add("docking.screen", "&Virtual Screening queue / resume…", planned_phase=9,
              tip="Importing several ligands already docks them as a batch")
        # Analyse
        a.add("analyse.results", "Docking &Results", slot=lambda: (self.dock_results.show(), self.dock_results.raise_()))
        a.add("analyse.descriptors", "Ligand &Descriptors / Drug-likeness…", slot=self._show_descriptors)
        a.add("analyse.similarity", "&Similarity && Clustering…", slot=self._show_similarity,
              tip="Tanimoto/Morgan fingerprint similarity and clustering of the prepared ligands")
        a.add("analyse.interaction_diagram", "2D Interaction &Diagram", slot=self._show_interactions,
              tip="Protein-ligand interaction fingerprint (ProLIF) for the selected pose")
        # Tools
        a.add("tools.selftest", "Run Job System &Self-Test", slot=self._run_self_test)
        a.add("tools.cancel_all", "&Cancel All Jobs", slot=self.job_runner.cancel_all)
        a.add("tools.log_folder", "Open &Log Folder", slot=self._open_log_folder)
        a.add("tools.project_folder", "Open &Project Folder", slot=self._open_project_folder)
        # Window
        a.add("window.save_layout", "&Save Layout As…", slot=self._save_named_layout)
        a.add("window.reset_layout", "&Reset Layout", slot=lambda: self.layouts.reset())
        # Help
        a.add("help.notice", "&Scientific Interpretation Notice", slot=self._show_scientific_notice)
        a.add("help.about", "&About WaweDock", slot=self._show_about)

    def _build_menus(self) -> None:
        a, mb = self.actions_, self.menuBar()

        def menu(title: str, keys: list[str | None]) -> QMenu:
            m = mb.addMenu(title)
            for k in keys:
                m.addSeparator() if k is None else m.addAction(a[k])
            return m

        menu("&File", ["file.new", "file.open", "file.save", "file.save_as", "file.archive", None,
                       "file.open_structure", "file.fetch_pdb", "file.import_ligands", "file.export_image", None,
                       "file.exit"])
        menu("&Edit", ["edit.undo", "edit.redo", "edit.delete", None,
                       "edit.select_around", "edit.clear_selection"])
        m = mb.addMenu("&View")
        theme_menu = m.addMenu("&Theme")
        group = QActionGroup(self)
        for key in ("view.theme.light", "view.theme.dark"):
            group.addAction(a[key])
            theme_menu.addAction(a[key])
        theme_menu.addSeparator()
        theme_menu.addAction(a["view.theme.toggle"])
        self.toolbars_menu = m.addMenu("Tool&bars")
        m.addSeparator()
        m.addAction(a["view.reset_camera"])
        m.addAction(a["view.label_residues"])
        meas = m.addMenu("&Measure")
        for k in ("view.measure.distance", "view.measure.angle", "view.measure.dihedral", "view.measure.off",
                  "view.measure.clear"):
            meas.addAction(a[k])
        m.addSeparator()
        m.addAction(a["view.fullscreen"])
        menu("&Build", ["build.from_smiles"])
        menu("&Prepare", ["prepare.protein", "prepare.ligands"])
        menu("&Docking", ["docking.box", "docking.run", "docking.cancel", None, "docking.screen"])
        menu("&Analyse", ["analyse.results", "analyse.descriptors", "analyse.interaction_diagram",
                          "analyse.similarity"])
        m = menu("&Tools", ["tools.project_folder", "tools.log_folder", None, "tools.cancel_all"])
        m.addMenu("&Diagnostics").addAction(a["tools.selftest"])
        m = mb.addMenu("&Window")
        panels = m.addMenu("&Panels")
        for dock in (self.dock_project, self.dock_docking, self.dock_properties,
                     self.dock_selection, self.dock_results, self.dock_report, self.dock_log,
                     self.dock_jobs):
            panels.addAction(dock.toggleViewAction())
        m.addSeparator()
        m.addAction(a["window.save_layout"])
        self.restore_layout_menu = m.addMenu("Restore &Layout")
        self.delete_layout_menu = m.addMenu("&Delete Layout")
        self.restore_layout_menu.aboutToShow.connect(
            lambda: self._fill_layout_menu(self.restore_layout_menu, self.layouts.restore_named))
        self.delete_layout_menu.aboutToShow.connect(
            lambda: self._fill_layout_menu(self.delete_layout_menu, self.layouts.delete_named))
        m.addAction(a["window.reset_layout"])
        menu("&Help", ["help.notice", None, "help.about"])

    def _build_toolbars(self) -> None:
        a, sp = self.actions_, QStyle.StandardPixmap
        icons = {"docking.box": sp.SP_FileDialogDetailedView,
                 "docking.run": sp.SP_MediaPlay, "docking.cancel": sp.SP_MediaStop,
                 "view.reset_camera": sp.SP_BrowserReload,
                 "view.fullscreen": sp.SP_TitleBarMaxButton}
        for key, icon in icons.items():
            a[key].setIcon(self._icon(icon))
        tb = QToolBar("Workspace", self)
        tb.setObjectName("toolbarWorkspace")
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        for keys in (("file.open", "file.save", "file.open_structure", "file.fetch_pdb"),
                     ("docking.box", "docking.run", "docking.cancel"),
                     ("view.reset_camera", "view.fullscreen")):
            if tb.actions():
                tb.addSeparator()
            for key in keys:
                tb.addAction(a[key])
        self.addToolBar(tb)
        self.toolbars_menu.addAction(tb.toggleViewAction())

    def _connect_signals(self) -> None:
        c, ex = self.controller, self.project_explorer
        self.docking_panel.step_activated.connect(
            lambda _i: (self.dock_docking.show(), self.dock_docking.raise_()))
        self.docking_panel.results_requested.connect(self.actions_["analyse.results"].trigger)
        ex.item_selected.connect(self._on_item_selected)
        ex.visibility_changed.connect(lambda n, v: self._viewer_call("set_visible", n, v))
        ex.representation_requested.connect(lambda n, r: self._viewer_call("set_representation", n, r))
        ex.color_requested.connect(lambda n, s: self._viewer_call("set_color_scheme", n, s))
        ex.zoom_requested.connect(lambda n: self._viewer_call("zoom_to", n))
        ex.remove_requested.connect(self._remove_object)
        self.job_runner.active_count_changed.connect(self.status.set_active_jobs)
        self.job_runner.job_progress.connect(lambda _i, pct, _m: self.status.set_progress(pct))
        c.message.connect(lambda m: self.status.showMessage(m, 8000))
        c.error.connect(self._show_error)
        c.project_changed.connect(self._update_title)
        c.structure_loaded.connect(self._on_structure_loaded)
        c.structure_edited.connect(self._on_structure_edited)
        c.receptor_prepared.connect(self._on_receptor_prepared)
        c.ligands_prepared.connect(self._on_ligands_prepared)
        c.box_changed.connect(self._on_box_changed)
        c.docking_finished.connect(self._on_docking_finished)
        c.selection_expand_requested.connect(lambda: self._viewer_call("expand_selection_to_chain"))
        c.selection_invert_requested.connect(lambda: self._viewer_call("invert_selection"))
        c.selection_mode_requested.connect(lambda mode: self._viewer_call("set_selection_mode", mode))
        c.selection_radius_requested.connect(lambda radius: self._viewer_call("select_within_radius", radius))
        c.edit_history_changed.connect(self._on_edit_history)
        c.editing_view_requested.connect(self._show_editing_structure)
        c.selection_cleared.connect(lambda: self._viewer_call("clear_selection"))
        self.results_panel.pose_selected.connect(self._show_pose)
        if self.viewer is not None:
            self.viewer.visibility_changed.connect(ex.set_visible)
            self.viewer.atom_info.connect(self.status.set_selection)
            self.viewer.measurement.connect(lambda m: self.status.showMessage(m, 15000))
            self.viewer.selection_changed.connect(c.set_selection)
            self.viewer.edit_requested.connect(self._on_viewer_edit)

    def _on_edit_history(self, undo: bool, redo: bool) -> None:
        self.actions_["edit.undo"].setEnabled(undo)
        self.actions_["edit.redo"].setEnabled(redo)

    def _on_viewer_edit(self, action: str) -> None:
        handlers = {"delete": self.controller.delete_selected_residues,
                    "undo": self.controller.undo_structure_edit,
                    "redo": self.controller.redo_structure_edit}
        if action in handlers:
            handlers[action]()

    def _show_editing_structure(self) -> None:
        if self.controller.source_pdb is None:
            return
        for name in self._receptor_objects:
            self.project_explorer.set_visible(name, False)
            self._viewer_call("set_visible", name, False)
        source = self.controller.source_pdb.stem
        self.project_explorer.set_visible(source, True)
        self._viewer_call("set_visible", source, True)
        self._viewer_call("set_measure_mode", "none")

    # ============================================================ viewer glue
    def _viewer_call(self, method: str, *args) -> None:
        if self.viewer is None:
            self.status.showMessage("3D viewer is unavailable in this session", 5000)
            return
        getattr(self.viewer, method)(*args)

    def _object_props(self, **kw) -> dict:
        return {k: v for k, v in kw.items()}

    def _on_structure_loaded(self, name: str, path: str, *, preserve_view: bool = False) -> None:
        import gemmi
        from mds.chem.edit import atom_records
        if self._source_object and self._source_object != name:
            self._viewer_call("remove_object", self._source_object)
            self.project_explorer.remove_object(self._source_object)
        if not preserve_view:
            for overlay in self._receptor_objects:
                self._viewer_call("remove_object", overlay)
                self.project_explorer.remove_object(overlay)
            self._receptor_objects.clear()
        self._source_object = name
        st = gemmi.read_structure(path)
        # One canonical representation for PDB and mmCIF. Keep every model, chain,
        # hydrogen, alternate location, water and hetero group on import.
        fmt = "cif" if Path(path).suffix.lower() in (".cif", ".mmcif") else "pdb"
        text = st.make_mmcif_document().as_string() if fmt == "cif" else st.make_pdb_string()
        if self.viewer:
            self.viewer.add_object(name, text, fmt, "protein", atom_records(st), preserve_view)
        h = self.controller.hetero
        waters = sum(1 for x in h if x.kind == "water")
        ligs = [x.label for x in h if x.kind == "ligand"]
        self.project_explorer.upsert_object("Proteins", name, {
            "Type": "Protein structure (as imported)", "File": Path(path).name,
            "Waters": waters, "Hetero groups": ", ".join(ligs) or "none"}, kind="protein")
        self._viewer_call("set_background", self.current_theme() is Theme.DARK)
        self.dock_docking.raise_()

    def _on_structure_edited(self, name: str, path: str) -> None:
        """Reload the loaded protein after a delete-selection / keep-only edit."""
        if self.viewer:
            self.viewer.remove_object(name)
        self.project_explorer.remove_object(name)
        # Drop any 'prepared' object we had for the pre-edit structure: it is now stale.
        for prepared_name in self._receptor_objects:
            if self.viewer:
                self.viewer.remove_object(prepared_name)
            self.project_explorer.remove_object(prepared_name)
        self._receptor_objects.clear()
        self._viewer_call("clear_measurements")
        self._viewer_call("clear_residue_selection")
        self._on_structure_loaded(Path(path).stem, path, preserve_view=True)
        self.dock_selection.show()
        self.dock_selection.raise_()

    def _on_receptor_prepared(self, result) -> None:
        raw = self.controller.source_pdb.stem if self.controller.source_pdb else "receptor"
        name = f"{raw} (prepared)"
        self._receptor_objects.append(name)
        if self.viewer:
            self.viewer.add_object(name, result.prepared_pdb.read_text(), "pdb", "protein")
        self.project_explorer.upsert_object("Proteins", name, {
            "Type": "Prepared receptor (hydrogens added)", "Atoms": result.n_atoms,
            "PDBQT": result.pdbqt.name, "Warnings": len(result.report.warnings)}, kind="protein")
        if self.project_explorer.find(raw):
            self.project_explorer.set_visible(raw, False)
            self._viewer_call("set_visible", raw, False)
        for h in result.hetero:
            if h.kind == "ligand" and h.extracted_to and self.viewer:
                label = f"{h.label} (crystal)"
                self._receptor_objects.append(label)
                self.viewer.add_object(label, h.extracted_to.read_text(), "pdb", "ligand")
                self.project_explorer.upsert_object("Ligands", label, {
                    "Type": "Co-crystal ligand (extracted)", "Atoms": h.atoms}, kind="ligand")
        if self.controller.box:
            self._on_box_changed(self.controller.box)

    def _on_ligands_prepared(self, outcome) -> None:
        for lig in outcome.prepared:
            props = {"Type": "Prepared ligand", "SMILES": lig.smiles, **lig.properties}
            self.project_explorer.upsert_object("Ligands", lig.name, props, visible=False, kind="ligand")
            if self.viewer:
                self.viewer.add_object(lig.name, Chem.MolToMolBlock(lig.mol), "sdf", "ligand")
                self.viewer.set_visible(lig.name, False)

    def _on_box_changed(self, box) -> None:
        if box is None:
            self._viewer_call("clear_box")
            self.project_explorer.remove_object("Docking box")
            return
        self._viewer_call("show_box", box.center, box.size)
        self.project_explorer.upsert_object("Docking Boxes", "Docking box", {
            "Type": "AutoDock Vina search space", **{k: f"{v:.2f}" for k, v in box.as_dict().items()},
            "Volume (Å³)": f"{box.volume:,.0f}"}, kind="box")

    def _on_docking_finished(self, outcomes: list[LigandOutcome]) -> None:
        self.project_explorer.clear_category("Docking Poses")
        for o in outcomes:
            self.project_explorer.upsert_object("Docking Poses", f"{o.name}: {o.status}", {
                "Type": "Docking outcome", "Status": o.status, "Message": o.message,
                "Best score (kcal/mol)": o.best_score if o.best_score is not None else "-",
                "Poses": len(o.poses)}, visible=False, kind="result")
        self.dock_results.show()
        self.dock_results.raise_()

    def _show_pose(self, outcome: LigandOutcome, pose_index: int) -> None:
        if self.viewer is None:
            return
        if self._pose_object:
            self.viewer.remove_object(self._pose_object)
        name = f"pose: {outcome.name} #{outcome.poses[pose_index].rank}"
        self._pose_object = name
        self.viewer.add_object(name, Chem.MolToMolBlock(outcome.mols[pose_index]), "sdf", "ligand")
        rec = next((n for n in self.project_explorer.object_names("Proteins") if n.endswith("(prepared)")), None)
        if rec:
            self.viewer.show_residues_around(name, 4.5, rec)
        self.status.showMessage(f"{name}: score {outcome.poses[pose_index].score:.2f} kcal/mol (predicted)", 10000)

    def _remove_object(self, name: str) -> None:
        self._viewer_call("remove_object", name)
        self.project_explorer.remove_object(name)

    # ================================================================ slots
    def current_theme(self) -> Theme:
        return Theme(self.settings.value("appearance/theme", Theme.LIGHT.value))

    def set_theme(self, theme: Theme) -> None:
        app = QApplication.instance()
        if app is not None:
            apply_theme(app, theme)
        self.settings.setValue("appearance/theme", theme.value)
        self.actions_[f"view.theme.{theme.value}"].setChecked(True)
        self._viewer_call("set_background", theme is Theme.DARK)

    def _toggle_theme(self) -> None:
        self.set_theme(Theme.LIGHT if self.current_theme() is Theme.DARK else Theme.DARK)

    def _toggle_fullscreen(self, on: bool) -> None:
        if on and self._focus_state is None:
            self._focus_state = (self.saveState(), self.isMaximized())
            for widget in self.findChildren(QDockWidget) + self.findChildren(QToolBar):
                widget.hide()
            self.menuBar().hide()
            self.statusBar().hide()
            self.showFullScreen()
        elif not on and self._focus_state is not None:
            state, maximized = self._focus_state
            self._focus_state = None
            self.showMaximized() if maximized else self.showNormal()
            self.restoreState(state)
            self.menuBar().show()
            self.statusBar().show()

    def _on_item_selected(self, title: str, props: dict) -> None:
        self.properties.show_properties(title, props)
        if props:
            self.dock_properties.raise_()
        self.status.set_selection(title)

    def _show_error(self, title: str, text: str) -> None:
        log.error("%s: %s", title, text)
        # Non-blocking (window-modal, returns immediately) so background jobs keep reporting.
        box = QMessageBox(QMessageBox.Icon.Warning, title, text, QMessageBox.StandardButton.Ok, self)
        box.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        box.open()

    def _update_title(self) -> None:
        p = self.controller.project
        self.setWindowTitle(f"{p.name} — {APP_NAME}" if p else APP_NAME)

    # ---- file actions
    def _start_dir(self) -> str:
        return self.settings.value("paths/last_dir", QStandardPaths.writableLocation(
            QStandardPaths.StandardLocation.DocumentsLocation))

    def _remember_dir(self, path: str) -> None:
        self.settings.setValue("paths/last_dir", str(Path(path).parent))

    def _open_structure(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open structure", self._start_dir(),
                                              "Protein structures (*.pdb *.ent *.cif *.mmcif)")
        if path:
            self._remember_dir(path)
            self.controller.load_structure(Path(path))

    def _fetch_pdb(self) -> None:
        pdb_id, ok = QInputDialog.getText(self, "Fetch from RCSB PDB", "4-character PDB ID (e.g. 1HSG):")
        if ok and pdb_id.strip():
            self.controller.fetch_pdb(pdb_id)

    def _new_project(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Choose an empty folder for the new project", self._start_dir())
        if not path:
            return
        root = Path(path)
        if (root / "project.json").exists():
            QMessageBox.warning(self, "New Project", "That folder already contains a project; use Open Project.")
            return
        self._reset_workspace()
        self.controller.new_project(root)
        self.status.showMessage(f"Project created: {root}", 6000)

    def _open_project(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Open project folder", self._start_dir())
        if not path:
            return
        try:
            project = Project.open(Path(path))
        except (FileNotFoundError, ValueError) as exc:
            QMessageBox.warning(self, "Open Project", str(exc))
            return
        self._reset_workspace()
        self.controller.set_project(project)
        self.controller.load_project_state()
        self.docking_panel.apply_params(self.controller.params)
        self.status.showMessage(f"Opened {project.name}", 6000)

    def _reset_workspace(self) -> None:
        self._receptor_objects.clear()
        self._source_object = ""
        for cat in ("Proteins", "Ligands", "Docking Boxes", "Docking Poses"):
            for name in self.project_explorer.object_names(cat):
                self._viewer_call("remove_object", name)
            self.project_explorer.clear_category(cat)
        self._viewer_call("clear_box")
        if self._pose_object:
            self._viewer_call("remove_object", self._pose_object)
            self._pose_object = ""
        self.results_panel.show_outcomes([])
        self.controller.reset_state()

    def _save_project(self) -> None:
        self.controller.ensure_project()
        self.controller.save_project()
        self.status.showMessage(f"Project saved: {self.controller.project.root}", 5000)

    def _save_project_as(self) -> None:
        if self.controller.project is None:
            return
        path = QFileDialog.getExistingDirectory(self, "Choose an empty folder for the copy", self._start_dir())
        if path:
            try:
                copy = self.controller.project.duplicate(Path(path))
            except FileExistsError as exc:
                QMessageBox.warning(self, "Save As", str(exc))
                return
            self.controller.save_project()
            self.controller.set_project(copy)

    def _archive_project(self) -> None:
        if self.controller.project is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Archive project", f"{self.controller.project.name}.zip",
                                              "Zip archive (*.zip)")
        if path:
            self.controller.save_project()
            self.status.showMessage(f"Archive written: {self.controller.project.archive(Path(path))}", 8000)

    def _export_image(self) -> None:
        if self.viewer is None:
            return
        project = self.controller.ensure_project()
        path, _ = QFileDialog.getSaveFileName(self, "Export image (3x resolution)",
                                              str(project.root / "exports" / "view.png"), "PNG image (*.png)")
        if path:
            self.viewer.export_image(path, 3)
            self.status.showMessage(f"Image written to {path}", 6000)

    # ---- viewer actions
    def _measure(self, mode: str) -> None:
        self._viewer_call("set_measure_mode", mode)

    def _select_around(self) -> None:
        ligs = self.project_explorer.object_names("Ligands")
        prots = [n for n in self.project_explorer.object_names("Proteins")]
        if not ligs or not prots:
            QMessageBox.information(self, "Residues around ligand", "Load a protein and a ligand first.")
            return
        lig, ok = QInputDialog.getItem(self, "Residues around ligand", "Ligand:", ligs, 0, False)
        if not ok:
            return
        prot = next((p for p in prots if p.endswith("(prepared)")), prots[0])
        cutoff, ok = QInputDialog.getDouble(self, "Residues around ligand", "Distance cutoff (Å):", 4.5, 1.0, 12.0, 1)
        if ok:
            self.project_explorer.set_visible(lig, True)
            self._viewer_call("set_visible", lig, True)
            self._viewer_call("show_residues_around", lig, cutoff, prot)

    def _clear_around(self) -> None:
        self._viewer_call("clear_residue_selection")

    def _focus_site(self) -> None:
        self.docking_panel.show_step(2)
        self.dock_docking.show()
        self.dock_docking.raise_()
        self.docking_panel.site_method.setFocus()

    # ---- analysis
    def _show_descriptors(self) -> None:
        prepared = self.controller.prepared
        if not prepared:
            QMessageBox.information(self, "Ligand descriptors", "Prepare ligands first.")
            return
        cols = ["Name", "MW", "LogP", "TPSA", "HBD", "HBA", "RotB", "Ro5_violations", "Lipinski_pass",
                "Veber_pass", "FormalCharge", "Torsions(Vina)"]
        dlg = QDialog(self)
        dlg.setWindowTitle("Ligand descriptors and drug-likeness filters")
        table = QTableWidget(len(prepared), len(cols))
        table.setHorizontalHeaderLabels(cols)
        for r, lig in enumerate(prepared):
            for cidx, col in enumerate(cols):
                value = lig.name if col == "Name" else lig.properties.get(col, "")
                table.setItem(r, cidx, QTableWidgetItem(str(value)))
        table.resizeColumnsToContents()
        lay = QVBoxLayout(dlg)
        lay.addWidget(table)
        lay.addWidget(QLabel("<small>Lipinski/Veber are empirical rules of thumb for oral drug-likeness, "
                             "not predictions of activity.</small>"))
        dlg.resize(980, 340)
        dlg.exec()

    def _show_interactions(self) -> None:
        if not self.controller.outcomes:
            QMessageBox.information(self, "Interaction diagram", "Run a docking calculation first, "
                                    "then select a pose in the Results panel.")
            return
        self.dock_results.show()
        self.dock_results.raise_()
        tabs = self.results_panel.findChild(QTabWidget)
        if tabs is not None:
            tabs.setCurrentWidget(self.results_panel.interactions_panel)
        self.results_panel.interactions_panel._analyse()

    def _show_similarity(self) -> None:
        if len(self.controller.prepared) < 2:
            QMessageBox.information(self, "Similarity & Clustering",
                                    "Prepare at least two ligands first.")
            return
        if self._similarity_dialog is None:
            self._similarity_dialog = SimilarityDialog(self.controller, self)
        self._similarity_dialog.show()
        self._similarity_dialog.raise_()
        self._similarity_dialog.activateWindow()

    # ---- misc
    def _run_self_test(self) -> None:
        self.job_runner.submit(JobSystemSelfTest())
        self.dock_jobs.show()
        self.dock_jobs.raise_()

    def _open_log_folder(self) -> None:
        if self.log_dir is not None:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.log_dir)))

    def _open_project_folder(self) -> None:
        if self.controller.project:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.controller.project.root)))

    def _save_named_layout(self) -> None:
        name, ok = QInputDialog.getText(self, "Save Layout", "Layout name:")
        name = name.strip().replace("/", "-")
        if ok and name:
            self.layouts.save_named(name)

    def _fill_layout_menu(self, menu: QMenu, handler) -> None:
        menu.clear()
        names = self.layouts.names()
        if not names:
            menu.addAction("(no saved layouts)").setEnabled(False)
        for name in names:
            menu.addAction(name, lambda n=name: handler(n))

    def _show_scientific_notice(self) -> None:
        QMessageBox.information(self, "Scientific Interpretation Notice", SCIENTIFIC_NOTICE)

    def _show_about(self) -> None:
        box = QMessageBox(self)
        box.setWindowTitle(f"About {APP_NAME}")
        box.setIconPixmap(QPixmap(str(assets_dir() / "wavedock.png")).scaledToWidth(96, Qt.TransformationMode.SmoothTransformation))
        box.setText(f"<h3>{APP_NAME} {__version__}</h3><p>Open-source molecular docking workbench: structure "
                    "preparation, AutoDock Vina docking, 3D visualisation and result analysis.</p>"
                    f"<p>Python {platform.python_version()} · Qt {qVersion()} · PySide6 {pyside_version}</p>"
                    "<p><small>AutoDock Vina (Apache-2.0), RDKit (BSD), Meeko (LGPL), Gemmi (MPL-2.0), "
                    "PDB2PQR/PROPKA, ProLIF (Apache-2.0), 3Dmol.js (BSD-3).<br>"
                    "Interaction analysis additionally needs Open Babel (GPL-2.0), installed "
                    "separately — it is not included in this program.</small></p>" + SCIENTIFIC_NOTICE)
        box.exec()

    # ======================================================= drag & drop
    @staticmethod
    def _paths(event: QDragEnterEvent | QDropEvent) -> list[Path]:
        if not event.mimeData().hasUrls():
            return []
        paths = [Path(u.toLocalFile()) for u in event.mimeData().urls() if u.isLocalFile()]
        return [p for p in paths if p.suffix.lower() in STRUCTURE_SUFFIXES | LIGAND_SUFFIXES]

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        event.acceptProposedAction() if self._paths(event) else event.ignore()

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        for p in self._paths(event):
            if p.suffix.lower() in STRUCTURE_SUFFIXES:
                self.controller.load_structure(p)
            else:
                self.controller.import_ligand_file(p)
        event.acceptProposedAction()

    # ============================================================ close
    def closeEvent(self, event: QCloseEvent) -> None:  # noqa: N802
        active = self.job_runner.active_count()
        if active:
            answer = QMessageBox.question(self, "Jobs Running",
                                          f"{active} background job(s) are still running. Cancel them and exit?")
            if answer != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.job_runner.cancel_all()
            self.job_runner.wait_for_done(8000)
        self.actions_["view.fullscreen"].setChecked(False)
        self.controller.save_project()
        self.layouts.save_session()
        self.settings.sync()
        logging.getLogger("mds").removeHandler(self.log_panel.handler)
        log.info("%s closed", APP_NAME)
        event.accept()


def default_log_dir() -> Path:
    base = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppLocalDataLocation)
    return Path(base or os.getcwd()) / "logs"
