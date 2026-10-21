"""Collapsible docking workflow backed by the existing controller and jobs."""

from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGroupBox, QHBoxLayout, QInputDialog, QLabel, QLineEdit,
                               QListWidget, QListWidgetItem, QPushButton, QScrollArea, QSpinBox,
                               QSizePolicy, QStyle, QVBoxLayout, QWidget)

from mds.chem.box import DockingBox
from mds.chem.ligand import SUPPORTED_SUFFIXES
from mds.chem.receptor import HeteroResidue, ReceptorOptions
from mds.chem.workflow import DockingParameters
from mds.gui.controller import Controller

SITE_METHODS = ["Co-crystal ligand", "Selected residues", "Whole protein (blind docking)", "Manual"]


def _spin(minimum: float, maximum: float, value: float, step: float = 1.0, decimals: int = 1) -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setRange(minimum, maximum)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(value)
    return s


STEP_TITLES = ("Receptor Preparation", "Ligand Preparation", "Docking Box Configuration",
               "Docking Execution", "Results Review")


class DockingPanel(QScrollArea):
    step_activated = Signal(int)
    results_requested = Signal()

    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        self.setWidgetResizable(True)
        self._running = controller.docking_running
        self._execution_status = ""
        self.current_step = -1
        self.step_buttons: list[QPushButton] = []
        self.step_labels: list[QLabel] = []
        self.step_arrows: list[QLabel] = []
        self.step_pages: list[QGroupBox] = []
        self.step_status: list[str] = [""] * len(STEP_TITLES)
        self.setFrameShape(QScrollArea.Shape.NoFrame)
        body = QWidget()
        lay = QVBoxLayout(body)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)
        for index, page in enumerate((self._receptor_group(), self._ligand_group(),
                                      self._site_group(), self._docking_group(), self._results_group())):
            page.setTitle("")
            page.setFlat(True)
            button = QPushButton()
            button.setObjectName(f"dockingStep{index + 1}")
            header = QHBoxLayout(button)
            header.setContentsMargins(8, 6, 8, 6)
            arrow = QLabel()
            arrow.setFixedWidth(16)
            label = QLabel()
            label.setWordWrap(True)
            for child in (arrow, label):
                child.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
            header.addWidget(arrow)
            header.addWidget(label, 1)
            self.step_labels.append(label)
            self.step_arrows.append(arrow)
            button.setCheckable(True)
            button.setMinimumHeight(self.fontMetrics().lineSpacing() * 3 + 12)
            button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            button.clicked.connect(lambda _on=False, i=index: self.show_step(
                -1 if self.current_step == i else i))
            self.step_buttons.append(button)
            self.step_pages.append(page)
            lay.addWidget(button)
            lay.addWidget(page)
        lay.addStretch()
        self.setWidget(body)
        self.show_step(0)

        c = self.c
        c.project_changed.connect(self._on_project_changed)
        c.structure_loaded.connect(self._on_structure)
        c.structure_edited.connect(lambda _old, path: self._on_structure(Path(path).stem, path))
        c.receptor_prepared.connect(self._on_receptor)
        c.hetero_found.connect(self._on_hetero)
        c.ligands_changed.connect(self._refresh_ligands)
        c.ligands_prepared.connect(lambda _o: self._refresh_ligands())
        c.box_changed.connect(self._on_box)
        c.reference_changed.connect(self._on_reference)
        c.docking_started.connect(lambda: self._set_running(True))
        c.docking_finished.connect(self._on_docking_finished)
        c.runner.job_cancelled.connect(lambda _i: self._on_job_stopped("Cancelled"))
        c.runner.job_failed.connect(lambda _i, _t: self._on_job_stopped("Failed"))
        for signal in (c.receptor_prepared, c.ligands_changed, c.ligands_prepared,
                       c.box_changed, c.docking_finished):
            signal.connect(self._refresh_steps)
        self._refresh_ligands()
        self._site_changed()
        self._refresh_steps()

    def show_step(self, index: int) -> None:
        """Expand one step, or collapse all with -1, retaining every control value."""
        if index < -1 or index >= len(self.step_pages):
            raise ValueError(f"Unknown docking step: {index}")
        self.current_step = index
        for i, (button, page) in enumerate(zip(self.step_buttons, self.step_pages)):
            expanded = i == index
            button.setChecked(expanded)
            self.step_arrows[i].setPixmap(self.style().standardIcon(
                QStyle.StandardPixmap.SP_ArrowDown if expanded else QStyle.StandardPixmap.SP_ArrowRight
            ).pixmap(12, 12))
            page.setVisible(expanded)
        if index >= 0:
            self.ensureWidgetVisible(self.step_buttons[index])
            self.step_activated.emit(index)

    def _results_group(self) -> QGroupBox:
        group = QGroupBox()
        layout = QVBoxLayout(group)
        self.results_summary = QLabel("Run docking to review poses and scores.")
        self.results_summary.setWordWrap(True)
        layout.addWidget(self.results_summary)
        self.review_btn = QPushButton("Review Results")
        self.review_btn.clicked.connect(lambda: self.results_requested.emit())
        layout.addWidget(self.review_btn)
        return group

    def _refresh_steps(self, *_args) -> None:
        c = self.c
        ready_names = {lig.name for lig in c.prepared}
        entry_names = {lig.name for lig in c.entries}
        ligands_ready = bool(entry_names) and entry_names <= ready_names
        valid_box = c.box is not None and not c.validate_box().errors
        successes = [o for o in c.outcomes if o.status == "ok" and o.poses]
        poses = sum(len(o.poses) for o in successes)
        failures = sum(o.status != "ok" or not o.poses for o in c.outcomes)
        statuses = [
            "Prepared" if c.receptor else ("Not prepared" if c.source_pdb else "Load a structure"),
            (f"{len(ready_names)}/{len(entry_names)} prepared" if entry_names else "Add ligands"),
            "Ready" if valid_box else ("Needs attention" if c.box else "Define a box"),
            "Running" if self._running else (self._execution_status or (
                "Needs attention" if failures else ("Completed" if successes else "Not run"))),
            f"{poses} poses available" if poses else "No poses yet",
        ]
        completed = [bool(c.receptor), ligands_ready, valid_box,
                     bool(successes) and not failures and not self._running and not self._execution_status, False]
        self.step_status = statuses
        for i, button in enumerate(self.step_buttons):
            marker = "\u2713 " if completed[i] else ""
            self.step_labels[i].setText(f"{i + 1}  {STEP_TITLES[i]}\n{marker}{statuses[i]}")
            button.setToolTip(f"{STEP_TITLES[i]}: {statuses[i]}. Click to expand or collapse.")
            button.setAccessibleName(f"Step {i + 1}: {STEP_TITLES[i]}. {statuses[i]}")
        self.results_summary.setText(
            f"{len(successes)} ligand(s), {poses} pose(s) available. "
            f"{failures} ligand(s) failed or produced no poses."
            if c.outcomes else "Run docking to review poses and scores.")
        self.review_btn.setEnabled(bool(c.outcomes))
        missing = []
        if c.receptor is None:
            missing.append("prepare a receptor")
        if not c.prepared:
            missing.append("prepare ligands")
        if not valid_box:
            missing.append("define a valid docking box")
        self.run_readiness.setText("Before docking: " + "; ".join(missing) + "."
                                   if missing else "Ready to dock with the current project inputs.")

    def _on_project_changed(self) -> None:
        self._execution_status = ""
        if self.c.source_pdb is None:
            self.rec_label.setText("No structure loaded (File > Open Structure or Fetch from PDB)")
            self.rec_status.clear()
        self._on_hetero(self.c.hetero)
        self._refresh_ligands()
        self._on_reference(self.c.reference_name)
        self._on_box(self.c.box)
        self._set_running(self.c.docking_running)
        self._refresh_steps()

    # ------------------------------------------------------------ 1. receptor
    def _receptor_group(self) -> QGroupBox:
        g = QGroupBox("1  Receptor")
        f = QFormLayout(g)
        f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.rec_label = QLabel("<i>No structure loaded (File ▸ Open Structure or Fetch from PDB)</i>")
        self.rec_label.setWordWrap(True)
        self.chains = QLineEdit()
        self.chains.setPlaceholderText("all chains  (or e.g. A,B)")
        self.keep_metals = QCheckBox("Keep metal ions")
        self.keep_metals.setToolTip("Retain Zn, Mg, Ca, Mn, Fe and Cu ions during receptor preparation.")
        self.keep_metals.setChecked(True)
        self.ph = _spin(4.0, 10.0, 7.4, 0.1)
        self.prep_rec_btn = QPushButton("Prepare Receptor")
        self.prep_rec_btn.setToolTip("Remove waters, rebuild missing atoms, assign protonation at the "
                                     "chosen pH (PDB2PQR/PROPKA), build PDBQT (Meeko)")
        self.prep_rec_btn.clicked.connect(self._prepare_receptor)
        self.rec_status = QLabel("")
        self.rec_status.setWordWrap(True)
        f.addRow(self.rec_label)
        f.addRow("Chains", self.chains)
        f.addRow(self.keep_metals)
        f.addRow("Protonation pH", self.ph)
        f.addRow(self.prep_rec_btn)
        f.addRow(self.rec_status)
        return g

    def _prepare_receptor(self) -> None:
        self.show_step(0)
        chains = [x.strip() for x in self.chains.text().split(",") if x.strip()] or None
        self.rec_status.setText("Preparing…")
        self.c.prepare_receptor(ReceptorOptions(chains=chains, keep_metals=self.keep_metals.isChecked(),
                                                ph=self.ph.value()))

    def _on_structure(self, name: str, path: str) -> None:
        self.rec_label.setText(f"Structure: <b>{name}</b> ({Path(path).name})")
        self.rec_status.setText("Not prepared yet")
        self._refresh_steps()

    def _on_receptor(self, result) -> None:
        w = len(result.report.warnings)
        self.rec_status.setText(f"✔ Prepared: {result.n_atoms} atoms, {w} warning(s). See the "
                                "Preparation Report tab.")

    # ------------------------------------------------------------- 2. ligands
    def _ligand_group(self) -> QGroupBox:
        g = QGroupBox("2  Ligands")
        v = QVBoxLayout(g)
        self.lig_list = QListWidget()
        self.lig_list.setMinimumHeight(110)
        v.addWidget(self.lig_list)
        row = QHBoxLayout()
        for text, slot in (("SMILES…", self._add_smiles), ("Import…", self._import_file),
                           ("Remove", self._remove_ligand)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        v.addLayout(row)
        self.redock_btn = QPushButton("Redock co-crystal ligand")
        self.redock_btn.setToolTip("Adds the crystal ligand (re-built in 3D from its chemistry, not its "
                                   "crystal pose), sets the box around it, and reports RMSD to the "
                                   "experimental pose after docking.")
        self.redock_btn.clicked.connect(self._redock)
        v.addWidget(self.redock_btn)
        f = QFormLayout()
        f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.ionize = QComboBox()
        self.ionize.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.ionize.setMinimumContentsLength(12)
        self.ionize.addItem("Keep input protonation", "none")
        self.ionize.addItem("Simple pH 7.4 rules (acids, amines, amidines)", "simple")
        self.confs = QSpinBox()
        self.confs.setRange(1, 20)
        self.lig_seed = QSpinBox()
        self.lig_seed.setRange(0, 2_000_000_000)
        self.lig_seed.setValue(42)
        f.addRow("Ionisation", self.ionize)
        f.addRow("Input conformers", self.confs)
        f.addRow("Random seed", self.lig_seed)
        v.addLayout(f)
        self.prep_lig_btn = QPushButton("Prepare Ligands")
        self.prep_lig_btn.clicked.connect(self._prepare_ligands)
        v.addWidget(self.prep_lig_btn)
        return g

    def _add_smiles(self) -> None:
        self.show_step(1)
        text, ok = QInputDialog.getText(self, "Add ligand", "SMILES (optionally followed by a name):")
        if not ok or not text.strip():
            return
        parts = text.strip().split(None, 1)
        self.c.add_smiles(parts[0], parts[1] if len(parts) > 1 else "")

    def _import_file(self) -> None:
        self.show_step(1)
        exts = " ".join(f"*{s}" for s in SUPPORTED_SUFFIXES)
        path, _ = QFileDialog.getOpenFileName(self, "Import ligands", "", f"Ligand files ({exts})")
        if path:
            self.c.import_ligand_file(Path(path))

    def _remove_ligand(self) -> None:
        for item in self.lig_list.selectedItems():
            self.c.remove_ligand(item.data(Qt.ItemDataRole.UserRole))

    def _redock(self) -> None:
        candidates = [h for h in self.c.hetero if h.kind == "ligand" and h.extracted_to]
        if not candidates:
            self.c.error.emit("No co-crystal ligand", "Prepare the receptor first; extracted hetero "
                              "groups are listed after preparation.")
            return
        if len(candidates) == 1:
            self.c.redock_cocrystal(candidates[0])
            return
        labels = [h.label for h in candidates]
        label, ok = QInputDialog.getItem(self, "Redock co-crystal ligand", "Ligand:", labels, 0, False)
        if ok:
            self.c.redock_cocrystal(candidates[labels.index(label)])

    def _prepare_ligands(self) -> None:
        self.show_step(1)
        self.c.prepare_ligands(ionize=self.ionize.currentData(), n_conformers=self.confs.value(),
                               seed=self.lig_seed.value())

    def _refresh_ligands(self) -> None:
        self.lig_list.clear()
        ready = {p.name for p in self.c.prepared}
        for e in self.c.entries:
            mark = "✔ " if e.name in ready else "○ "
            item = QListWidgetItem(f"{mark}{e.name}")
            item.setData(Qt.ItemDataRole.UserRole, e.name)
            self.lig_list.addItem(item)

    # ----------------------------------------------------------- 3. binding site
    def _site_group(self) -> QGroupBox:
        g = QGroupBox("3  Binding site / docking box")
        f = QFormLayout(g)
        f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.site_method = QComboBox()
        self.site_method.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.site_method.setMinimumContentsLength(12)
        self.site_method.addItems(SITE_METHODS)
        self.site_method.currentIndexChanged.connect(self._site_changed)
        self.site_hetero = QComboBox()
        self.site_hetero.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.site_hetero.setMinimumContentsLength(12)
        self.residues = QLineEdit()
        self.residues.setPlaceholderText("e.g. A:25, A:27, B:50   (chain:number)")
        self.padding = _spin(2, 30, 8.0, 0.5)
        self.site_btn = QPushButton("Create Box")
        self.site_btn.clicked.connect(self._make_box)
        self.cx, self.cy, self.cz = (_spin(-9999, 9999, 0, 0.5, 2) for _ in range(3))
        self.sx, self.sy, self.sz = (_spin(1, 126, 20, 1, 1) for _ in range(3))
        for s in (self.cx, self.cy, self.cz, self.sx, self.sy, self.sz):
            s.setMinimumWidth(60)
            s.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            s.valueChanged.connect(self._manual_edit)
        self._updating_box = False
        self.box_status = QLabel("")
        self.box_status.setWordWrap(True)
        f.addRow("Method", self.site_method)
        f.addRow("Ligand", self.site_hetero)
        f.addRow("Residues", self.residues)
        f.addRow("Padding (Å)", self.padding)
        f.addRow(self.site_btn)
        for label, widgets in (("Centre x/y/z", (self.cx, self.cy, self.cz)),
                               ("Size x/y/z (Å)", (self.sx, self.sy, self.sz))):
            row = QHBoxLayout()
            for w in widgets:
                row.addWidget(w)
            f.addRow(label, row)
        f.addRow(self.box_status)
        return g

    def _site_changed(self) -> None:
        m = self.site_method.currentIndex()
        self.site_hetero.setEnabled(m == 0)
        self.residues.setEnabled(m == 1)
        self.padding.setEnabled(m in (0, 1))
        self.site_btn.setEnabled(m != 3)

    def _on_hetero(self, hetero: list[HeteroResidue]) -> None:
        self.site_hetero.clear()
        for h in hetero:
            if h.kind == "ligand" and h.extracted_to:
                self.site_hetero.addItem(f"{h.label} ({h.atoms} atoms)", h)

    def _make_box(self) -> None:
        m = self.site_method.currentIndex()
        try:
            if m == 0:
                h = self.site_hetero.currentData()
                if h is None:
                    raise ValueError("No co-crystal ligand available (prepare the receptor first).")
                self.c.box_from_hetero(h, self.padding.value())
            elif m == 1:
                residues = set()
                for tok in self.residues.text().replace(";", ",").split(","):
                    tok = tok.strip()
                    if tok:
                        chain, num = tok.split(":")
                        residues.add((chain.strip(), int(num)))
                if not residues:
                    raise ValueError("Enter residues as chain:number, e.g. A:25, A:27")
                self.c.box_from_residues(residues, self.padding.value())
            elif m == 2:
                self.c.blind_box()
        except (ValueError, OSError) as exc:
            self.c.error.emit("Cannot create box", str(exc))

    def _manual_edit(self) -> None:
        if self._updating_box:
            return
        self.c.set_box(DockingBox((self.cx.value(), self.cy.value(), self.cz.value()),
                                  (self.sx.value(), self.sy.value(), self.sz.value())))

    def _on_box(self, box: DockingBox | None) -> None:
        if box is None:
            self.box_status.setText("")
            return
        self._updating_box = True
        for spin, val in zip((self.cx, self.cy, self.cz, self.sx, self.sy, self.sz), (*box.center, *box.size)):
            spin.setValue(val)
        self._updating_box = False
        rep = self.c.validate_box()
        icon = {"error": "✖", "warning": "⚠", "info": "ℹ"}
        self.box_status.setText("<br>".join(f"{icon[e.level.value]} {e.text}" for e in rep.entries))

    # ----------------------------------------------------------- 4. docking
    def _docking_group(self) -> QGroupBox:
        g = QGroupBox("4  AutoDock Vina")
        f = QFormLayout(g)
        f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.scoring = QComboBox()
        self.scoring.addItem("Vina", "vina")
        self.scoring.addItem("Vinardo", "vinardo")
        self.exh = QSpinBox()
        self.exh.setRange(1, 1024)
        self.exh.setValue(16)
        self.exh.setToolTip("Search effort. 8 = Vina default; use 32+ for large boxes or flexible ligands.")
        self.modes = QSpinBox()
        self.modes.setRange(1, 100)
        self.modes.setValue(9)
        self.erange = _spin(0.5, 10.0, 3.0, 0.5)
        self.cpu = QSpinBox()
        self.cpu.setRange(0, os.cpu_count() or 64)
        self.cpu.setSpecialValueText("all cores")
        self.seed = QSpinBox()
        self.seed.setRange(0, 2_000_000_000)
        self.seed.setValue(42)
        f.addRow("Scoring function", self.scoring)
        f.addRow("Exhaustiveness", self.exh)
        f.addRow("Number of poses", self.modes)
        f.addRow("Energy range (kcal/mol)", self.erange)
        f.addRow("CPU threads", self.cpu)
        f.addRow("Random seed", self.seed)
        self.run_readiness = QLabel()
        self.run_readiness.setWordWrap(True)
        f.addRow(self.run_readiness)
        self.run_btn = QPushButton("▶  Run Docking")
        self.run_btn.setStyleSheet("font-weight: 600; padding: 6px;")
        self.run_btn.clicked.connect(self._run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self.c.cancel_docking)
        self.ref_label = QLabel("")
        self.ref_label.setWordWrap(True)
        row = QHBoxLayout()
        row.addWidget(self.run_btn, 3)
        row.addWidget(self.cancel_btn, 1)
        f.addRow(row)
        f.addRow(self.ref_label)
        note = QLabel("<small>Scores are Vina predictions in kcal/mol, not measured affinities.</small>")
        note.setWordWrap(True)
        f.addRow(note)
        return g

    def _run(self) -> None:
        self.show_step(3)
        self.c.run_docking(DockingParameters(self.exh.value(), self.modes.value(), self.erange.value(),
                                             self.cpu.value(), self.seed.value(), self.scoring.currentData()))

    def _set_running(self, running: bool) -> None:
        self._running = running
        if running:
            self._execution_status = ""
        self.run_btn.setEnabled(not running)
        self.cancel_btn.setEnabled(running)
        self._refresh_steps()

    def _on_docking_finished(self, _outcomes) -> None:
        self._execution_status = ""
        self._set_running(False)

    def _on_job_stopped(self, status: str) -> None:
        if self._running and not self.c.docking_running:
            self._execution_status = status
        self._set_running(self.c.docking_running)

    def _on_reference(self, name: str) -> None:
        self.ref_label.setText(f"Validation reference: <b>{name}</b> (RMSD to crystal pose will be "
                               "reported)" if name else "")

    def apply_params(self, p: DockingParameters) -> None:
        self.exh.setValue(p.exhaustiveness)
        self.modes.setValue(p.num_modes)
        self.erange.setValue(p.energy_range)
        self.cpu.setValue(p.cpu)
        self.seed.setValue(p.seed)
        self.scoring.setCurrentIndex(max(0, self.scoring.findData(p.scoring)))
