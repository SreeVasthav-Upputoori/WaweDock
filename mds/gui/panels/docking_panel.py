"""Docking setup: receptor -> ligands -> binding site -> parameters -> run."""

from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGroupBox, QHBoxLayout, QInputDialog, QLabel, QLineEdit,
                               QListWidget, QListWidgetItem, QPushButton, QScrollArea, QSpinBox,
                               QVBoxLayout, QWidget)

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


class DockingPanel(QScrollArea):
    def __init__(self, controller: Controller, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.c = controller
        self.setWidgetResizable(True)
        body = QWidget()
        lay = QVBoxLayout(body)
        lay.addWidget(self._receptor_group())
        lay.addWidget(self._ligand_group())
        lay.addWidget(self._site_group())
        lay.addWidget(self._docking_group())
        lay.addStretch()
        self.setWidget(body)

        c = self.c
        c.structure_loaded.connect(self._on_structure)
        c.receptor_prepared.connect(self._on_receptor)
        c.hetero_found.connect(self._on_hetero)
        c.ligands_changed.connect(self._refresh_ligands)
        c.ligands_prepared.connect(lambda _o: self._refresh_ligands())
        c.box_changed.connect(self._on_box)
        c.reference_changed.connect(self._on_reference)
        c.docking_started.connect(lambda: self._set_running(True))
        c.docking_finished.connect(lambda _o: self._set_running(False))
        c.runner.job_cancelled.connect(lambda _i: self._set_running(c.docking_running))
        c.runner.job_failed.connect(lambda _i, _t: self._set_running(c.docking_running))
        self._refresh_ligands()
        self._site_changed()

    # ------------------------------------------------------------ 1. receptor
    def _receptor_group(self) -> QGroupBox:
        g = QGroupBox("1  Receptor")
        f = QFormLayout(g)
        self.rec_label = QLabel("<i>No structure loaded (File ▸ Open Structure or Fetch from PDB)</i>")
        self.rec_label.setWordWrap(True)
        self.chains = QLineEdit()
        self.chains.setPlaceholderText("all chains  (or e.g. A,B)")
        self.keep_metals = QCheckBox("Keep catalytic/structural metal ions (Zn, Mg, Ca, Mn, Fe, Cu)")
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
        chains = [x.strip() for x in self.chains.text().split(",") if x.strip()] or None
        self.rec_status.setText("Preparing…")
        self.c.prepare_receptor(ReceptorOptions(chains=chains, keep_metals=self.keep_metals.isChecked(),
                                                ph=self.ph.value()))

    def _on_structure(self, name: str, path: str) -> None:
        self.rec_label.setText(f"Structure: <b>{name}</b> ({Path(path).name})")
        self.rec_status.setText("Not prepared yet")

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
        for text, slot in (("Add SMILES…", self._add_smiles), ("Import file…", self._import_file),
                           ("Remove", self._remove_ligand)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        v.addLayout(row)
        self.redock_btn = QPushButton("Redock co-crystal ligand (validation)")
        self.redock_btn.setToolTip("Adds the crystal ligand (re-built in 3D from its chemistry, not its "
                                   "crystal pose), sets the box around it, and reports RMSD to the "
                                   "experimental pose after docking.")
        self.redock_btn.clicked.connect(self._redock)
        v.addWidget(self.redock_btn)
        f = QFormLayout()
        self.ionize = QComboBox()
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
        text, ok = QInputDialog.getText(self, "Add ligand", "SMILES (optionally followed by a name):")
        if not ok or not text.strip():
            return
        parts = text.strip().split(None, 1)
        self.c.add_smiles(parts[0], parts[1] if len(parts) > 1 else "")

    def _import_file(self) -> None:
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
        self.site_method = QComboBox()
        self.site_method.addItems(SITE_METHODS)
        self.site_method.currentIndexChanged.connect(self._site_changed)
        self.site_hetero = QComboBox()
        self.residues = QLineEdit()
        self.residues.setPlaceholderText("e.g. A:25, A:27, B:50   (chain:number)")
        self.padding = _spin(2, 30, 8.0, 0.5)
        self.site_btn = QPushButton("Create Box")
        self.site_btn.clicked.connect(self._make_box)
        self.cx, self.cy, self.cz = (_spin(-9999, 9999, 0, 0.5, 2) for _ in range(3))
        self.sx, self.sy, self.sz = (_spin(1, 126, 20, 1, 1) for _ in range(3))
        for s in (self.cx, self.cy, self.cz, self.sx, self.sy, self.sz):
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
        self.c.run_docking(DockingParameters(self.exh.value(), self.modes.value(), self.erange.value(),
                                             self.cpu.value(), self.seed.value(), self.scoring.currentData()))

    def _set_running(self, running: bool) -> None:
        self.run_btn.setEnabled(not running)
        self.cancel_btn.setEnabled(running)

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
