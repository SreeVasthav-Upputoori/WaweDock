"""Application controller: holds the docking workflow state and drives background jobs.

Panels never run science themselves; they call methods here and listen to the signals.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QObject, QStandardPaths, Signal
from rdkit import Chem

from mds.chem import box as boxlib
from mds.chem.box import DockingBox
from mds.chem.contacts import Contact, cached_protein_atoms, pose_contacts
from mds.chem.interactions import InteractionResult
from mds.chem.ligand import LigandEntry, PreparedLigand, entry_from_smiles, read_ligands
from mds.chem.poses import fetch_ccd_template, reference_from_pdb_ligand
from mds.chem.receptor import HeteroResidue, ReceptorOptions, ReceptorResult, classify_hetero, read_structure
from mds.chem.report import PreparationError, Report
from mds.chem.workflow import DockingParameters, LigandOutcome, results_rows, write_csv
from mds.jobs.latest import LatestJobQueue
from mds.jobs import Job, JobRunner
from mds.jobs.science_jobs import (DockingJob, FetchPdbJob, InteractionJob, LigandPrepOutcome,
                                   PrepareLigandsJob, PrepareReceptorJob, InspectReceptorJob, ContactsJob, SimilarityJob)
from mds.project import Project

log = logging.getLogger("mds.controller")


class _ReferenceJob(Job):
    """Download the chemical definition of a co-crystal ligand and build the reference pose."""

    def __init__(self, hetero: HeteroResidue, cache: Path) -> None:
        self.hetero, self.cache = hetero, cache
        self.name = f"Reference ligand {hetero.name}"

    def run(self, ctx) -> tuple[Chem.Mol, Chem.Mol]:
        ctx.progress(20, "Fetching ligand definition from RCSB")
        template = fetch_ccd_template(self.hetero.name, self.cache)
        ref = reference_from_pdb_ligand(self.hetero.extracted_to, template)
        return template, ref


class Controller(QObject):
    project_changed = Signal()
    structure_loaded = Signal(str, str)         # object name, path
    structure_edited = Signal(str, str)         # object name, path (after delete/keep-only)
    hetero_found = Signal(list)                 # list[HeteroResidue]
    selection_changed = Signal(dict, str)       # detail, summary
    selection_cleared = Signal()
    selection_expand_requested = Signal()       # main window forwards to viewer
    selection_invert_requested = Signal()
    selection_mode_requested = Signal(str)
    selection_radius_requested = Signal(float)
    edit_history_changed = Signal(bool, bool)
    editing_view_requested = Signal()
    receptor_checked = Signal(object)           # read-only Report
    receptor_check_busy = Signal(bool)
    receptor_prepared = Signal(object)          # ReceptorResult
    ligands_changed = Signal()
    ligands_prepared = Signal(object)           # LigandPrepOutcome
    box_changed = Signal(object)                # DockingBox | None
    reference_changed = Signal(str)             # name or ""
    docking_started = Signal()
    docking_finished = Signal(list)             # list[LigandOutcome]
    contacts_ready = Signal(object, int, object)
    contacts_failed = Signal(object, int, str)
    interactions_started = Signal()
    interactions_ready = Signal(object, int, object)   # LigandOutcome, pose index, InteractionResult
    similarity_ready = Signal(object)            # SimilarityResult
    report = Signal(str, str)                   # title, text
    message = Signal(str)
    error = Signal(str, str)                    # title, text

    def __init__(self, runner: JobRunner, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.runner = runner
        self._contact_queue = LatestJobQueue(runner, self)
        self._box_geometry = boxlib.BoxGeometryCache()
        self.project: Project | None = None
        self.source_pdb: Path | None = None
        self.hetero: list[HeteroResidue] = []
        self.receptor: ReceptorResult | None = None
        self.entries: list[LigandEntry] = []
        self.prepared: list[PreparedLigand] = []
        self.ligand_options = {"ionize": "none", "n_conformers": 1, "seed": 42}
        self._box_ligand_path: Path | None = None
        self.box: DockingBox | None = None
        self.reference: Chem.Mol | None = None
        self.reference_name = ""
        self.params = DockingParameters(exhaustiveness=16)
        self.outcomes: list[LigandOutcome] = []
        self.run_dir: Path | None = None
        self._callbacks: dict[int, Callable[[object], None]] = {}
        self._docking_job: int | None = None
        self._receptor_check_job: int | None = None
        self._interaction_job: int | None = None
        self._cached_protein_rdkit = None             # (receptor_pdb path, RDKit Mol)
        self.selection: dict[str, list[dict]] = {}    # object_name -> [{"chain", "resi"}]
        self._undo_structures: list[Path] = []
        self._redo_structures: list[Path] = []
        runner.job_finished.connect(self._job_finished)
        runner.job_failed.connect(self._job_failed)
        runner.job_cancelled.connect(self._job_cancelled)

    # ------------------------------------------------------------------ jobs
    def _submit(self, job: Job, on_done: Callable[[object], None]) -> int:
        job_id = self.runner.submit(job)
        self._callbacks[job_id] = on_done
        return job_id

    def _job_finished(self, job_id: int, result: object) -> None:
        cb = self._callbacks.pop(job_id, None)
        if cb:
            try:
                cb(result)
            except Exception as exc:  # noqa: BLE001
                log.exception("Post-processing of job %d failed", job_id)
                self.error.emit("Error", f"{type(exc).__name__}: {exc}")

    def _finish_receptor_check(self, job_id: int) -> None:
        if job_id == self._receptor_check_job:
            self._receptor_check_job = None
            self.receptor_check_busy.emit(False)

    def _job_failed(self, job_id: int, tb: str) -> None:
        self._finish_receptor_check(job_id)
        if self._callbacks.pop(job_id, None) is None:
            return
        last = tb.strip().splitlines()[-1] if tb.strip() else "unknown error"
        if job_id == self._docking_job:
            self._docking_job = None
        if job_id == self._interaction_job:
            self._interaction_job = None
        self.error.emit(self.runner.name(job_id) + " failed", last)

    def _job_cancelled(self, job_id: int) -> None:
        self._finish_receptor_check(job_id)
        self._callbacks.pop(job_id, None)
        if job_id == self._docking_job:
            self._docking_job = None
            self.message.emit("Docking cancelled")
        if job_id == self._interaction_job:
            self._interaction_job = None

    # --------------------------------------------------------------- project
    def ensure_project(self) -> Project:
        if self.project is None:
            docs = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.DocumentsLocation)
            from datetime import datetime
            root = Path(docs) / "WaweDock Projects" / f"Project-{datetime.now():%Y%m%d-%H%M%S}"
            self.set_project(Project.create(root))
        return self.project  # type: ignore[return-value]

    def set_project(self, project: Project) -> None:
        self.project = project
        self.project_changed.emit()

    def new_project(self, root: Path) -> None:
        self.reset_state()
        self.set_project(Project.create(root))

    def reset_state(self) -> None:
        self._contact_queue.clear()
        self._box_geometry.clear()
        self._box_ligand_path = None
        self.source_pdb, self.receptor, self.box = None, None, None
        self.hetero, self.entries, self.prepared, self.outcomes = [], [], [], []
        self.reference, self.reference_name, self.run_dir = None, "", None
        self._cached_protein_rdkit = None
        self._reset_edit_history()

    def _reset_edit_history(self) -> None:
        self._undo_structures.clear()
        self._redo_structures.clear()
        self.clear_selection()
        self.edit_history_changed.emit(False, False)

    # ------------------------------------------------------------- structure
    def load_structure(self, path: Path) -> None:
        project = self.ensure_project()
        path = Path(path)
        try:
            st = read_structure(path)
        except PreparationError as exc:
            self.error.emit("Cannot open structure", str(exc))
            return
        copy = project.import_original(path)
        self._contact_queue.clear()
        self._box_geometry.clear()
        self._box_ligand_path = None
        self.source_pdb = copy
        self.receptor = None
        self._cached_protein_rdkit = None
        self._reset_edit_history()
        self.hetero = classify_hetero(st)
        self.hetero_found.emit(self.hetero)
        self.structure_loaded.emit(path.stem, str(copy))
        self.message.emit(f"Loaded {path.name}")

    def fetch_pdb(self, pdb_id: str) -> None:
        project = self.ensure_project()
        self._submit(FetchPdbJob(pdb_id, project.root / "originals"),
                     lambda p: self.load_structure(p))  # type: ignore[arg-type]

    # -------------------------------------------------------------- receptor
    def inspect_receptor(self, options: ReceptorOptions) -> None:
        if self.source_pdb is None:
            self.error.emit("No structure", "Open a protein structure first.")
            return
        if self._receptor_check_job is not None:
            return
        source, project = self.source_pdb, self.project
        def done(result):
            self._receptor_check_job = None
            self.receptor_check_busy.emit(False)
            if self.source_pdb != source or self.project is not project:
                return
            self.report.emit("Receptor check", result.render())
            self.receptor_checked.emit(result)
        self._receptor_check_job = self._submit(InspectReceptorJob(source, options), done)
        self.receptor_check_busy.emit(True)

    def prepare_receptor(self, options: ReceptorOptions) -> None:
        if self.source_pdb is None:
            self.error.emit("No structure", "Open a PDB/mmCIF structure first (File > Open Structure).")
            return
        project = self.ensure_project()
        out = project.prepared_dir("receptor")
        source = self.source_pdb
        self._submit(PrepareReceptorJob(source, out, options),
                     lambda result: self._receptor_done(result) if self.source_pdb == source else None)

    def _receptor_done(self, result: ReceptorResult) -> None:
        self._contact_queue.clear()
        self._cached_protein_rdkit = None
        self.receptor = result
        self.hetero = result.hetero or self.hetero
        self.report.emit("Receptor preparation", result.report.render())
        self.receptor_prepared.emit(result)
        self.hetero_found.emit(self.hetero)
        self.save_project()
        n_warn = len(result.report.warnings)
        self.message.emit(f"Receptor prepared ({result.n_atoms} atoms, {n_warn} warning(s))")

    # --------------------------------------------------------------- ligands
    def add_smiles(self, smiles: str, name: str) -> bool:
        entry = entry_from_smiles(smiles, name or f"ligand_{len(self.entries) + 1}")
        if entry.mol is None:
            self.error.emit("Invalid structure", entry.error)
            return False
        self._add_entries([entry])
        return True

    def import_ligand_file(self, path: Path) -> None:
        project = self.ensure_project()
        try:
            entries = read_ligands(Path(path))
        except PreparationError as exc:
            self.error.emit("Cannot read ligands", str(exc))
            return
        project.import_original(Path(path))
        bad = [e for e in entries if e.mol is None]
        if bad:
            self.report.emit("Ligand import", "\n".join(f"[ERROR  ] {e.name}: {e.error}" for e in bad))
        self._add_entries([e for e in entries if e.mol is not None])
        self.message.emit(f"Imported {len(entries) - len(bad)} ligand(s), {len(bad)} rejected")

    def _add_entries(self, entries: list[LigandEntry]) -> None:
        names = {e.name for e in self.entries}
        for e in entries:
            base, n = e.name, 2
            while e.name in names:
                e.name, n = f"{base}_{n}", n + 1
            names.add(e.name)
            self.entries.append(e)
        self.prepared = []
        self.ligands_changed.emit()

    def remove_ligand(self, name: str) -> None:
        self.entries = [e for e in self.entries if e.name != name]
        self.prepared = [p for p in self.prepared if p.name != name]
        self.ligands_changed.emit()

    def prepare_ligands(self, *, ionize: str, n_conformers: int, seed: int) -> None:
        if not self.entries:
            self.error.emit("No ligands", "Add a SMILES string or import a ligand file first.")
            return
        self.ligand_options = {"ionize": ionize, "n_conformers": n_conformers, "seed": seed}
        out = self.ensure_project().prepared_dir("ligands")
        self._submit(PrepareLigandsJob(list(self.entries), out, seed=seed, n_conformers=n_conformers,
                                       ionize=ionize), self._ligands_done)

    def _ligands_done(self, outcome: LigandPrepOutcome) -> None:
        self.prepared = outcome.prepared
        lines = []
        for lig in outcome.prepared:
            lines.append(f"=== {lig.name} ({lig.smiles}) ===\n{lig.report.render()}")
        for name, rep in outcome.failures:
            lines.append(f"=== {name} - FAILED ===\n{rep.render()}")
        self.report.emit("Ligand preparation", "\n\n".join(lines))
        self.ligands_prepared.emit(outcome)
        self.save_project()
        self.message.emit(f"{len(outcome.prepared)} ligand(s) prepared, {len(outcome.failures)} failed")

    # ---------------------------------------------------------------- box
    def set_box(self, box: DockingBox | None) -> None:
        self.box = box
        self.box_changed.emit(box)

    def box_from_hetero(self, hetero: HeteroResidue, padding: float) -> DockingBox:
        assert hetero.extracted_to is not None
        self._box_ligand_path = hetero.extracted_to
        box = boxlib.box_from_ligand_file(hetero.extracted_to, padding)
        self.set_box(box)
        return box

    def box_from_residues(self, residues: set[tuple[str, int]], padding: float) -> DockingBox:
        if self.receptor is None:
            raise ValueError("Prepare the receptor first.")
        box = boxlib.box_from_residues(self.receptor.prepared_pdb, residues, padding)
        self.set_box(box)
        return box

    def blind_box(self) -> DockingBox:
        if self.receptor is None:
            raise ValueError("Prepare the receptor first.")
        box = boxlib.blind_box(self.receptor.prepared_pdb)
        self.set_box(box)
        return box

    def set_box_reference(self, path: Path | None) -> None:
        self._box_ligand_path = path
        if self.box is not None:
            self.box_changed.emit(self.box)

    def validate_box(self) -> Report:
        if self.box is None:
            rep = Report()
            rep.error("No docking box defined.")
            return rep
        extent = None
        if self.prepared:
            import numpy as np
            extent = max(self._box_geometry.extent(
                np.asarray(p.mol.GetConformer().GetPositions())[
                    [a.GetAtomicNum() != 1 for a in p.mol.GetAtoms()]]) for p in self.prepared)
        pdbqt = self.receptor.pdbqt if self.receptor else None
        receptor_xyz = self._box_geometry.file_coords(pdbqt) if pdbqt else None
        coords, label = None, "Reference ligand"
        try:
            if self._box_ligand_path is not None:
                coords = self._box_geometry.file_coords(self._box_ligand_path, ligand=True)
                label = f"Co-crystal ligand {self._box_ligand_path.stem}"
            elif self.reference is not None and self.reference.GetNumConformers():
                coords = Chem.RemoveHs(self.reference).GetConformer().GetPositions()
                label = self.reference_name or label
        except (OSError, ValueError, RuntimeError) as exc:
            result = boxlib.validate_box(self.box, ligand_extent=extent, receptor_xyz=receptor_xyz)
            result.warn(f"Cannot check reference ligand fit: {exc}")
            return result
        return boxlib.validate_box(self.box, ligand_extent=extent, ligand_coords=coords,
                                   ligand_label=label, receptor_xyz=receptor_xyz)

    # ---------------------------------------------------------- reference
    def set_reference(self, hetero: HeteroResidue | None) -> None:
        if hetero is None or hetero.extracted_to is None:
            self.reference, self.reference_name = None, ""
            self.reference_changed.emit("")
            return
        cache = self.ensure_project().prepared_dir("cache")
        self._submit(_ReferenceJob(hetero, cache), lambda r, h=hetero: self._reference_done(h, r))

    def _reference_done(self, hetero: HeteroResidue, result: tuple[Chem.Mol, Chem.Mol]) -> None:
        template, ref = result
        self.reference, self.reference_name = ref, hetero.label
        self.reference_changed.emit(hetero.label)
        self.message.emit(f"Reference pose set from {hetero.label} (RMSD will be reported)")

    def redock_cocrystal(self, hetero: HeteroResidue) -> None:
        """Queue the co-crystal ligand (re-embedded from its chemistry) for a redocking test."""
        def done(result: tuple[Chem.Mol, Chem.Mol]) -> None:
            template, ref = result
            self.reference, self.reference_name = ref, hetero.label
            entry = LigandEntry(f"{hetero.name}_redock", Chem.Mol(template), "", "RCSB definition")
            self._add_entries([entry])
            self.box_from_hetero(hetero, 8.0)
            self.reference_changed.emit(hetero.label)
            self.message.emit(f"{hetero.name} queued for redocking; box set around the crystal ligand")
        cache = self.ensure_project().prepared_dir("cache")
        self._submit(_ReferenceJob(hetero, cache), done)

    # ------------------------------------------------------------- docking
    def run_docking(self, params: DockingParameters) -> None:
        if self.docking_running:
            self.message.emit("Docking is already running.")
            return
        problems = []
        if self.receptor is None:
            problems.append("Prepare the receptor first.")
        if not self.prepared:
            problems.append("Prepare at least one ligand first.")
        if self.box is None:
            problems.append("Define a docking box.")
        if problems:
            self.error.emit("Cannot start docking", "\n".join(problems))
            return
        box_report = self.validate_box()
        if box_report.errors:
            self.error.emit("Invalid docking box", "\n".join(box_report.errors))
            return
        self.params = params
        project = self.ensure_project()
        self.run_dir = project.new_run_dir()
        job = DockingJob(self.receptor.pdbqt, list(self.prepared), self.box, params, self.run_dir,  # type: ignore[union-attr]
                         self.reference, self.source_pdb.name if self.source_pdb else "")
        self._docking_job = self._submit(job, self._docking_done)
        self.docking_started.emit()
        self.message.emit("Docking started")

    def cancel_docking(self) -> None:
        if self._docking_job is not None:
            self.runner.cancel(self._docking_job)

    def _docking_done(self, outcomes: list[LigandOutcome]) -> None:
        self._docking_job = None
        self.outcomes = outcomes
        project = self.ensure_project()
        if self.run_dir is not None:
            runs = project.state.setdefault("runs", [])
            rel = project.rel(self.run_dir)
            if rel not in runs:
                runs.append(rel)
        self.docking_finished.emit(outcomes)
        self.save_project()
        ok = sum(1 for o in outcomes if o.status == "ok")
        self.message.emit(f"Docking finished: {ok}/{len(outcomes)} ligand(s) succeeded")

    @property
    def docking_running(self) -> bool:
        return self._docking_job is not None

    # -------------------------------------------------------------- analysis
    def contacts_for(self, outcome: LigandOutcome, pose_index: int) -> list[Contact]:
        if self.receptor is None:
            return []
        atoms = cached_protein_atoms(self.receptor.prepared_pdb)
        return pose_contacts(outcome.mols[pose_index], atoms)

    def request_contacts(self, outcome: LigandOutcome, pose_index: int) -> None:
        if self.receptor is None:
            self._contact_queue.clear()
            self.contacts_ready.emit(outcome, pose_index, [])
            return
        receptor, project = self.receptor, self.project
        source = receptor.prepared_pdb

        def current():
            return self.project is project and self.receptor is receptor

        def done(contacts):
            if current():
                self.contacts_ready.emit(outcome, pose_index, contacts)

        def failed(message):
            if current():
                lines = message.strip().splitlines()
                self.contacts_failed.emit(outcome, pose_index, lines[-1] if lines else "Unknown error")

        self._contact_queue.submit(ContactsJob(source, outcome.mols[pose_index]), done, failed)

    def clear_contact_requests(self) -> None:
        self._contact_queue.clear()

    def analyze_interactions(self, outcome: LigandOutcome, pose_index: int) -> None:
        """ProLIF interaction fingerprint + 2D diagram for one pose, run in the background.

        The receptor's one-time RDKit rebuild (Open Babel bond perception) is cached across
        calls for the same prepared receptor, so later poses analyse in well under a second.
        """
        if self.receptor is None:
            self.error.emit("No receptor", "Prepare the receptor first.")
            return
        if self._interaction_job is not None:
            return  # a previous request is still running; its result will arrive shortly
        cached = None
        if self._cached_protein_rdkit is not None and self._cached_protein_rdkit[0] == self.receptor.prepared_pdb:
            cached = self._cached_protein_rdkit[1]
        self.interactions_started.emit()
        job = InteractionJob(self.receptor.prepared_pdb, outcome.mols[pose_index], cached_protein=cached)
        self._interaction_job = self.runner.submit(job)
        self._callbacks[self._interaction_job] = lambda r, o=outcome, i=pose_index: self._interactions_done(o, i, r)

    def _interactions_done(self, outcome: LigandOutcome, pose_index: int, result: tuple) -> None:
        self._interaction_job = None
        protein, interaction_result = result
        self._cached_protein_rdkit = (self.receptor.prepared_pdb, protein)
        self.interactions_ready.emit(outcome, pose_index, interaction_result)

    def compute_similarity(self, cutoff: float, reference_smiles: str = "") -> None:
        """Tanimoto similarity + clustering of the prepared ligands, run in the background."""
        if len(self.prepared) < 2:
            self.error.emit("Not enough ligands", "Prepare at least two ligands first.")
            return
        reference = None
        reference_smiles = reference_smiles.strip()
        if reference_smiles:
            reference = Chem.MolFromSmiles(reference_smiles)
            if reference is None:
                self.error.emit("Invalid reference", f"Could not parse SMILES: {reference_smiles}")
                return
        mols = [p.mol for p in self.prepared]
        names = [p.name for p in self.prepared]
        job = SimilarityJob(mols, names, cutoff, reference, reference_smiles)
        self._submit(job, self.similarity_ready.emit)

    # -------------------------------------------------------------- selection
    def set_selection(self, detail: dict, summary: str) -> None:
        """Called when the viewer emits a selection event."""
        self.selection = detail or {}
        self.selection_changed.emit(self.selection, summary)

    def clear_selection(self) -> None:
        self.selection = {}
        self.selection_changed.emit({}, "")
        self.selection_cleared.emit()

    def expand_selection_to_chain(self) -> None:
        """Hand control to the viewer: it knows which chains each selected residue is in."""
        self.selection_expand_requested.emit()

    def invert_selection(self) -> None:
        self.selection_invert_requested.emit()

    def _edit_source_pdb(self, residues: list[dict], *, keep: bool) -> None:
        """Write an immutable structural revision and activate it in the viewer."""
        from uuid import uuid4
        from mds.chem.edit import edit_atoms
        if self.project is None or self.source_pdb is None or not residues:
            return
        src = self.source_pdb
        edits_dir = self.project.prepared_dir("edited_structure")
        suffix = ".cif" if src.suffix.lower() in (".cif", ".mmcif") else ".pdb"
        dst = edits_dir / f"structure_{uuid4().hex[:12]}_edited{suffix}"
        try:
            n = edit_atoms(src, dst, residues, keep=keep)
        except (ValueError, RuntimeError, OSError) as exc:
            self.error.emit("Cannot edit structure", str(exc))
            return
        if not n:
            self.message.emit("No atoms changed.")
            return
        self._undo_structures.append(src)
        self._redo_structures.clear()
        self._activate_structure(dst)
        self.message.emit(f"Removed {n} atom(s). Undo is available; prepare the receptor again before docking.")

    def _activate_structure(self, path: Path) -> None:
        import gemmi
        previous = self.source_pdb
        st = gemmi.read_structure(str(path))
        self.source_pdb = path
        self.receptor = None
        self._cached_protein_rdkit = None
        self.reference, self.reference_name = None, ""
        self.reference_changed.emit("")
        self.set_box(None)
        self.hetero = classify_hetero(st) if len(st) else []
        self.clear_selection()
        self.hetero_found.emit(self.hetero)
        self.structure_edited.emit(previous.stem if previous else "", str(path))
        self.edit_history_changed.emit(bool(self._undo_structures), bool(self._redo_structures))
        self.save_project()

    def undo_structure_edit(self) -> None:
        if self._undo_structures and self.source_pdb:
            self._redo_structures.append(self.source_pdb)
            self._activate_structure(self._undo_structures.pop())

    def redo_structure_edit(self) -> None:
        if self._redo_structures and self.source_pdb:
            self._undo_structures.append(self.source_pdb)
            self._activate_structure(self._redo_structures.pop())

    def cleanup_structure(self, *, waters=False, hetero=False, ligands=False,
                          remove_chains=(), retain_selected=True) -> None:
        from mds.chem.edit import cleanup_selection
        if self.source_pdb is None:
            return
        retained = self.selection.get(self.source_pdb.stem, []) if retain_selected else []
        try:
            victims = cleanup_selection(self.source_pdb, waters=waters, hetero=hetero,
                                        ligands=ligands, remove_chains=remove_chains, retain=retained)
        except (ValueError, RuntimeError, OSError) as exc:
            self.error.emit("Cannot clean up structure", str(exc))
            return
        if victims:
            self._edit_source_pdb(victims, keep=False)
        else:
            self.message.emit("No components match the cleanup options.")

    def delete_selected_residues(self) -> None:
        source_name = self.source_pdb.stem if self.source_pdb else ""
        residues: list[dict] = []
        for name, res_list in self.selection.items():
            if name == source_name:
                residues.extend(res_list)
        if not residues:
            self.error.emit("Delete selection", "Select residues on the loaded protein first.")
            return
        self._edit_source_pdb(residues, keep=False)

    def keep_only_selected_residues(self) -> None:
        source_name = self.source_pdb.stem if self.source_pdb else ""
        residues: list[dict] = []
        for name, res_list in self.selection.items():
            if name == source_name:
                residues.extend(res_list)
        if not residues:
            self.error.emit("Keep only selection", "Select residues on the loaded protein first.")
            return
        self._edit_source_pdb(residues, keep=True)

    # --------------------------------------------------------------- export
    def export_csv(self, path: Path) -> None:
        write_csv(path, results_rows(self.outcomes))

    def export_excel(self, path: Path) -> None:
        import pandas as pd
        rows = results_rows(self.outcomes)
        note = pd.DataFrame({"Note": [
            "Scores are predicted by AutoDock Vina's empirical scoring function (kcal/mol).",
            "They are not measured binding affinities and are not evidence of biological activity."]})
        with pd.ExcelWriter(path) as xl:
            pd.DataFrame(rows).to_excel(xl, sheet_name="Docking results", index=False)
            note.to_excel(xl, sheet_name="Read me", index=False)

    def export_sdf(self, path: Path) -> None:
        from mds.chem.poses import write_poses_sdf
        mols, props = [], []
        for o in self.outcomes:
            for i, (p, m) in enumerate(zip(o.poses, o.mols)):
                r = o.reference_rmsd[i] if i < len(o.reference_rmsd) else None
                mols.append(m)
                props.append({"Name": o.name, "Rank": p.rank, "Score_kcal_per_mol": p.score,
                              "RMSD_to_reference": r if r is not None else ""})
        write_poses_sdf(Path(path), mols, props)

    # ----------------------------------------------------- project persistence
    def save_project(self) -> None:
        if self.project is None:
            return
        pr, st = self.project, self.project.state
        st["source_structure"] = pr.rel(self.source_pdb) if self.source_pdb else None
        st["receptor"] = ({"pdbqt": pr.rel(self.receptor.pdbqt),
                           "prepared_pdb": pr.rel(self.receptor.prepared_pdb),
                           "n_atoms": self.receptor.n_atoms} if self.receptor else None)
        st["ligands"] = [{"name": e.name, "smiles": Chem.MolToSmiles(Chem.RemoveHs(e.mol))}
                         for e in self.entries if e.mol is not None]
        st["ligand_options"] = self.ligand_options
        st["box"] = self.box.as_dict() if self.box else None
        st["parameters"] = vars(self.params)
        pr.save()

    def load_project_state(self) -> None:
        """Rebuild the workflow from ``project.json`` (called after opening a project)."""
        if self.project is None:
            return
        pr, st = self.project, self.project.state
        self.reset_state()
        if st.get("source_structure"):
            src = pr.abs(st["source_structure"])
            if src.exists():
                import gemmi
                self.source_pdb = src
                structure = gemmi.read_structure(str(src))
                self.hetero = classify_hetero(structure) if len(structure) else []
                self.structure_loaded.emit(src.stem, str(src))
        rec = st.get("receptor")
        if rec and pr.abs(rec["pdbqt"]).exists():
            lig_dir = pr.abs("prepared/receptor/extracted_ligands")
            hetero = list(self.hetero)
            for h in hetero:
                f = lig_dir / f"{h.name}_{h.chain}{h.number}.pdb"
                if h.kind == "ligand" and f.exists():
                    h.extracted_to = f
            self.receptor = ReceptorResult(pr.abs(rec["pdbqt"]), pr.abs(rec["prepared_pdb"]), Report(),
                                           hetero, rec.get("n_atoms", 0))
            self.receptor_prepared.emit(self.receptor)
        self.hetero_found.emit(self.hetero)
        for item in st.get("ligands", []):
            e = entry_from_smiles(item["smiles"], item["name"])
            if e.mol is not None:
                self.entries.append(e)
        self.ligands_changed.emit()
        if st.get("box"):
            b = st["box"]
            self.set_box(DockingBox((b["center_x"], b["center_y"], b["center_z"]),
                                    (b["size_x"], b["size_y"], b["size_z"])))
        if st.get("parameters"):
            self.params = DockingParameters(**st["parameters"])
        self.ligand_options = st.get("ligand_options", self.ligand_options)
        runs = st.get("runs", [])
        if runs:
            self._load_run(pr.abs(runs[-1]))
        self.project_changed.emit()

    def _load_run(self, run_dir: Path) -> None:
        import csv
        csv_path = run_dir / "results.csv"
        if not csv_path.exists():
            return
        self.run_dir = run_dir
        outcomes: dict[str, LigandOutcome] = {}
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        from mds.chem.workflow import _safe
        from mds.engines.base import Pose
        for row in rows:
            o = outcomes.setdefault(row["Ligand"], LigandOutcome(row["Ligand"], row.get("Status", "ok"),
                                                                 row.get("Message", ""),
                                                                 row.get("SMILES", ""),
                                                                 folder=run_dir / _safe(row["Ligand"])))
            if row.get("Status") == "ok" and row.get("Pose"):
                o.poses.append(Pose(int(row["Pose"]), float(row["Score (kcal/mol)"]),
                                    float(row.get("RMSD lb") or 0), float(row.get("RMSD ub") or 0), "",
                                    int(row.get("Conformer", 1)) - 1))
                r = row.get("RMSD to reference (A)")
                o.reference_rmsd.append(float(r) if r else None)
        for o in outcomes.values():
            sdf = (o.folder or run_dir) / "poses.sdf"
            if sdf.exists():
                o.mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]
        self.outcomes = [o for o in outcomes.values() if o.status != "ok" or len(o.mols) == len(o.poses)]
        self.docking_finished.emit(self.outcomes)
