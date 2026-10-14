# WaweDock — Architecture & Roadmap

## 1. Guiding principles

1. **Science before decoration.** No feature is presented as operational until it performs a real, tested computation. Unbuilt features are visibly labelled "Planned — Phase N" and disabled.
2. **Layered architecture.** Scientific code never imports Qt widgets; GUI code never contains scientific algorithms.
3. **Adapters for anything replaceable** — docking engines, renderers, interaction analysers, structure-repair tools.
4. **Reproducibility by default.** Every computation records inputs, parameters, software versions and random seeds.
5. **The GUI thread never computes.** Anything longer than ~100 ms runs as a background job.

## 2. Layers

```
┌──────────────────────────────────────────────────────────────┐
│ GUI (mds.gui)            PySide6 widgets, docks, dialogs     │
│   └─ talks to services only through signals / plain calls    │
├──────────────────────────────────────────────────────────────┤
│ Application services     workflows: "prepare receptor",      │
│ (mds.services)           "run docking", "screen library"     │
├──────────────┬───────────────┬───────────────┬───────────────┤
│ Scientific   │ Engine        │ Rendering     │ Persistence   │
│ core         │ adapters      │ adapters      │ (mds.project) │
│ (mds.chem)   │ (mds.engines) │ (mds.render)  │ SQLite + files│
├──────────────┴───────────────┴───────────────┴───────────────┤
│ Infrastructure: jobs (mds.jobs), logging (mds.core), I/O     │
└──────────────────────────────────────────────────────────────┘
```

### Planned package map (packages are created only when their phase starts)

| Package            | Responsibility                                             | Phase |
|--------------------|------------------------------------------------------------|-------|
| `mds.core`         | logging, configuration, shared types, errors               | 1     |
| `mds.jobs`         | background job runner, cancellation, progress              | 1     |
| `mds.gui`          | main window, panels, dialogs, themes, layouts              | 1+    |
| `mds.render`       | `ViewerBackend` interface + concrete backend               | 2     |
| `mds.model`        | molecular object model (structures, selections, boxes)     | 3     |
| `mds.io`           | readers/writers: PDB, mmCIF, SDF, MOL2, SMILES, PDBQT      | 3     |
| `mds.project`      | project folder, SQLite metadata, provenance records        | 3     |
| `mds.prep.protein` | receptor preparation pipeline + report                     | 4     |
| `mds.prep.ligand`  | ligand preparation pipeline (RDKit/Meeko)                  | 5     |
| `mds.engines`      | `DockingEngine` interface; `VinaEngine` (later Smina/Gnina)| 6–7   |
| `mds.analysis`     | interactions, fingerprints, RMSD, descriptors, clustering  | 8, 10 |
| `mds.screening`    | queues, campaigns, resume logic                            | 9     |
| `mds.validation`   | redocking benchmarks, reproducibility tests                | 11    |

## 3. Technology compatibility strategy (verified on this machine, Python 3.12, Windows 11)

| Component          | Source                         | Status / decision                                                    | Licence      |
|--------------------|--------------------------------|----------------------------------------------------------------------|--------------|
| PySide6            | pip                            | Installed (6.12). GUI framework.                                     | LGPL-3       |
| RDKit              | pip wheel                      | Available. Ligand chemistry, descriptors, fingerprints.              | BSD-3        |
| Meeko              | pip                            | Available (0.8). Ligand + receptor PDBQT preparation.                | LGPL-2.1     |
| Gemmi              | pip wheel                      | Available. Fast PDB/mmCIF parsing, structure inspection.             | MPL-2.0      |
| BioPython          | pip wheel                      | Available. Secondary structure helpers, sequence handling.           | BSD-like     |
| **AutoDock Vina**  | **official `vina.exe` release**| **No Windows wheel for cp312 on PyPI.** Run the executable via a subprocess job → also gives clean cancellation (kill the process). | Apache-2.0 |
| ProLIF             | pip                            | In core deps. Interaction fingerprints + 2D diagram (mds/chem/interactions.py). | Apache-2.0 |
| Open Babel         | pip (`openbabel-wheel`)        | **Optional extra** (`wavedock[interactions]`). Bond perception for ProLIF's protein model; excluded from the PyInstaller build. | GPL-2.0 ⚠ |
| PLIP               | —                               | Not used (ProLIF covers interaction fingerprinting; avoids a second GPL dependency). | GPL-2.0 ⚠ |
| PyMOL open-source  | pip wheel                      | Available. Candidate viewer backend (evaluated in Phase 2).          | BSD-like     |
| VTK                | pip wheel                      | Available. Alternative viewer backend.                               | BSD-3        |
| Pandas / NumPy     | pip                            | Available.                                                           | BSD-3        |
| PyQtGraph          | pip                            | Available. Fast in-GUI plots.                                        | MIT          |
| PyInstaller        | pip                            | Phase 12.                                                            | GPL + exception |

⚠ **Licence note:** PLIP and Open Babel are GPL. Bundling them into a distributed executable makes the distribution GPL-bound. ProLIF (Apache-2.0) is the interaction analyser; the one GPL dependency it needs from this app — Open Babel, for building a chemically-correct protein model — is an optional extra (`pip install wavedock[interactions]`) that is explicitly excluded from the PyInstaller build (`wavedock.spec`). Without it, interaction analysis reports a clear message instead of silently failing or faking a result.

Optional tools are detected at runtime (`shutil.which`, import probes) and the corresponding features are disabled with an explanatory tooltip when absent — never silently faked.

## 4. Key design decisions

### 4.1 Background jobs (`mds.jobs`)
* `Job` (abstract) — implement `run(ctx)`; report progress with `ctx.progress()`; call `ctx.check_cancelled()` at safe points.
* `JobRunner` — owns a `QThreadPool`; re-emits job events as Qt signals on the GUI thread.
* External programs (Vina, structure repair tools) run as subprocesses inside a job; cancellation terminates the process.
* CPU-heavy pure-Python work will move to process pools when needed (GIL); RDKit/Gemmi release the GIL for most heavy operations.

### 4.2 Docking engine adapters (Phase 6–7)
```python
class DockingEngine(Protocol):
    name: str
    def version(self) -> str: ...
    def validate(self, request: DockingRequest) -> list[ValidationIssue]: ...
    def run(self, request: DockingRequest, ctx: JobContext) -> DockingResult: ...
```
Results always carry `engine`, `engine_version`, `seed`, `score_units="kcal/mol"` and `score_kind="Vina scoring function"` — scores are never labelled as experimental affinities or as other engines' scores (e.g. CDOCKER).

### 4.3 Rendering (Phase 2)
GUI code talks to a `ViewerBackend` interface (load object, set representation, colour, select, measure, draw box, screenshot). The concrete backend is chosen after a Phase-2 evaluation of embedded PyMOL vs VTK.

### 4.4 Projects (Phase 3)
```
MyProject.mdsproj/
  project.sqlite        metadata, provenance, parameters, annotations
  originals/            imported files — read-only, never overwritten
  prepared/             receptors / ligands produced by preparation
  runs/<run-id>/        engine inputs, logs, poses
  reports/              exports
```

## 5. Development phases

| # | Phase                                                   | Exit criterion                                            |
|---|---------------------------------------------------------|-----------------------------------------------------------|
| 1 | Desktop foundation                                      | Window, docks, themes, layouts, logging, jobs — tested    |
| 2 | Embedded 3D viewer                                      | Load PDB, rotate/zoom, representations, colouring          |
| 3 | File import, object model, project saving               | Round-trip a project with proteins + ligands               |
| 4 | Protein preparation                                     | Receptor PDBQT + preparation report                       |
| 5 | Ligand preparation                                      | Batch ligand PDBQT with descriptors & warnings            |
| 6 | Docking box                                             | Box from ligand/residues, interactive editing, validation |
| 7 | AutoDock Vina integration                               | Redock a known complex, RMSD < 2 Å                        |
| 8 | Results & interaction analysis                          | Pose table, interactions, 2D diagram, exports             |
| 9 | Batch docking & virtual screening                       | Resumable multi-ligand campaign                           |
| 10| Descriptors, drug-likeness, analytics                   | Ro5/Veber, similarity, clustering                          |
| 11| Validation, testing, performance                        | Benchmark suite passes                                    |
| 12| Windows packaging                                       | PyInstaller build runs on a clean machine                 |

## 6. Phase 1 roadmap

| Step | Deliverable                                    | Files                                   |
|------|------------------------------------------------|-----------------------------------------|
| 1.1  | Environment, package skeleton, entry point     | `pyproject.toml`, `mds/__init__.py`, `mds/__main__.py`, `mds/app.py` |
| 1.2  | Structured logging + crash handler             | `mds/core/logging_setup.py`             |
| 1.3  | Background job framework                       | `mds/jobs/runner.py`, `mds/jobs/diagnostics.py` |
| 1.4  | Main window, menus, toolbars, shared actions   | `mds/gui/main_window.py`, `mds/gui/actions.py` |
| 1.5  | Dock panels                                    | `mds/gui/panels/*.py`                   |
| 1.6  | Central viewport host (filled in Phase 2)      | `mds/gui/viewport.py`                   |
| 1.7  | Status bar with job progress                   | `mds/gui/status_bar.py`                 |
| 1.8  | Light / dark themes                            | `mds/gui/theme.py`                      |
| 1.9  | Layout save / restore / named layouts          | `mds/gui/layouts.py`                    |
| 1.10 | Drag-and-drop file reception                   | `main_window.py`                        |
| 1.11 | Automated tests                                | `tests/`                                |

## 7. Status (2026-10-14)

Delivered ahead of the original phase order at the user's request:

| Area | Implementation | Module |
|---|---|---|
| 3D viewer | 3Dmol.js in QtWebEngine behind `ViewerBackend` | `mds/render/` |
| Receptor prep | Gemmi → PDB2PQR/PROPKA → Meeko, full report | `mds/chem/receptor.py` |
| Ligand prep | RDKit (ETKDGv3, MMFF94) → Meeko | `mds/chem/ligand.py` |
| Docking box | ligand / residues / blind / manual + validation | `mds/chem/box.py` |
| Vina | `vina.exe` subprocess adapter, cancellable | `mds/engines/vina.py` |
| Poses & RMSD | Meeko → RDKit, symmetry-aware RMSD | `mds/chem/poses.py` |
| Contacts | distance-based candidates (conservative) | `mds/chem/contacts.py` |
| Interactions (Phase 8) | ProLIF fingerprint (H-bond, hydrophobic, pi-stacking, salt bridge, ...) + 2D network diagram, fully offline | `mds/chem/interactions.py` |
| Projects | folder + `project.json`, run folders with `run.json` | `mds/project.py` |
| Packaging | PyInstaller one-folder build | `wavedock.spec` |

Environment findings: OpenMM/PDBFixer and RDKit's `rdMolStandardize` DLLs are blocked by Smart
App Control on the development PC, so they are not used. Open Babel (GPL) is needed to build a
chemically-correct protein model for ProLIF (Gemmi's bundled monomer templates don't cover
PDB2PQR's protonation-variant residue names like HID/HIE/ASH/CYX) — kept as an optional,
separately-installed extra, excluded from the PyInstaller build; see §3's licence note. Still
open: screening queues/resume (Phase 9), similarity/clustering (Phase 10), flexible side chains,
code signing (Phase 12).
