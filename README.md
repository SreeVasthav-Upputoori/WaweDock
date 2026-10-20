# WaweDock

Open-source molecular docking workbench for Windows: protein preparation, ligand preparation,
AutoDock Vina docking, interactive 3D visualisation and result analysis.
Architecture and roadmap: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Start

* **Desktop shortcut** `WaweDock` (also `WaweDock.lnk` in this folder), or
* from PowerShell:

```powershell
cd "D:\Tools\Docking tool"
.\.venv\Scripts\Activate.ps1
python -m mds
```

## Typical workflow

1. **File ▸ Fetch from PDB** (e.g. `1HSG`) or **Open Structure** (PDB/mmCIF, or drag-and-drop).
2. **Docking Setup ▸ 1 Receptor ▸ Prepare Receptor** — waters removed, metals kept, missing atoms
   rebuilt, protonation at pH 7.4 (PDB2PQR/PROPKA), PDBQT built (Meeko). Every decision is listed in
   the *Preparation Report* tab.
3. **2 Ligands** — add SMILES or import SDF/MOL/MOL2/SMI/CSV → **Prepare Ligands** (RDKit 3D + MMFF,
   Meeko PDBQT, descriptors, stereo/ionisation warnings). *Redock co-crystal ligand* sets up a validation run.
4. **3 Binding site** — box from co-crystal ligand, residues, whole protein, or manual; drawn in 3D and validated.
5. **4 AutoDock Vina ▸ Run Docking** (F5) — runs in the background, cancellable.
6. **Results** — ranked poses (kcal/mol), RMSD to crystal pose when a reference exists, distance-based
   contact candidates, 3D pose view; export CSV / Excel / SDF. Every run folder contains `run.json`
   (versions, seeds, parameters, commands) for reproducibility.
7. **Results ▸ Interactions (ProLIF)** (or **Analyse ▸ 2D Interaction Diagram**) — classifies contacts
   for the selected pose (H-bond, hydrophobic, pi-stacking, salt bridge, metal coordination, ...) and
   draws a 2D interaction diagram, entirely offline. **Needs Open Babel** (see below); unavailable
   otherwise, with a clear message rather than a silent failure or a faked result.
8. **Analyse ▸ Similarity & Clustering** — Tanimoto similarity on Morgan fingerprints and Butina
   clustering across the prepared ligands, with an optional reference compound to rank against.

Docking scores are predictions from Vina's empirical scoring function — not measured affinities.
Interaction types are geometric/chemical classifications of a predicted pose (ProLIF), not
experimentally verified interactions.

## Interactive structure editing

Opening a PDB or mmCIF displays all chains, waters, ligands, ions, hydrogens and alternate
locations. Cleanup happens only when requested; receptor preparation remains a separate step.

Use the **Selection & Cleanup** panel on the left:

- Choose atom, residue, chain, or component selection. Components select a hetero residue
  (such as a ligand or ion), or a protein chain.
- Click to select; Shift/Ctrl-click toggles additional components. **Alt-drag** draws a freehand
  lasso. **Shift+Alt-drag** adds to the selection. The lasso selects through the structure's depth;
  rotate the view or use a 3D radius to isolate a spatial region.
- Press **Delete**, use the viewport's right-click menu, or click **Delete selected atoms**.
  **Keep only selected** removes everything else. These operations change the structure files
  used by subsequent receptor preparation, with each revision stored in the project.
- Set a radius in Å and click **Select within radius of selection** to select around the chosen
  atoms in 3D. The selection mode determines whether hits expand to whole residues or chains.
- Cleanup can remove waters, HETATM components, non-protein ligands, or named chains. Enable
  **Retain selected atoms / components** to protect selected items during cleanup.
- **Undo/Redo** buttons and **Ctrl+Z / Ctrl+Y** reverse structural edits during the current session.
  Original imported files remain unchanged. After editing, prepare the receptor again before docking.

Prepared receptors and docking poses are analysis overlays. Click **Show imported structure for
editing** to return to the editable source after preparation.

## Validation (redocking)

```powershell
python -m mds --validate tests\data\1HSG.pdb MK1 "<SMILES>" out_dir
```

| Complex | Ligand | Top-pose RMSD to crystal |
|---|---|---|
| 3PTB trypsin | benzamidine (amidinium) | 0.37 Å |
| 1HSG HIV protease | indinavir, neutral input | 1.16 Å |
| 1HSG HIV protease | indinavir, simple pH 7.4 rules | 2.04 Å (borderline) |

## Optional: interaction analysis (Open Babel)

The "Interactions (ProLIF)" feature needs Open Babel to build a chemically-correct protein model.
It is GPL-2.0 licensed, so it is an optional extra, not a core dependency, and is **excluded from
the packaged `.exe`** (see Licensing below):

```powershell
pip install wavedock[interactions]      # or: pip install openbabel-wheel
```

## Tests

```powershell
pytest          # 32 tests: preparation, Vina adapter, redocking, reproducibility, interactions, GUI workflow
```

## Windows executable

```powershell
pip install pyinstaller
pyinstaller --noconfirm wavedock.spec     # -> dist\WaweDock\WaweDock.exe (+ _internal\)
```

The build bundles AutoDock Vina 1.2.7, the monomer library and the 3D viewer, but deliberately
**not** Open Babel (see Licensing), so "Interactions (ProLIF)" is unavailable in this build.
**Note:** Windows 11 *Smart App Control* blocks new unsigned executables (this machine blocks
`WaweDock.exe` for that reason). Sign the executable with a code-signing certificate before
distributing it, or run from source as above.

Logs: `%LOCALAPPDATA%\WaweDock\WaweDock\logs\mds.log` (*Tools ▸ Open Log Folder*).
Logo source: `tools/make_logo.py`.

## Third-party components

AutoDock Vina (Apache-2.0), RDKit (BSD-3), Meeko (LGPL-2.1), Gemmi (MPL-2.0), PDB2PQR (BSD-3),
PROPKA (LGPL-2.1), ProLIF (Apache-2.0), 3Dmol.js (BSD-3), vis-network (Apache-2.0/MIT),
PySide6/Qt (LGPL-3), CCP4 monomer library (LGPL-3).

**Licensing note:** Open Babel (GPL-2.0) is the one optional dependency that is GPL — it powers
protein bond-perception for interaction analysis only. Bundling GPL code into a distributed
executable would make that distribution GPL-bound, so it is never included in `wavedock.spec` and
must be installed separately for that feature to work from source.
