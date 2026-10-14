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

Docking scores are predictions from Vina's empirical scoring function — not measured affinities.

## Validation (redocking)

```powershell
python -m mds --validate tests\data\1HSG.pdb MK1 "<SMILES>" out_dir
```

| Complex | Ligand | Top-pose RMSD to crystal |
|---|---|---|
| 3PTB trypsin | benzamidine (amidinium) | 0.37 Å |
| 1HSG HIV protease | indinavir, neutral input | 1.16 Å |
| 1HSG HIV protease | indinavir, simple pH 7.4 rules | 2.04 Å (borderline) |

## Tests

```powershell
pytest          # 30 tests: preparation, Vina adapter, redocking, reproducibility, GUI workflow
```

## Windows executable

```powershell
pip install pyinstaller
pyinstaller --noconfirm wavedock.spec     # -> dist\WaweDock\WaweDock.exe (+ _internal\)
```

The build bundles AutoDock Vina 1.2.7, the monomer library and the 3D viewer.
**Note:** Windows 11 *Smart App Control* blocks new unsigned executables (this machine blocks
`WaweDock.exe` for that reason). Sign the executable with a code-signing certificate before
distributing it, or run from source as above.

Logs: `%LOCALAPPDATA%\WaweDock\WaweDock\logs\mds.log` (*Tools ▸ Open Log Folder*).
Logo source: `tools/make_logo.py`.

## Third-party components

AutoDock Vina (Apache-2.0), RDKit (BSD-3), Meeko (LGPL-2.1), Gemmi (MPL-2.0), PDB2PQR (BSD-3),
PROPKA (LGPL-2.1), 3Dmol.js (BSD-3), PySide6/Qt (LGPL-3), CCP4 monomer library (LGPL-3).
