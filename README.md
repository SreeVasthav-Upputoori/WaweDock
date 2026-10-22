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

On Windows, WaweDock automatically adds itself to **Open with** on first launch,
for PDB/mmCIF and supported molecular ligand files. This works in the packaged
app and source/pip installations; it uses the installation's own executable or
Python environment, with no developer-machine paths. Existing default apps stay
unchanged. Reopening WaweDock repairs registration after moving its folder.

Right-click a supported file and choose **Open with > WaweDock** (it may first
appear under **Choose another app**). You can also launch directly with
`python -m mds "C:\path\protein.pdb"`. To repair registration manually, run
`python -m mds --register-file-associations`.

## Receptor checks

In **Prepare > Receptor Preparation**, use **Check Receptor** before preparing.
The read-only background check reviews the selected chains for incomplete backbone
atoms, chain breaks and alternate conformations, and describes water removal,
metal retention, excluded hetero groups/cofactors and the selected protonation pH.
Open **Window > Panels > Preparation Report** for full details. Changing preparation
options clears the previous check summary. The check does not rebuild atoms or
predict pKa; review the prepared receptor as well. Explicit-water retention is not
yet supported and is now rejected rather than silently ignored.

## Docking-box preview

The box preview follows manual centre/size edits and the selected co-crystal ligand.
The sidebar reports the number of heavy atoms outside the box, maximum overrun,
and minimum face clearance for an enclosed reference pose. Amber indicates a
warning and red a validation error; green means no reported issues. A reference
with less than 2 A face clearance gets an advisory, not an automatic box change.
Prepared ligands are checked by size, without treating their arbitrary starting
coordinates as binding-site poses. Without a positioned reference, the report
explicitly says containment was not checked.

## Workspace controls

The molecular viewer occupies the main workspace, with a compact Project Explorer
on the left and Docking Setup on the right. Panels can be resized, closed, or
detached by dragging their title bars. Reopen panels through **Window > Panels**;
Preparation Report, Properties, Selection & Cleanup, Log and Jobs are
available on demand. **Window > Reset Layout** restores this arrangement.

Right-click a molecule directly in the viewer for representation, colour, labels,
hydrogen/water/ligand visibility and molecular operations. These controls apply to
the clicked structure immediately. Custom molecular and background colours use a
colour picker. Right-click empty space for background settings, camera reset and
**Show / Hide Structures**, which can restore a hidden molecule. Project Explorer
also retains its representation and colour menus. There is no Appearance sidebar.
Preparation, measurements and image export remain in the top-level menus.
Ligands and docking poses have **Color Options > Carbon Color** (preserves other
element colours) and **Solid Ligand Color**, with preset and custom colours plus
reset. For ligands embedded in a protein, use **Ligand Color Options** to colour
the bound ligands independently of the protein and waters. Selection remains
yellow until cleared, then the chosen ligand colours return.

Changing the viewport background (black, white, dark grey or a custom colour)
does not change the application's interface colours.

Single left-click selects one amino-acid residue or one water molecule;
double-click selects the whole protein (including modified amino acids) or all
waters in the clicked structure. Selection recolours the existing representation
opaque yellow without adding spheres or changing its size, including on prepared
read-only overlays. Deselecting restores the original colours. Shift/Ctrl-click
toggles additional residues. The Selection & Cleanup mode controls lasso and
expansion tools. Turning off
**Show Labels** removes molecular labels, including when nearby-residue display
is active; measurement labels remain until measurements are cleared.

Right-clicking preserves the current selection. **Select Clicked…** selects the
clicked atom, residue, chain or component; Shift/Ctrl-click and Alt-drag lasso still
work as before. **Delete Selected Atoms / Residues** edits the imported structure
and supports Undo/Redo. Prepared receptors and docking poses remain read-only.
Hydrogen, water and ligand visibility switches only change the view; they do not
remove atoms from the project or change docking inputs.

Docking Setup is a five-step workflow: **Receptor Preparation**, **Ligand
Preparation**, **Docking Box Configuration**, **Docking Execution**, and **Results
Review**. Click a step header to expand it; click it again to collapse it. Only one
step shows its detailed settings at a time, and changing steps preserves inputs.
Headers report completion, partial ligand preparation, invalid boxes, running
jobs, cancellation and failure from the current project state. The execution step
lists missing prerequisites; Review Results opens Results & Analysis.
Preparation and docking continue to run in background jobs.

Use the top-level **Prepare | Dock | Results & Analysis** switcher to change
workspaces. Results & Analysis places a searchable, numerically sortable pose
table beside the same molecular viewer, with a resizable divider. Selecting a row
shows its actual docking pose. Ctrl/Shift-select rows and open **Compare Poses**
to compare their scores and RMSD metrics; **Score Comparison** charts the best
score per ligand, and clicking a bar selects that pose. Contacts and ProLIF
interaction analysis remain available in the analysis tabs. The **Export** menu
provides CSV, Excel, SDF and a standalone HTML metrics report; exports include all
results regardless of the search filter. Workspace switching preserves the
viewer, camera and molecular selection, and restores the setup panel arrangement
when returning to Prepare or Dock. Workspace and splitter layouts are saved with
the session and named layouts. There is no bottom Results dock.

Use **View > Full Screen** for a molecular viewing mode with panels and toolbars
hidden. Press **Escape** to restore the previous layout and window state.
This redesign updates the saved layout version, so older dock arrangements reset
once; project files and molecular data are unaffected.

The redesign is being delivered in separate steps: (1) main workspace and layout,
(2) molecular viewer context menus, (3) docking workflow, (4) Results & Analysis
workspace, and (5) visual design consistency. Each step is committed and pushed
individually. The current implementation completes steps 1 through 4.

Technical hardening (step 6) prevents overlapping docking runs, keeps Run/Cancel
actions in sync with background jobs, and retains every pose in CSV, Excel and
SDF exports when optional reference RMSD values are unavailable. The application
UI colour scheme is unchanged.

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
pytest          # preparation, docking, interactions, GUI workflows and background jobs
```

## Windows executable

```powershell
pip install pyinstaller
pyinstaller --noconfirm wavedock.spec     # -> dist\WaweDock\WaweDock.exe (+ _internal\)
```

For an installer, install Inno Setup 6 and compile the setup recipe after building
the executable:

```powershell
ISCC.exe installer/wavedock.iss          # -> dist/installer/WaweDock-Setup-0.1.0.exe
```

Setup installs per user, creates desktop/Start menu shortcuts, and registers
**Open with** during installation. Uninstall removes only this installation's
registration, preserving other applications and newer WaweDock installations.
Portable builds register automatically when first opened. No admin rights or
manual registration script is required by users.

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
