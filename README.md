# WaweDock

Open-source molecular docking workbench for Windows: protein preparation, ligand preparation,
AutoDock Vina docking, interactive 3D visualisation and result analysis.
Architecture and roadmap: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Start

* **Desktop shortcut** `WaweDock` (also `WaweDock.lnk` in this folder), or
* from PowerShell:

```powershell
cd "D:\Tools\WaweDock"
.\.venv\Scripts\Activate.ps1
python -m mds
```

For a source checkout, `WaweDock.cmd` also launches from its own folder without
activating the environment or relying on an old shortcut working directory.

Startup logs are created **before** importing the GUI and chemistry libraries.
If startup raises an exception, a Windows error dialog identifies the problem and
the diagnostic report at `%LOCALAPPDATA%\WaweDock\WaweDock\logs\startup-error.txt`.
The same folder contains `mds.log`, windowed-launch console output and native crash
traces. If application data is unwritable, diagnostics use `%TEMP%\WaweDock\logs`.

Run `.\.venv\Scripts\python.exe -m mds --startup-check` to verify the main window,
3D viewer and the availability report using isolated settings. A usable protein
workspace can pass with chemistry unavailable. `--startup-check-full` additionally
requires the actual RDKit/Meeko backend and exits nonzero when it cannot load.
Neither check repairs missing/blocked dependencies or changes Windows security policy.

Windows Smart App Control can block an unsigned dependency even when the Python
code is correct. The observed failure was RDKit's `RDKitRDInchiLib` DLL. Its policy
event was `VerifiedAndReputableDesktop` (Code Integrity event 3077). Microsoft
[does not provide per-app Smart App Control exceptions](https://support.microsoft.com/en-us/windows/security/threat-malware-protection/smart-app-control-frequently-asked-questions).
A supported distribution needs trusted signatures/reputation for its native
libraries; signing only the outer executable does not ensure that its dependent
DLLs will load. See [Microsoft's signing guidance](https://learn.microsoft.com/en-us/windows/apps/develop/smart-app-control/code-signing-for-smart-app-control).
Optional chemistry is now checked in a separate local process. If it is missing,
blocked, crashes or times out, the application opens its protein workspace and
clearly disables ligand preparation, receptor PDBQT preparation, docking and
chemistry-dependent analysis. Protein viewing, selection, cleanup, diagnostics and
independently available PDBFixer/OpenMM repair operations remain usable. Project
saves preserve ligand/run data that the unavailable backend could not load. The
current machine passes the core startup check and correctly fails the full-feature
check while its RDKit DLL remains blocked; security settings are unchanged.

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

## Molecular appearance defaults

New views use a **black viewport** and **sequence-spectrum proteins** (blue at the
low residue-number end, red at the high end). The rest of the app keeps its current
UI theme. Atom/HETATM deletion, cleanup, undo and redo preserve the current viewport
background, camera, representation and colour settings. The spectrum range is kept
across edits so deleting terminal residues does not recolour the remaining protein.
When editing switches from a visible prepared receptor to its source, its appearance
is carried over as well.

## Pose-by-pose review

In **Results & Analysis**, choose the crystal reference, check **Lock reference**,
and check **Hide protein**. Use **Fit view** once, then **Previous / Next** or the
pose table to browse. The reference keeps its coordinates and stays visible;
protein objects are hidden without deletion and return to their prior visibility
when unchecked. **Keep camera** is on by default and enforced while the reference
is locked. Docked pose colours and representation carry across pose changes.
Navigation follows the current table sort/filter. Fit View never realigns atoms.

The [one-job pose-review audit](docs/audits/2026-10-22-pose-review.md) records what
was exercised and the remaining RDKit-related integration-test limitation.

## Responsiveness

Box validation reuses unchanged protein/reference coordinates and ligand dimensions;
file rewrites and coordinate changes invalidate the cached values. Contact previews
run in the background. Rapid pose changes keep only the latest queued request, and
late results cannot replace the currently selected pose's contacts. Protein contact
geometry is cached, with spatial preprocessing and bounded distance blocks to avoid
large temporary arrays. These optimizations preserve the existing distance criteria.

Developers can reproduce the box/contact timings with
`python tools/benchmark_responsiveness.py` from a source checkout.

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

## Installation profiles and Windows releases

The base Python installation provides the protein workspace. Chemistry and repair
backends are explicit extras:

```powershell
python -m pip install -e ".[docking,repair,dev]"  # full development environment
python -m pip install -e ".[repair,dev]"          # protein tools and supported repairs
```

Installing an extra does not override Windows trust decisions. When optional
chemistry cannot load, the existing window shows which tools are unavailable;
there are no simulated docking or repair results. Restart after changing installed
backends. Optional-backend tests are reported as skipped when their native library
cannot load; that is not a passing full scientific validation.

For a local development executable (not a public installer):

```powershell
$env:WAVEDOCK_BUILD_EDITION = "Protein"   # or Full
python -m PyInstaller --noconfirm wavedock.spec
```

**Protein** excludes RDKit, Meeko and associated docking-only packages from the
bundle and registers only PDB/mmCIF formats. **Full** includes docking dependencies
and enables their tools only after a successful capability check. Both reuse the
same UI, selection, viewer, project format and optional repair backend. Qt's Windows
baseline is Windows 10 1809 or later / Windows 11; these installer builds target
x64. See [Qt platform requirements](https://doc.qt.io/qt-6.10/supported-platforms.html).

Public installer builds must use a valid trusted code-signing certificate with a
private key in the current user's Windows certificate store, Windows SDK SignTool
and Inno Setup 6 on PATH:

```powershell
.\scripts\build_windows_release.ps1 -CertificateThumbprint YOUR_CERTIFICATE_THUMBPRINT -Edition Full
.\scripts\build_windows_release.ps1 -CertificateThumbprint YOUR_CERTIFICATE_THUMBPRINT -Edition Protein
```

Full builds first require a working chemistry backend and the real scientific/GUI
docking tests on the approved build host. The script then signs native EXE/DLL/PYD files that do not already have valid signatures,
verifies every native signature, launches the frozen executable outside the checkout
with the appropriate startup check, and then builds a signed installer/uninstaller.
A Full build cannot pass on protein-only fallback. The installer recipe refuses
unsigned direct builds. Signing credentials are not stored in the repository and
no certificate is created or trusted automatically. A self-signed certificate is
not a substitute for a certificate Windows trusts.

Setup installs per user and creates desktop/Start menu shortcuts and file
associations. Output names include the edition. Open Babel remains separate for
licensing reasons, so its interaction-analysis feature is not bundled.

**Verification boundary:** source/native protein-workspace behavior is tested on
this machine with Smart App Control unchanged. No production-signed installer or
clean-machine compatibility audit has been produced in this increment: a trusted
signing certificate is not installed here, and local PowerShell policy prevents
release-script execution. The script syntax is verified; its signing/build stages
require an approved build environment and were not executed here. No execution-policy
bypass is applied. The workflow refuses to build a public installer without signing
credentials. Before publishing, run the signed installers on
clean supported Windows systems with Smart App Control enabled; this code does
not claim compatibility with every architecture, driver or organisation policy.

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

## Protein health and repair

Open **Docking Setup > Receptor Preparation > Protein Health & Repair** for an asynchronous, read-only health check, issue navigation in the existing viewer and JSON report export. Checks include canonical heavy-atom completeness, within-segment backbone connectivity, broad geometry/overlap screens, alternate conformers and sequence-backed unmodelled regions. Existing structure colours, selection and preparation remain available.

**Repair Protein** now proposes selectable conformer choices, supported missing-atom rebuilding and opt-in hydrogen addition. Build a separate candidate, compare it with the original, then accept it after review. Existing coordinates and non-protein components are retained; Undo/Redo and PDB/mmCIF export with provenance are available. Install backends with `python -m pip install -e ".[repair]"`.

Sequence-backed internal short loops (up to 8 residues) can now be modelled locally with PDBFixer and restrained Amber14/OpenMM refinement. They are separate, opt-in computational hypotheses with independent geometry gates and fixed original coordinates; no MODELLER setup is needed. General refinement and statistical Ramachandran/rotamer validation are not provided. Read [coverage, repair workflow and scientific limits](docs/protein-health.md) before interpreting findings.
