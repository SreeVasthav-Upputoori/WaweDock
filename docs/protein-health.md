# Protein Health & Repair

Open **Docking Setup > Receptor Preparation > Protein Health & Repair**. The existing preparation controls and viewer remain available. The new dock can be collapsed, tabbed, hidden or reopened from Window > Panels.

**Check Protein** screens all chains in the loaded source structure's first model in a background job. **Validate Structure** repeats these same supported screens against the current source coordinates. Neither action prepares a receptor, changes coordinates, establishes docking readiness or certifies that the structure is healthy.

The issue list contains severity and location; selecting an issue displays the finding, recommended correction and review status, selects its exact available atoms and focuses the existing viewer. Unmodelled regions focus their observed flanking residues when available. Findings with no coordinates cannot be selected in 3D. Right-click an existing editable selection and choose **Inspect Protein Issue** to locate a matching diagnosed issue.

Optional small crosses mark errors in red, warnings in orange and notes in blue. Molecular selections remain yellow. These markers do not replace or recolour the protein representation. At most 250 residue markers are drawn, prioritising errors and warnings; the issue list is authoritative. Markers follow structure/component visibility and can be disabled. Issue counts are capped at 5,000 for pathological input, with truncation disclosed in the report.

## Checks and boundaries

- Canonical amino-acid backbone and side-chain heavy-atom completeness, checked per alternate conformer with shared blank-altloc atoms. OXT is optional; present OXT is included in bonded-neighbour exclusions.
- Peptide C-N connectivity within a continuous segment. Explicit PDB TER boundaries, separate chains and mmCIF subchain boundaries are respected. Numbering jumps alone never establish missing sequence.
- Unmodelled sequence regions require a supplied entity sequence and a consistent mapping of observed residues. Sequence-only alignments must be unique; ambiguous repeated sequences are not classified. Terminal unmodelled sequence is informational. No claim of experimental disorder is made without experimental metadata.
- Broad backbone N-CA, CA-C and C-O length screens and N-CA-C angle screening. These are conservative flags, not restraint-library Z-scores.
- Severe heavy-atom overlaps (van der Waals overlap greater than 0.8 Angstrom), using a spatial grid and standard amino-acid topology exclusions through two bonds. Incompatible altloc pairs and plausible S-S disulfide pairs are excluded. Metals, symmetry mates and hydrogen clashes are outside this screen. Nonstandard covalent links require manual review.
- Alternate conformations, protein hydrogen presence and nonstandard/component inventory. Finding some hydrogens does not establish hydrogen completeness or correct protonation.

Gemmi supplies parsing, entity information and sequence alignment; see its [structure and sequence documentation](https://gemmi.readthedocs.io/en/latest/mol.html). The existing exact atom identity bridge is reused for navigation. Native parser traversal and original atom indices are retained even when an additional TER-aware parse is required to identify segment boundaries.

**Not performed:** Ramachandran/rotamer validation, full ligand stereochemistry, charge/force-field checks, pKa prediction or reconstruction. All findings require review. Input coordinates may themselves be experimental or computational; this module cannot establish their provenance. The diagnostic actions do not reconstruct coordinates. Supported, explicitly selected repairs are described below.

## Export and lifecycle

Export the JSON report from the dock. It includes source SHA-256, backend/version, findings, exact atom selectors, coverage and limitations. Structure edits, Undo/Redo, new imports and project changes invalidate the report and cancel outdated work. Cancel does not change coordinates. The input is never written by the diagnostics module.

## Supported automatic and manual repairs

Install the optional repair backends into the application environment:

```powershell
python -m pip install -e ".[repair]"
```

The tested development environment uses official PDBFixer 1.12.0, OpenMM 8.6.1 and Biopython 1.88. Without these packages, conformer selection remains functional and the proposal explains why template reconstruction/hydrogens are unavailable. The PyInstaller recipe collects installed repair backends, templates, metadata and OpenMM native libraries. A new executable has not been built or audited by this increment; build it from an environment with the repair extra installed.

1. Click **Repair Protein** to generate a proposal. **All protein residues** covers the first model's canonical protein residues; **Viewer selection** uses the existing viewer's exact selected residues. Change the scope/pH and click Repair Protein again to refresh the proposal. The right-click **Repair Missing Atoms...** action opens the same manual workflow.
2. Review the checkboxes. Alternate conformers are separate choices, with the highest mean occupancy recommended; choose only one per residue. Missing side-chain atoms and backbone O can be reconstructed when N/CA/C anchors are present. Hydrogen addition is optional and initially unchecked. No coordinates change at this point.
3. Click **Build selected repairs**. PDBFixer rebuilds only the selected atom names with seed 42; terminal oxygens, loops and residue replacements are disabled. Existing coordinates are fixed. All original waters, cofactors, metals, ligands, nonstandard residues, other models and PDB CONECT records are preserved in the output. Alternate conformers with explicit connectivity are not modified automatically.
4. Inspect **Before: original / After: reconstructed**. The original remains loaded; the candidate is a read-only overlay using the same camera and appearance. Reconstructed atoms use thin cyan sticks; original atoms in the same bonded residue are shown as element-coloured sticks to expose their anchors in cartoon/ribbon view. Only parsed bonds are displayed, subject to the existing component-visibility settings. Ordinary selection stays yellow. The issue report and reconstruction provenance remain available.
5. **Use repaired structure** activates only an eligible candidate after review, invalidates previous receptor preparation and adds a revision to the existing Undo/Redo history. Blocked candidates cannot be activated or exported as repaired structures. Original files are immutable. Undo restores the previous source; redo restores the candidate and its reconstructed-atom markings.

Hydrogen placement uses Amber14/ff14SB force-field templates, reproducible initialisation (seed 42) and default residue pKas at the selected pH. It is available only for a complete canonical protein without existing protein hydrogens, unresolved internal gaps or missing required terminal OXT. Missing heavy atoms and alternate conformers must be resolved by the selected changes before hydrogens are added. Explicit input disulfides are carried into the backend and their geometry is checked. Non-protein components are excluded from protonation but retained in the output and included in reconstructed-atom contact screening. This does not assign docking charges, predict pKa with PROPKA, or settle metal-site/binding-site protonation. Re-run receptor preparation with appropriate protonation/charge settings before docking.

Native repair computation runs in an isolated worker process with two CPU threads. The interface and molecular viewer remain responsive; cancellation stops the native worker. A stale source/proposal/candidate cannot replace the active structure. An abandoned candidate may remain in the project's prepared/protein_repair folder, but is never activated automatically.

### Geometry gates and uncertainty

The real merged coordinates are re-read and checked, including unchanged identities/coordinates of surviving atoms, selected-residue completeness, reconstructed topology bonds and broad bond-length/angle bounds, severe overlaps with the full retained structure (including hydrogens), and new screening findings in modified regions. Plausible donor-H...acceptor contacts are distinguished from severe overlaps using donor identity, distance and angle. These conservative screens do not replace full stereochemical validation. Failed geometry gates block acceptance even if PDBFixer/minimisation completed. Existing abnormalities elsewhere may remain; accepted candidates are labelled **Review required**, never universally healthy or docking-ready.

### Provenance and exports

**Export repaired structure** writes PDB or mmCIF plus a `.repair.json` sidecar containing backend versions, source/candidate/export SHA-256, selected corrections, per-change status, added/removed atom identities, inherited reconstructed atoms, before/after reports, blockers and warnings. The structure also carries explicit atom-origin annotations (PDB REMARK 950 or a custom mmCIF category). Existing mmCIF metadata categories, including cofactor bond dictionaries, are retained during mmCIF repair. PDB CONECT bonds are represented as standard struct_conn records on mmCIF export; mmCIF component/explicit bonds are represented as CONECT on PDB export. Re-import recognises these annotations; later repair revisions retain prior reconstructed atoms' provenance. Atom origin is never guessed from B factors. PDB export uses standard coordinate precision; prefer mmCIF when preserving additional decimal precision matters.

Loop reconstruction, missing-residue modelling, rebuilding N/CA/C anchors, general selected-region refinement and Ramachandran/rotamer validation remain for the advanced increment. There are no placeholder refinement/reconstruction buttons. See [PDBFixer's documented operations](https://github.com/openmm/pdbfixer/blob/master/Manual.html) and [OpenMM's hydrogen-placement interface](https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html).

## Validation record (2026-10-22)

The new diagnostic tests cover unchanged input files, separate chains, terminal OXT, explicit TER boundaries, insertion codes, per-conformer missing atoms, unambiguous/ambiguous sequence evidence, geometric flags, report export and cancellation. Real MainWindow/QtWebEngine tests cover asynchronous completion, issue-to-atom navigation, exact representation styles, background/atom-count preservation, marker visibility and invalidation of in-flight results after source changes. A native application preview was also inspected.

The broader native Windows regression run passed 81 checks and failed two existing tests. Both failures were reproduced against the unchanged preceding commit (`2c17e5b`) in an isolated archived copy:

- `test_molecular_fullscreen_restores_panels_and_maximised_window`: immediate `isMaximized()` assertion after restoring fullscreen.
- `test_results_workspace_displays_selected_pose_and_keeps_viewer_state`: expects the former `pose: ligand B #1` object name instead of the existing stable pose-review object.

These baseline failures were not changed or hidden by this increment.

## Supported-repair validation record (2026-10-22)

The repair increment is covered by 60 focused scientific/native-application checks. These included real PDBFixer missing-atom reconstruction and Amber14 hydrogen placement on 3PTB, preserving every existing atom's coordinates/occupancy/B factor/serial and ligand/water/metal content; preserving CONECT and explicit disulfides; conformer selection; immutable source files; blocked geometry and stale candidate checks; PDB/mmCIF export/re-import; and provenance inherited by later repairs. Native MainWindow/QtWebEngine tests exercised the real worker process with an active GUI heartbeat timer, cancellation, exact manual residue scope, staged comparisons, acceptance, Undo/Redo, source colours/background and narrow panel layouts. Final comparison lifecycle tests also verify scene restoration when the health panel is hidden. Additional export checks verify actual (non-empty) CONECT preservation, both bond-conversion directions and source mmCIF metadata retention. The packaged executable has not been rebuilt or audited in this increment.


## Repair audit correction (2026-10-23)

The current backend completes atoms in existing residues; it does **not** reconstruct loops. The panel displays this limitation before building, reports zero loops rebuilt, and shows the first rejection reasons beside the candidate rather than relying on warning counts.

An audit of the user's 1XBB candidate reproduced 72 added atoms and severe side-chain overlaps, including GLU A:600 / TYR A:629 at approximately 1.42 Angstrom. The candidate remains blocked. Its unresolved internal regions A:406?410 and A:529?532 were not reconstructed. The original file contains sequence and missing-residue records, but earlier cleanup used a minimal PDB writer that discarded sequence evidence. Subsequent edits now retain PDB sequence metadata. Previously edited files are not silently patched from another file; use the original import to retain its evidence.

Both heavy-atom and hydrogen backend snapshots now separate geometrically disconnected backbone segments rather than constructing artificial peptide bonds across them. Completing a side chain no longer makes an unchanged backbone break appear to be a new validation failure: comparisons use its C/N endpoint identities and coordinates. Existing breaks remain reported. These changes do not relax the clash or geometry gates, provide a loop-modelling backend, or establish docking readiness.
