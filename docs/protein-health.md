# Protein Health & Repair: diagnostic increment

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

**Not performed:** Ramachandran/rotamer validation, full ligand stereochemistry, charge/force-field checks, pKa prediction or reconstruction. All findings require review. Input coordinates may themselves be experimental or computational; this module cannot establish their provenance. No coordinates are reconstructed by this increment.

## Export and lifecycle

Export the JSON report from the dock. It includes source SHA-256, backend/version, findings, exact atom selectors, coverage and limitations. Structure edits, Undo/Redo, new imports and project changes invalidate the report and cancel outdated work. Cancel does not change coordinates. The input is never written by the diagnostics module.

## Following increments

Supported repairs will be added separately, with an explicit selectable change proposal, preserved cofactors/nonstandard components, provenance, before/after checks and the existing Undo history. PDBFixer, OpenMM and Biopython are not installed in the current development environment; no inactive repair buttons or simulated repairs are included. [PDBFixer's documented operations](https://github.com/openmm/pdbfixer/blob/master/Manual.html) provide a reference for evaluating the next repair adapter.

Sequence-guided loop reconstruction and restrained refinement require a tested modelling backend plus geometry validation. They must not be represented as successful merely because a backend or minimiser finished. Repaired-structure export belongs to that repair increment; this read-only increment exports validation reports only.

## Validation record (2026-10-22)

The new diagnostic tests cover unchanged input files, separate chains, terminal OXT, explicit TER boundaries, insertion codes, per-conformer missing atoms, unambiguous/ambiguous sequence evidence, geometric flags, report export and cancellation. Real MainWindow/QtWebEngine tests cover asynchronous completion, issue-to-atom navigation, exact representation styles, background/atom-count preservation, marker visibility and invalidation of in-flight results after source changes. A native application preview was also inspected.

The broader native Windows regression run passed 81 checks and failed two existing tests. Both failures were reproduced against the unchanged preceding commit (`2c17e5b`) in an isolated archived copy:

- `test_molecular_fullscreen_restores_panels_and_maximised_window`: immediate `isMaximized()` assertion after restoring fullscreen.
- `test_results_workspace_displays_selected_pose_and_keeps_viewer_state`: expects the former `pose: ligand B #1` object name instead of the existing stable pose-review object.

These baseline failures were not changed or hidden by this increment.
