# Pose review audit - 22 October 2026

One real AutoDock Vina docking job completed using copied, previously prepared
1XBB receptor and caffeine ligand inputs. The extracted STI crystal ligand was
kept as the visual reference. Original project files were left unchanged.

This audits the review workflow; caffeine and STI are different compounds, so
this is **not** a redocking accuracy or reference-RMSD benchmark. Scores below
are Vina predictions in kcal/mol, not experimental affinities.

Parameters: exhaustiveness 8, six modes requested, energy range 5 kcal/mol,
CPU 2, seed 42. All six returned poses were loaded into the native Windows
QtWebEngine/3Dmol viewer through the shared pose-review session and controls.

| Pose | Predicted score (kcal/mol) |
| --- | --- |
| 1 | -5.319 |
| 2 | -5.116 |
| 3 | -5.077 |
| 4 | -5.038 |
| 5 | -4.983 |
| 6 | -4.965 |

## Verified

- All six pose coordinate sets are distinct and can be visited with Next/Previous.
- Reference atom coordinates and camera view stay unchanged on pose replacement.
- Lock Reference prevents hiding the selected reference and hides other references.
- Hide Protein changes visibility only; restoring it retains the model.
- Pose colours and representation persist across replacements.
- Fit View frames the locked reference and current pose together.
- Clearing the current pose leaves the protein and reference intact.

## Gaps fixed

Added explicit reference-lock, keep-camera, hide-protein and navigation controls.
Pose replacement now uses one stable scene object with preserved appearance.
Results completion preserves that object's Project Explorer entry. Hidden proteins
are excluded from automatic nearby-residue styling while poses are browsed.

Thirty independent regression checks passed, including native Windows viewer tests,
reference-lock state, reversible visibility, job handling and numerical geometry.

## Validation limit

Windows App Control still rejects RDKit's `rdmolfiles` DLL. Fresh preparation and
the complete MainWindow/controller workflow could not be executed in this
session. The real engine job used existing prepared inputs; viewer checks used
the same review-session class and control widget wired into MainWindow. Full
workspace regression tests were added for execution once RDKit is available.

## Local artifacts

`test/pose-review-audit/` contains copied inputs, their SHA-256 fingerprints in
`audit-job.json`, the Vina output/log, all six pose files, `viewer-audit.json`,
and screenshots under `screenshots/`. These scratch artifacts are not committed.
Recheck that run with `python tools/audit_pose_review.py test/pose-review-audit`.
