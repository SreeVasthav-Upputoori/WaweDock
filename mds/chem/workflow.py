"""Docking workflow: ligands x receptor -> scored, exported, reproducible run folder.

Run folder layout (``runs/<run-id>/``)::

    run.json            parameters, software versions, seeds, commands, timings
    receptor/           receptor PDBQT used
    <ligand>/           ligand PDBQT, vina.log, poses.pdbqt, poses.sdf
    results.csv         one row per pose
"""

from __future__ import annotations

import csv
import json
import logging
import platform
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

from rdkit import Chem

from mds import __version__
from mds.chem.box import DockingBox
from mds.chem.ligand import PreparedLigand
from mds.chem.poses import poses_to_mols, symmetric_rmsd, write_poses_sdf
from mds.engines.base import (CancelFn, DockingCancelled, DockingEngine, DockingError,
                              DockingRequest, DockingResult, Pose, ProgressFn)

log = logging.getLogger("mds.chem.workflow")

SCORE_DISCLAIMER = ("Scores are predicted by an empirical scoring function in kcal/mol. They are "
                    "not measured binding affinities and are not evidence of biological activity.")


@dataclass
class DockingParameters:
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0
    cpu: int = 0
    seed: int = 42
    scoring: str = "vina"


@dataclass
class LigandOutcome:
    name: str
    status: str                          # "ok" | "failed" | "cancelled"
    message: str = ""
    smiles: str = ""
    poses: list[Pose] = field(default_factory=list)
    mols: list[Chem.Mol] = field(default_factory=list)
    reference_rmsd: list[float | None] = field(default_factory=list)
    result: DockingResult | None = None
    folder: Path | None = None

    @property
    def best_score(self) -> float | None:
        return self.poses[0].score if self.poses else None


def run_docking(
    engine: DockingEngine,
    receptor_pdbqt: Path,
    ligands: list[PreparedLigand],
    box: DockingBox,
    params: DockingParameters,
    run_dir: Path,
    *,
    reference: Chem.Mol | None = None,
    progress: ProgressFn | None = None,
    cancelled: CancelFn | None = None,
    receptor_name: str = "",
) -> list[LigandOutcome]:
    """Dock every prepared ligand; failures are recorded per ligand, not fatal for the batch."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "receptor").mkdir(exist_ok=True)
    rec_copy = run_dir / "receptor" / "receptor.pdbqt"
    rec_copy.write_text(Path(receptor_pdbqt).read_text())

    started = datetime.now().isoformat(timespec="seconds")
    outcomes: list[LigandOutcome] = []
    commands: dict[str, list[list[str]]] = {}
    total_jobs = sum(max(1, len(l.conformer_pdbqts)) for l in ligands)
    done = 0

    for lig in ligands:
        folder = run_dir / _safe(lig.name)
        folder.mkdir(exist_ok=True)
        outcome = LigandOutcome(lig.name, "ok", smiles=lig.smiles, folder=folder)
        all_poses: list[tuple[Pose, int]] = []
        last_result: DockingResult | None = None
        try:
            for ci, text in enumerate(lig.conformer_pdbqts or [lig.pdbqt]):
                if cancelled and cancelled():
                    raise DockingCancelled()
                lig_file = folder / f"ligand_conf{ci + 1}.pdbqt"
                lig_file.write_text(text)
                request = DockingRequest(rec_copy, lig_file, box, params.exhaustiveness,
                                         params.num_modes, params.energy_range, params.cpu,
                                         params.seed, params.scoring)

                def sub_progress(pct: float, msg: str, _d=done) -> None:
                    if progress:
                        progress((_d + pct / 100) / total_jobs * 100, f"{lig.name}: {msg}")

                result = engine.run(request, folder / f"conf{ci + 1}", sub_progress, cancelled)
                commands.setdefault(lig.name, []).append(result.command)
                last_result = result
                all_poses += [(p, ci) for p in result.poses]
                done += 1
            all_poses.sort(key=lambda t: t[0].score)
            ordered: list[Pose] = []
            conf_text = {i: t for i, t in enumerate(lig.conformer_pdbqts or [lig.pdbqt])}
            for rank, (pose, ci) in enumerate(all_poses[: max(params.num_modes, 1)], start=1):
                ordered.append(Pose(rank, pose.score, pose.rmsd_lb, pose.rmsd_ub, pose.pdbqt, ci))
            outcome.poses, outcome.result = ordered, last_result
            # Pose molecules (bond orders from the ligand's own REMARK SMILES)
            mols = []
            for pose in ordered:
                mols.append(poses_to_mols([pose], conf_text[pose.conformer])[0])
            outcome.mols = mols
            for m in mols:
                outcome.reference_rmsd.append(_safe_rmsd(m, reference) if reference is not None else None)
            _write_ligand_outputs(folder, outcome)
        except DockingCancelled:
            outcome.status, outcome.message = "cancelled", "Cancelled by user"
            outcomes.append(outcome)
            log.warning("Docking cancelled during %s", lig.name)
            break
        except (DockingError, Exception) as exc:  # noqa: BLE001 — record, continue batch
            outcome.status, outcome.message = "failed", str(exc)
            log.error("Docking failed for %s: %s", lig.name, exc)
        outcomes.append(outcome)

    _write_run_files(run_dir, engine, receptor_name, box, params, ligands, outcomes, commands,
                     started)
    return outcomes


def _safe_rmsd(mol: Chem.Mol, reference: Chem.Mol) -> float | None:
    try:
        return round(symmetric_rmsd(mol, reference), 3)
    except ValueError as exc:
        log.warning("Reference RMSD unavailable: %s", exc)
        return None


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name) or "ligand"


def _write_ligand_outputs(folder: Path, outcome: LigandOutcome) -> None:
    (folder / "poses.pdbqt").write_text("".join(f"MODEL {p.rank}\n{p.pdbqt}ENDMDL\n"
                                               for p in outcome.poses))
    props = []
    for i, p in enumerate(outcome.poses):
        r = outcome.reference_rmsd[i] if i < len(outcome.reference_rmsd) else None
        props.append({"Name": outcome.name, "Rank": p.rank, "Score_kcal_per_mol": p.score,
                      "RMSD_lb": p.rmsd_lb, "RMSD_ub": p.rmsd_ub,
                      "RMSD_to_reference": r if r is not None else "",
                      "Engine": outcome.result.engine if outcome.result else "",
                      "Note": "Predicted score, not a measured affinity"})
    write_poses_sdf(folder / "poses.sdf", outcome.mols, props)


def results_rows(outcomes: list[LigandOutcome]) -> list[dict]:
    rows = []
    for o in outcomes:
        if o.status != "ok":
            rows.append({"Ligand": o.name, "Pose": "", "Score (kcal/mol)": "", "Status": o.status,
                         "Message": o.message})
            continue
        for i, p in enumerate(o.poses):
            r = o.reference_rmsd[i] if i < len(o.reference_rmsd) else None
            rows.append({"Ligand": o.name, "Pose": p.rank, "Score (kcal/mol)": p.score,
                         "RMSD lb": p.rmsd_lb, "RMSD ub": p.rmsd_ub,
                         "RMSD to reference (A)": r if r is not None else "",
                         "Conformer": p.conformer + 1, "SMILES": o.smiles, "Status": "ok"})
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    cols: list[str] = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with Path(path).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def _write_run_files(run_dir: Path, engine: DockingEngine, receptor_name: str, box: DockingBox,
                     params: DockingParameters, ligands: list[PreparedLigand],
                     outcomes: list[LigandOutcome], commands: dict, started: str) -> None:
    import meeko
    import rdkit

    write_csv(run_dir / "results.csv", results_rows(outcomes))
    meta = {
        "application": f"WaweDock {__version__}",
        "started": started, "finished": datetime.now().isoformat(timespec="seconds"),
        "engine": engine.name, "engine_version": engine.version(),
        "receptor": receptor_name, "box": box.as_dict(), "parameters": asdict(params),
        "ligands": [{"name": l.name, "smiles": l.smiles, "properties": l.properties}
                    for l in ligands],
        "commands": commands,
        "software": {"python": platform.python_version(), "rdkit": rdkit.__version__,
                     "meeko": meeko.__version__, "platform": platform.platform()},
        "outcomes": {o.name: {"status": o.status, "message": o.message,
                              "best_score_kcal_per_mol": o.best_score} for o in outcomes},
        "disclaimer": SCORE_DISCLAIMER,
    }
    (run_dir / "run.json").write_text(json.dumps(meta, indent=2))

