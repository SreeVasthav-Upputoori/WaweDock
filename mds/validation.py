"""Redocking validation: can the pipeline reproduce an experimentally determined pose?

Usage (also works with the packaged executable)::

    WaweDock.exe --validate <structure.pdb> <LIGAND_RESNAME> <ligand SMILES> <output_dir>

The co-crystal ligand is rebuilt in 3D from its SMILES (its crystal coordinates are
NOT used as input), docked into the prepared receptor, and the heavy-atom RMSD of the
top-ranked pose to the crystal pose is reported.  RMSD < 2.0 A is the conventional
success criterion.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from rdkit import Chem

from mds.chem.box import box_from_ligand_file
from mds.chem.ligand import entry_from_smiles, prepare_ligand
from mds.chem.poses import reference_from_pdb_ligand
from mds.chem.receptor import prepare_receptor
from mds.chem.workflow import DockingParameters, run_docking
from mds.engines import VinaEngine


def redock(structure: Path, resname: str, smiles: str, out_dir: Path, *, exhaustiveness: int = 16,
           ionize: str = "simple", seed: int = 7) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    structure = Path(structure)
    rec = prepare_receptor(structure, out_dir / "receptor")
    resname = resname.upper()
    hetero = next((h for h in rec.hetero if h.name == resname.upper() and h.extracted_to), None)
    if hetero is None:
        raise ValueError(f"Ligand {resname} not found in {structure.name}")
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError("The co-crystal ligand SMILES could not be parsed")
    template = Chem.RemoveHs(template)
    if template.GetNumHeavyAtoms() == 0:
        raise ValueError("The co-crystal ligand has no heavy atoms")
    reference = reference_from_pdb_ligand(hetero.extracted_to, template)
    ligand = prepare_ligand(entry_from_smiles(smiles, resname), seed=seed, ionize=ionize)
    if ligand.mol is None:
        raise RuntimeError(f"Ligand {resname} could not be prepared: {ligand.report.render()}")
    box = box_from_ligand_file(hetero.extracted_to)
    engine = VinaEngine()
    outcome = run_docking(engine, rec.pdbqt, [ligand], box,
                          DockingParameters(exhaustiveness=exhaustiveness, seed=seed), out_dir / "run",
                          reference=reference)[0]
    if outcome.status != "ok":
        raise RuntimeError(outcome.message)
    rmsds = [rmsd for rmsd in outcome.reference_rmsd if rmsd is not None]
    top_rmsd = outcome.reference_rmsd[0] if outcome.reference_rmsd else None
    best_rmsd = min(rmsds) if rmsds else None
    result = {
        "structure": Path(structure).name, "ligand": resname, "engine": engine.version(),
        "top_score_kcal_per_mol": outcome.poses[0].score, "top_pose_rmsd_A": top_rmsd,
        "best_rmsd_any_pose_A": best_rmsd,
        "success_top_pose": bool(top_rmsd is not None and top_rmsd < 2.0),
        "exhaustiveness": exhaustiveness, "seed": seed, "elapsed_s": round(time.monotonic() - t0, 1),
    }
    (out_dir / "validation.json").write_text(json.dumps(result, indent=2))
    return result
