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
    t0 = time.monotonic()
    rec = prepare_receptor(Path(structure), out_dir / "receptor")
    hetero = next((h for h in rec.hetero if h.name == resname.upper() and h.extracted_to), None)
    if hetero is None:
        raise ValueError(f"Ligand {resname} not found in {Path(structure).name}")
    template = Chem.MolFromSmiles(smiles)
    reference = reference_from_pdb_ligand(hetero.extracted_to, Chem.RemoveHs(template))
    ligand = prepare_ligand(entry_from_smiles(smiles, resname), seed=seed, ionize=ionize)
    box = box_from_ligand_file(hetero.extracted_to)
    engine = VinaEngine()
    outcome = run_docking(engine, rec.pdbqt, [ligand], box,
                          DockingParameters(exhaustiveness=exhaustiveness, seed=seed), out_dir / "run",
                          reference=reference)[0]
    if outcome.status != "ok":
        raise RuntimeError(outcome.message)
    top_rmsd = outcome.reference_rmsd[0]
    result = {
        "structure": Path(structure).name, "ligand": resname, "engine": engine.version(),
        "top_score_kcal_per_mol": outcome.poses[0].score, "top_pose_rmsd_A": top_rmsd,
        "best_rmsd_any_pose_A": min(r for r in outcome.reference_rmsd if r is not None),
        "success_top_pose": bool(top_rmsd is not None and top_rmsd < 2.0),
        "exhaustiveness": exhaustiveness, "seed": seed, "elapsed_s": round(time.monotonic() - t0, 1),
    }
    (out_dir / "validation.json").write_text(json.dumps(result, indent=2))
    return result
