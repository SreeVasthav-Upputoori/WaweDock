"""Background jobs wrapping the scientific core (no widgets here)."""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path



from mds.chem.ligand import LigandEntry, PreparedLigand, prepare_ligand
from mds.chem.interactions import InteractionResult, compute_interactions, protein_to_rdkit
from mds.chem.receptor import ReceptorOptions, ReceptorResult, inspect_receptor, prepare_receptor
from mds.chem.report import PreparationError, Report
from mds.chem.similarity import SimilarityResult, cluster, compare_to_reference
from mds.chem.workflow import DockingParameters, LigandOutcome, run_docking
from mds.chem.box import DockingBox
from mds.engines import VinaEngine
from mds.engines.base import DockingCancelled
from mds.jobs.runner import Job, JobContext


class FetchPdbJob(Job):
    """Download a structure from the RCSB Protein Data Bank."""

    def __init__(self, pdb_id: str, dest_dir: Path) -> None:
        self.pdb_id = pdb_id.strip().upper()
        self.dest_dir = Path(dest_dir)
        self.name = f"Fetch {self.pdb_id} from RCSB PDB"

    def run(self, ctx: JobContext) -> Path:
        if not (len(self.pdb_id) == 4 and self.pdb_id.isalnum()):
            raise ValueError(f"'{self.pdb_id}' is not a valid 4-character PDB ID")
        ctx.progress(10, "Contacting files.rcsb.org")
        url = f"https://files.rcsb.org/download/{self.pdb_id}.pdb"
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310 — fixed https host
                data = resp.read()
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"Download failed for {self.pdb_id}: {exc}") from exc
        ctx.check_cancelled()
        self.dest_dir.mkdir(parents=True, exist_ok=True)
        path = self.dest_dir / f"{self.pdb_id}.pdb"
        path.write_bytes(data)
        ctx.progress(100, "Downloaded")
        return path


class ContactsJob(Job):
    def __init__(self, source: Path, pose) -> None:
        self.source, self.pose = Path(source), pose
        self.name = "Pose contact preview"

    def run(self, ctx: JobContext):
        from mds.chem.contacts import cached_protein_atoms, pose_contacts
        ctx.check_cancelled()
        protein = cached_protein_atoms(self.source)
        ctx.check_cancelled()
        result = pose_contacts(self.pose, protein)
        ctx.check_cancelled()
        return result


class InspectReceptorJob(Job):
    def __init__(self, source: Path, options: ReceptorOptions) -> None:
        self.source, self.options = Path(source), options
        self.name = f"Check receptor {self.source.name}"

    def run(self, ctx: JobContext) -> Report:
        ctx.check_cancelled()
        ctx.progress(10, "Checking receptor structure and preparation options")
        result = inspect_receptor(self.source, self.options)
        ctx.check_cancelled()
        ctx.progress(100, "Receptor check complete")
        return result


class PrepareReceptorJob(Job):
    def __init__(self, source: Path, out_dir: Path, options: ReceptorOptions) -> None:
        self.source, self.out_dir, self.options = Path(source), Path(out_dir), options
        self.name = f"Prepare receptor {self.source.name}"

    def run(self, ctx: JobContext) -> ReceptorResult:
        ctx.progress(5, "Cleaning structure and adding hydrogens")
        result = prepare_receptor(self.source, self.out_dir, self.options)
        ctx.progress(100, "Receptor ready")
        return result


@dataclass
class LigandPrepOutcome:
    prepared: list[PreparedLigand]
    failures: list[tuple[str, Report]]


class PrepareLigandsJob(Job):
    def __init__(self, entries: list[LigandEntry], out_dir: Path, *, seed: int, n_conformers: int,
                 ionize: str) -> None:
        self.entries, self.out_dir = entries, Path(out_dir)
        self.seed, self.n_conformers, self.ionize = seed, n_conformers, ionize
        self.name = f"Prepare {len(entries)} ligand(s)"

    def run(self, ctx: JobContext) -> LigandPrepOutcome:
        prepared: list[PreparedLigand] = []
        failures: list[tuple[str, Report]] = []
        for i, entry in enumerate(self.entries, start=1):
            ctx.check_cancelled()
            ctx.progress(100 * (i - 1) / len(self.entries), f"{entry.name} ({i}/{len(self.entries)})")
            try:
                lig = prepare_ligand(entry, seed=self.seed, n_conformers=self.n_conformers,
                                     ionize=self.ionize)
                lig.write(self.out_dir)
                prepared.append(lig)
            except PreparationError as exc:
                failures.append((entry.name, exc.report if exc.report.entries else _single(str(exc))))
        ctx.progress(100, "Ligands prepared")
        return LigandPrepOutcome(prepared, failures)


def _single(message: str) -> Report:
    rep = Report()
    rep.error(message)
    return rep


class DockingJob(Job):
    def __init__(self, receptor_pdbqt: Path, ligands: list[PreparedLigand], box: DockingBox,
                 params: DockingParameters, run_dir: Path, reference: Chem.Mol | None,
                 receptor_name: str) -> None:
        self.args = (receptor_pdbqt, ligands, box, params, run_dir)
        self.reference, self.receptor_name = reference, receptor_name
        self.name = f"AutoDock Vina: {len(ligands)} ligand(s)"

    def run(self, ctx: JobContext) -> list[LigandOutcome]:
        receptor_pdbqt, ligands, box, params, run_dir = self.args
        try:
            return run_docking(VinaEngine(), receptor_pdbqt, ligands, box, params, run_dir,
                               reference=self.reference, progress=ctx.progress,
                               cancelled=lambda: ctx.cancelled, receptor_name=self.receptor_name)
        except DockingCancelled:
            ctx.check_cancelled()
            raise


class InteractionJob(Job):
    """ProLIF interaction fingerprint + 2D diagram for one docked pose.

    Rebuilding the receptor as an RDKit molecule (Open Babel bond perception) is the slow part
    of this job (roughly a second for a few thousand atoms); the caller passes a cached result
    in when analysing several poses of the same receptor so it is only built once.
    """

    def __init__(self, receptor_pdb: Path, pose, *, cached_protein=None) -> None:
        self.receptor_pdb, self.pose, self.cached_protein = Path(receptor_pdb), pose, cached_protein
        self.name = "Interaction analysis (ProLIF)"

    def run(self, ctx: JobContext) -> tuple[object, InteractionResult]:
        ctx.progress(10, "Building protein model")
        protein = self.cached_protein if self.cached_protein is not None else protein_to_rdkit(self.receptor_pdb)
        ctx.progress(50, "Running ProLIF fingerprint")
        result = compute_interactions(protein, self.pose)
        ctx.progress(100, "Interaction analysis complete")
        return protein, result


class SimilarityJob(Job):
    """Tanimoto similarity + Butina clustering of the prepared ligands."""

    def __init__(self, mols: list[Chem.Mol], names: list[str], cutoff: float,
                 reference: Chem.Mol | None = None, reference_name: str = "") -> None:
        self.mols, self.names, self.cutoff = mols, names, cutoff
        self.reference, self.reference_name = reference, reference_name
        self.name = f"Similarity & clustering ({len(mols)} ligands)"

    def run(self, ctx: JobContext) -> SimilarityResult:
        ctx.progress(20, "Computing fingerprints")
        if self.reference is not None:
            result = compare_to_reference(self.mols, self.names, self.reference, self.reference_name,
                                          cutoff=self.cutoff)
        else:
            result = cluster(self.mols, self.names, cutoff=self.cutoff)
        ctx.progress(100, "Clustering complete")
        return result
