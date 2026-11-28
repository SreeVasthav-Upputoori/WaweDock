"""Frame-wise protein-ligand interaction occupancy for PDB trajectories.

This module intentionally analyses supplied PDB conformations. It does not run a molecular
dynamics simulation or claim that occupancy-derived classifications are experimental facts.
"""

from __future__ import annotations

import csv
import json
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class FrameOccupancy:
    frame: int
    residue_interactions: dict[str, tuple[str, ...]]
    n_interactions: int


@dataclass(frozen=True)
class ResidueOccupancy:
    residue: str
    n_frames: int
    frame_count: int
    occupancy: float
    interactions: dict[str, float] = field(default_factory=dict)
    interaction_counts: dict[str, int] = field(default_factory=dict)
    classification: str = "unclassified"


@dataclass(frozen=True)
class OccupancyResult:
    ligand: str
    n_frames: int
    frames: list[FrameOccupancy]
    residues: list[ResidueOccupancy]
    caveat: str


def _classify(occupancy: float, n_interaction_types: int) -> str:
    if occupancy >= 0.75:
        return "anchor"
    if occupancy >= 0.50 and n_interaction_types >= 2:
        return "hotspot"
    if occupancy >= 0.25:
        return "accessible"
    return "transient"


def split_pdb_frames(source: Path) -> list[str]:
    """Return PDB text blocks from a multi-model PDB or a directory of PDB files."""
    source = Path(source)
    if source.is_dir():
        paths = sorted(p for p in source.iterdir()
                       if p.is_file() and p.suffix.lower() in (".pdb", ".ent"))
        if not paths:
            raise ValueError(f"No PDB frames found in {source}")
        return [p.read_text(encoding="utf-8", errors="replace") for p in paths]

    if not source.is_file():
        raise ValueError(f"Trajectory source not found: {source}")

    text = source.read_text(encoding="utf-8", errors="replace")
    if "MODEL " not in text:
        return [text]

    frames: list[str] = []
    current: list[str] | None = None
    for line in text.splitlines():
        if line.startswith("MODEL "):
            if current is not None:
                frames.append("\n".join(current) + "\nEND\n")
            current = []
        elif line.startswith("ENDMDL"):
            if current is not None:
                frames.append("\n".join(current) + "\nEND\n")
                current = None
        elif current is not None:
            current.append(line)
    if not frames:
        raise ValueError(f"No complete MODEL/ENDMDL blocks found in {source}")
    return frames


def _protein_and_ligand_pdb(frame_text: str, ligand_resname: str) -> tuple[str, str]:
    ligand_resname = ligand_resname.upper()
    protein_lines: list[str] = []
    ligand_lines: list[str] = []

    for line in frame_text.splitlines():
        record = line[:6].strip()
        if record not in ("ATOM", "HETATM"):
            continue
        if record == "ATOM":
            protein_lines.append(line)
            continue
        resname = line[17:20].strip().upper()
        if resname == ligand_resname:
            ligand_lines.append(line)

    if not protein_lines:
        raise ValueError("Frame contains no ATOM protein records")
    if not ligand_lines:
        raise ValueError(f"Frame contains no HETATM records for ligand residue {ligand_resname}")

    protein = "\n".join(protein_lines) + "\nEND\n"
    ligand = "\n".join(ligand_lines) + "\nEND\n"
    return protein, ligand


def classify_residue(occupancy: float, n_interaction_types: int) -> str:
    return _classify(occupancy, n_interaction_types)


def analyze_occupancy(source: Path, ligand_resname: str, ligand_smiles: str = "",
                      *, progress=None) -> OccupancyResult:
    """Analyze interaction occupancy across PDB frame(s).

    With a SMILES string, ligand coordinates are read from each frame and mapped to the supplied
    chemistry. Without SMILES, the matching ideal ligand definition is fetched from RCSB.
    """
    from rdkit import Chem

    from mds.chem.interactions import compute_interactions, protein_to_rdkit
    from mds.chem.poses import assign_bond_orders_from_template, fetch_ccd_template

    ligand_resname = ligand_resname.upper()
    frames = split_pdb_frames(source)
    frame_results: list[FrameOccupancy] = []
    residue_types: dict[str, dict[str, int]] = {}

    template = Chem.MolFromSmiles(ligand_smiles) if ligand_smiles else None
    if ligand_smiles and template is None:
        raise ValueError("Ligand SMILES could not be parsed")

    for index, frame_text in enumerate(frames):
        if progress:
            progress(index, len(frames), f"Analyzing frame {index + 1} of {len(frames)}")

        protein_pdb, ligand_pdb = _protein_and_ligand_pdb(frame_text, ligand_resname)
        with tempfile.TemporaryDirectory(prefix="wavedock-occupancy-") as folder:
            folder = Path(folder)
            receptor_path = folder / "receptor.pdb"
            ligand_path = folder / "ligand.pdb"
            receptor_path.write_text(protein_pdb, encoding="utf-8")
            ligand_path.write_text(ligand_pdb, encoding="utf-8")

            protein = protein_to_rdkit(receptor_path)
            if template is None:
                template = fetch_ccd_template(ligand_resname)

            raw = Chem.MolFromPDBFile(str(ligand_path), removeHs=True, sanitize=False,
                                      proximityBonding=True)
            if raw is None:
                raise ValueError(f"Could not read ligand coordinates in frame {index + 1}")
            ligand = assign_bond_orders_from_template(template, raw)

            result = compute_interactions(protein, ligand, diagram=False)

        residue_map: dict[str, tuple[str, ...]] = {}
        for residue in result.residues:
            labels = tuple(sorted(set(residue.interactions)))
            residue_map[residue.residue] = labels
            counts = residue_types.setdefault(residue.residue, {})
            for label in labels:
                counts[label] = counts.get(label, 0) + 1

        frame_results.append(FrameOccupancy(
            frame=index + 1,
            residue_interactions=residue_map,
            n_interactions=result.n_interactions,
        ))

    n_frames = len(frame_results)
    residues: list[ResidueOccupancy] = []
    for residue, counts in sorted(residue_types.items()):
        frame_count = sum(
            1 for frame in frame_results if residue in frame.residue_interactions
        )
        occupancy = frame_count / n_frames if n_frames else 0.0
        interaction_counts = dict(sorted(counts.items()))
        interactions = {
            label: round(count / n_frames, 6) if n_frames else 0.0
            for label, count in interaction_counts.items()
        }
        classification = classify_residue(occupancy, len(interaction_counts))
        residues.append(ResidueOccupancy(
            residue=residue,
            n_frames=n_frames,
            frame_count=frame_count,
            occupancy=round(occupancy, 6),
            interactions=interactions,
            interaction_counts=interaction_counts,
            classification=classification,
        ))

    caveat = (
        "Occupancy and residue classifications are heuristics calculated from the supplied PDB "
        "frames. They describe frame frequency, not measured affinity or experimental hotspot status."
    )
    return OccupancyResult(ligand_resname, n_frames, frame_results, residues, caveat)


def write_occupancy_csv(result: OccupancyResult, path: Path) -> Path:
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Residue", "Frame count", "Occupancy", "Classification",
                         "Interaction types", "Interaction occupancies"])
        for residue in result.residues:
            writer.writerow([
                residue.residue,
                residue.frame_count,
                f"{residue.occupancy:.4f}",
                residue.classification,
                "; ".join(sorted(residue.interaction_counts)),
                "; ".join(f"{label}={value:.4f}" for label, value in residue.interactions.items()),
            ])
    return path


def write_occupancy_json(result: OccupancyResult, path: Path) -> Path:
    path = Path(path)
    path.write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
    return path
