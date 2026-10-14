"""Engine-independent docking data types and the adapter interface."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from mds.chem.box import DockingBox
from mds.chem.report import Report


class DockingError(Exception):
    """The docking engine failed or was configured incorrectly."""


class DockingCancelled(Exception):
    """The user cancelled the calculation."""


@dataclass(frozen=True)
class DockingRequest:
    receptor_pdbqt: Path
    ligand_pdbqt: Path
    box: DockingBox
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0            # kcal/mol
    cpu: int = 0                         # 0 = all cores
    seed: int = 42
    scoring: str = "vina"                # "vina" | "vinardo"


@dataclass
class Pose:
    rank: int
    score: float                         # predicted binding free energy, kcal/mol (engine score!)
    rmsd_lb: float                       # engine-reported RMSD vs. best pose (lower bound)
    rmsd_ub: float
    pdbqt: str                           # coordinates of this pose
    conformer: int = 0                   # which input conformer produced it


@dataclass
class DockingResult:
    engine: str
    engine_version: str
    scoring: str
    seed: int
    poses: list[Pose]
    command: list[str]
    log: str
    elapsed_s: float
    score_units: str = "kcal/mol"
    score_note: str = ("Predicted by the engine's empirical scoring function. This is NOT an "
                       "experimentally measured binding affinity.")


ProgressFn = Callable[[float, str], None]
CancelFn = Callable[[], bool]


class DockingEngine(Protocol):
    name: str

    def version(self) -> str: ...

    def validate(self, request: DockingRequest) -> Report: ...

    def run(self, request: DockingRequest, workdir: Path, progress: ProgressFn | None = None,
            cancelled: CancelFn | None = None) -> DockingResult: ...

