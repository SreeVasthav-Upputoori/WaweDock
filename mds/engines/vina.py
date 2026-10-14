"""AutoDock Vina adapter (official ``vina.exe`` run as a subprocess).

Running the engine as a child process gives real cancellation (the process is killed),
isolates crashes from the GUI, and works identically from source and in the packaged app.
Results are labelled as AutoDock Vina results — never as any other program's.
"""

from __future__ import annotations

import logging
import re
import subprocess
import threading
import time
from pathlib import Path

from mds.chem.report import Report
from mds.core.paths import vina_executable
from mds.engines.base import (CancelFn, DockingCancelled, DockingError, DockingRequest,
                              DockingResult, Pose, ProgressFn)

log = logging.getLogger("mds.engines.vina")

_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
_CREATE_NO_WINDOW = 0x08000000


class VinaEngine:
    name = "AutoDock Vina"

    def __init__(self, executable: Path | None = None) -> None:
        self.executable = Path(executable) if executable else vina_executable()
        self._version: str | None = None

    # ------------------------------------------------------------------
    def available(self) -> bool:
        return self.executable.is_file()

    def version(self) -> str:
        if self._version is None:
            if not self.available():
                raise DockingError(f"AutoDock Vina executable not found: {self.executable}")
            out = subprocess.run([str(self.executable), "--version"], capture_output=True,
                                 text=True, timeout=30, creationflags=_CREATE_NO_WINDOW)
            self._version = (out.stdout or out.stderr).strip().splitlines()[0]
        return self._version

    def validate(self, request: DockingRequest) -> Report:
        rep = Report()
        if not self.available():
            rep.error(f"AutoDock Vina executable not found at {self.executable}")
        for label, path in (("receptor", request.receptor_pdbqt), ("ligand", request.ligand_pdbqt)):
            if not Path(path).is_file():
                rep.error(f"The {label} PDBQT file does not exist: {path}")
        if request.exhaustiveness < 1:
            rep.error("Exhaustiveness must be at least 1.")
        elif request.exhaustiveness < 8:
            rep.warn("Exhaustiveness below 8 (Vina default) reduces reproducibility of the best pose.")
        if not 1 <= request.num_modes <= 100:
            rep.error("Number of poses must be between 1 and 100.")
        if request.energy_range <= 0:
            rep.error("Energy range must be positive (kcal/mol).")
        if request.scoring not in ("vina", "vinardo"):
            rep.error(f"Unsupported scoring function '{request.scoring}'.")
        if any(s <= 0 for s in request.box.size):
            rep.error("Box dimensions must be positive.")
        return rep

    # ------------------------------------------------------------------
    def command(self, request: DockingRequest, out_path: Path) -> list[str]:
        b = request.box
        cmd = [
            str(self.executable),
            "--receptor", str(request.receptor_pdbqt), "--ligand", str(request.ligand_pdbqt),
            "--center_x", f"{b.center[0]:.3f}", "--center_y", f"{b.center[1]:.3f}",
            "--center_z", f"{b.center[2]:.3f}",
            "--size_x", f"{b.size[0]:.3f}", "--size_y", f"{b.size[1]:.3f}",
            "--size_z", f"{b.size[2]:.3f}",
            "--exhaustiveness", str(request.exhaustiveness),
            "--num_modes", str(request.num_modes),
            "--energy_range", f"{request.energy_range:g}",
            "--seed", str(request.seed),
            "--scoring", request.scoring,
            "--out", str(out_path),
        ]
        if request.cpu > 0:
            cmd += ["--cpu", str(request.cpu)]
        return cmd

    def run(self, request: DockingRequest, workdir: Path, progress: ProgressFn | None = None,
            cancelled: CancelFn | None = None) -> DockingResult:
        report = self.validate(request)
        if report.errors:
            raise DockingError("; ".join(report.errors))
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        out_path = workdir / "poses.pdbqt"
        out_path.unlink(missing_ok=True)
        cmd = self.command(request, out_path)
        log.info("Running Vina: %s", " ".join(cmd))

        start = time.monotonic()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                creationflags=_CREATE_NO_WINDOW, cwd=workdir)
        chunks: list[bytes] = []
        stars = [0]

        def _reader() -> None:
            assert proc.stdout is not None
            while True:
                b = proc.stdout.read(1)
                if not b:
                    break
                chunks.append(b)
                if b == b"*":
                    stars[0] += 1
                    if progress:
                        progress(min(99.0, stars[0] / 51 * 100), "Vina search in progress")

        reader = threading.Thread(target=_reader, daemon=True)
        reader.start()
        try:
            while proc.poll() is None:
                if cancelled and cancelled():
                    proc.kill()
                    proc.wait()
                    raise DockingCancelled()
                time.sleep(0.1)
        finally:
            if proc.poll() is None:
                proc.kill()
            reader.join(timeout=5)
        log_text = b"".join(chunks).decode("utf-8", "replace")
        (workdir / "vina.log").write_text(log_text)

        if proc.returncode != 0 or not out_path.exists():
            tail = "\n".join(log_text.strip().splitlines()[-6:])
            raise DockingError(f"AutoDock Vina failed (exit code {proc.returncode}):\n{tail}")

        poses = parse_poses(out_path.read_text())
        if not poses:
            raise DockingError("Vina finished but produced no poses (is the box too small or "
                               "does it miss the binding site?).")
        if progress:
            progress(100, "Vina finished")
        return DockingResult(self.name, self.version(), request.scoring, request.seed, poses, cmd,
                             log_text, time.monotonic() - start)


def parse_poses(pdbqt_text: str) -> list[Pose]:
    """Split a Vina multi-model PDBQT into :class:`Pose` objects (scores in kcal/mol)."""
    poses: list[Pose] = []
    block: list[str] = []
    score = lb = ub = 0.0
    for line in pdbqt_text.splitlines():
        if line.startswith("MODEL"):
            block = []
        elif line.startswith("ENDMDL"):
            poses.append(Pose(len(poses) + 1, score, lb, ub, "\n".join(block) + "\n"))
        else:
            m = _RESULT_RE.match(line)
            if m:
                score, lb, ub = (float(m.group(i)) for i in (1, 2, 3))
            block.append(line)
    return poses

