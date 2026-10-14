"""Project folders: inputs, prepared files, runs and a JSON state file.

::

    <project>/
      project.json     state: receptor, ligands, box, parameters, runs
      originals/       imported files, never overwritten (new names get a numeric suffix)
      prepared/        receptor PDBQT, prepared ligands, preparation reports
      runs/            one folder per docking run (see mds.chem.workflow)
      exports/         CSV / SDF / Excel / images
"""

from __future__ import annotations

import json
import logging
import shutil
from datetime import datetime
from pathlib import Path

from mds import __version__

log = logging.getLogger("mds.project")
SUBFOLDERS = ("originals", "prepared", "runs", "exports")
STATE_FILE = "project.json"


class Project:
    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.state: dict = {}

    # -------------------------------------------------------------- lifecycle
    @classmethod
    def create(cls, root: Path) -> Project:
        root = Path(root)
        if (root / STATE_FILE).exists():
            raise FileExistsError(f"A project already exists in {root}")
        for sub in SUBFOLDERS:
            (root / sub).mkdir(parents=True, exist_ok=True)
        project = cls(root)
        project.state = {"application": "WaweDock", "version": __version__,
                         "created": datetime.now().isoformat(timespec="seconds"), "runs": []}
        project.save()
        log.info("Project created: %s", root)
        return project

    @classmethod
    def open(cls, root: Path) -> Project:
        root = Path(root)
        state_path = root / STATE_FILE
        if not state_path.is_file():
            raise FileNotFoundError(f"No {STATE_FILE} in {root}; not a WaweDock project folder")
        project = cls(root)
        project.state = json.loads(state_path.read_text(encoding="utf-8"))
        for sub in SUBFOLDERS:
            (root / sub).mkdir(exist_ok=True)
        log.info("Project opened: %s", root)
        return project

    def save(self) -> None:
        self.state["saved"] = datetime.now().isoformat(timespec="seconds")
        tmp = self.root / (STATE_FILE + ".tmp")
        tmp.write_text(json.dumps(self.state, indent=2), encoding="utf-8")
        tmp.replace(self.root / STATE_FILE)  # atomic: a crash never leaves half a file

    @property
    def name(self) -> str:
        return self.root.name

    # ------------------------------------------------------------------ files
    def import_original(self, source: Path) -> Path:
        """Copy ``source`` into ``originals/`` without overwriting anything."""
        source = Path(source)
        dest = self.root / "originals" / source.name
        n = 1
        while dest.exists():
            if dest.read_bytes() == source.read_bytes():
                return dest                      # identical file already imported
            dest = dest.with_name(f"{source.stem}_{n}{source.suffix}")
            n += 1
        shutil.copy2(source, dest)
        return dest

    def prepared_dir(self, name: str) -> Path:
        path = self.root / "prepared" / name
        path.mkdir(parents=True, exist_ok=True)
        return path

    def new_run_dir(self) -> Path:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        path = self.root / "runs" / f"run-{stamp}"
        n = 1
        while path.exists():
            path = self.root / "runs" / f"run-{stamp}-{n}"
            n += 1
        path.mkdir(parents=True)
        return path

    def rel(self, path: Path) -> str:
        try:
            return str(Path(path).resolve().relative_to(self.root.resolve())).replace("\\", "/")
        except ValueError:
            return str(path)

    def abs(self, rel: str) -> Path:
        p = Path(rel)
        return p if p.is_absolute() else self.root / p

    # ---------------------------------------------------------------- archive
    def archive(self, destination: Path) -> Path:
        """Zip the whole project (for archival / sharing)."""
        base = Path(destination).with_suffix("")
        return Path(shutil.make_archive(str(base), "zip", self.root))

    def duplicate(self, new_root: Path) -> Project:
        new_root = Path(new_root)
        if new_root.exists() and any(new_root.iterdir()):
            raise FileExistsError(f"{new_root} is not empty")
        shutil.copytree(self.root, new_root, dirs_exist_ok=True)
        return Project.open(new_root)
