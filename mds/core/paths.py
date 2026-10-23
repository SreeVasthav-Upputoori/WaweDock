"""Locations of bundled resources, valid both from source and in a PyInstaller build."""

from __future__ import annotations

import sys
from pathlib import Path


def package_root() -> Path:
    """Directory that contains the ``mds`` package data (``mds/``)."""
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent)) / "mds"
    return Path(__file__).resolve().parent.parent


def data_dir() -> Path:
    return package_root() / "data"


def monomer_library_dir() -> Path:
    return data_dir() / "monomers"


def vina_executable() -> Path:
    return data_dir() / "vina" / "vina.exe"


def assets_dir() -> Path:
    return package_root() / "assets"


def build_edition() -> str:
    """Source installs default to Full; frozen builds declare their installed features."""
    import json
    path=package_root()/'build-info.json'
    if not path.exists(): return 'Full'
    data=json.loads(path.read_text(encoding='utf-8'))
    edition=data.get('edition')
    if edition not in ('Full','Protein'): raise ValueError('Invalid application build edition.')
    return edition
