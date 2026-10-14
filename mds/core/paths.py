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
