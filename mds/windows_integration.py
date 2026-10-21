"""Per-user Explorer integration for source, pip, and frozen Windows builds."""
from __future__ import annotations

import ctypes
from pathlib import Path
import subprocess
import sys

from mds.core.paths import assets_dir

EXTENSIONS = (".pdb", ".ent", ".cif", ".mmcif", ".sdf", ".mol", ".mol2", ".smi", ".smiles")
PROGID = "WaweDock.MolecularFile"
CLASSES = r"Software\Classes"
CAPABILITIES = r"Software\WaweDock\Capabilities"


def launch_command() -> str:
    if getattr(sys, "frozen", False):
        parts = [sys.executable]
    else:
        pythonw = Path(sys.executable).with_name("pythonw.exe")
        launcher = Path(__file__).with_name("windows_launcher.py")
        parts = [str(pythonw), str(launcher)]
    return subprocess.list2cmdline(parts) + ' "%1"'


def registration_entries() -> list[tuple[str, str, str]]:
    command = launch_command()
    icon = f'"{assets_dir() / "binding-pocket.ico"}",0'
    prog = CLASSES + "\\" + PROGID
    application = CLASSES + r"\Applications\WaweDock.exe"
    entries = [
        (prog, "", "WaweDock molecular file"),
        (prog, "FriendlyTypeName", "Molecular structure"),
        (prog + r"\DefaultIcon", "", icon),
        (prog + r"\Application", "ApplicationName", "WaweDock"),
        (prog + r"\Application", "ApplicationIcon", icon),
        (prog + r"\Application", "AppUserModelID", "WaweDock.Desktop"),
        (prog + r"\shell\open", "FriendlyAppName", "WaweDock"),
        (prog + r"\shell\open\command", "", command),
        (application, "FriendlyAppName", "WaweDock"),
        (application + r"\DefaultIcon", "", icon),
        (application + r"\shell\open\command", "", command),
        (CAPABILITIES, "ApplicationName", "WaweDock"),
        (CAPABILITIES, "ApplicationDescription", "Molecular structure viewing and docking"),
        (CAPABILITIES, "ApplicationIcon", icon),
        (r"Software\RegisteredApplications", "WaweDock", CAPABILITIES),
    ]
    for extension in EXTENSIONS:
        entries.extend([
            (CLASSES + "\\" + extension + r"\OpenWithProgids", PROGID, ""),
            (application + r"\SupportedTypes", extension, ""),
            (CAPABILITIES + r"\FileAssociations", extension, PROGID),
        ])
    return entries


def _read(registry, key: str, name: str):
    try:
        with registry.OpenKey(registry.HKEY_CURRENT_USER, key) as handle:
            return registry.QueryValueEx(handle, name)[0]
    except FileNotFoundError:
        return None


def register() -> bool:
    """Repair changed paths, but avoid registry writes and shell refreshes otherwise."""
    if sys.platform != "win32":
        return False
    import winreg
    changed = False
    for key, name, value in registration_entries():
        if _read(winreg, key, name) == value:
            continue
        with winreg.CreateKey(winreg.HKEY_CURRENT_USER, key) as handle:
            winreg.SetValueEx(handle, name, 0, winreg.REG_SZ, value)
        changed = True
    if changed:
        ctypes.windll.shell32.SHChangeNotify(0x08000000, 0, None, None)
    return changed


def unregister() -> bool:
    """Remove only this installation's values; preserve newer installations/defaults."""
    if sys.platform != "win32":
        return False
    import winreg
    command_key = CLASSES + "\\" + PROGID + r"\shell\open\command"
    if _read(winreg, command_key, "") != launch_command():
        return False
    entries = registration_entries()
    for key, name, value in reversed(entries):
        if _read(winreg, key, name) == value:
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, key, 0, winreg.KEY_SET_VALUE) as handle:
                winreg.DeleteValue(handle, name)
    # Delete empty keys only. Never remove a key containing another app's values.
    keys = {key for key, _, _ in entries if key != r"Software\RegisteredApplications"}
    owned = (CLASSES + "\\" + PROGID, CLASSES + r"\Applications\WaweDock.exe", r"Software\WaweDock")
    for original in list(keys):
        key = original.rsplit("\\", 1)[0]
        while any(key == root or key.startswith(root + "\\") for root in owned):
            keys.add(key)
            key = key.rsplit("\\", 1)[0]
    for key in sorted(keys, key=lambda k: k.count("\\"), reverse=True):
        try:
            winreg.DeleteKey(winreg.HKEY_CURRENT_USER, key)
        except OSError:
            pass
    ctypes.windll.shell32.SHChangeNotify(0x08000000, 0, None, None)
    return True
