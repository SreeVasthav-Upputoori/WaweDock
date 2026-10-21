"""Register WaweDock in the current user's Windows Open with menu.

Run with the project venv: python tools/register_open_with.py
This adds supported file types without changing their default application.
"""
from pathlib import Path
import ctypes
import subprocess
import sys
import winreg

ROOT = Path(__file__).resolve().parent.parent
EXTENSIONS = (".pdb", ".ent", ".cif", ".mmcif", ".sdf", ".mol", ".mol2", ".smi", ".smiles")


def register() -> None:
    pythonw = Path(sys.executable).with_name("pythonw.exe")
    launcher = ROOT / "tools" / "wavedock_open.pyw"
    icon = ROOT / "mds" / "assets" / "binding-pocket.ico"
    for path in (pythonw, launcher, icon):
        if not path.is_file():
            raise FileNotFoundError(path)
    command = subprocess.list2cmdline([str(pythonw), str(launcher)]) + ' "%1"'

    def value(key, name, text):
        with winreg.CreateKey(winreg.HKEY_CURRENT_USER, key) as handle:
            winreg.SetValueEx(handle, name, 0, winreg.REG_SZ, text)

    classes = r"Software\Classes"
    prog = classes + r"\WaweDock.MolecularFile"
    value(prog, "", "WaweDock molecular file")
    value(prog, "FriendlyTypeName", "Molecular structure")
    value(prog + r"\DefaultIcon", "", f'"{icon}",0')
    value(prog + r"\shell\open", "FriendlyAppName", "WaweDock")
    value(prog + r"\shell\open\command", "", command)
    application = classes + r"\Applications\WaweDock.exe"
    value(application, "FriendlyAppName", "WaweDock")
    value(application + r"\DefaultIcon", "", f'"{icon}",0')
    value(application + r"\shell\open\command", "", command)
    capabilities = r"Software\WaweDock\Capabilities"
    value(capabilities, "ApplicationName", "WaweDock")
    value(capabilities, "ApplicationDescription", "Molecular structure viewing and docking")
    value(capabilities, "ApplicationIcon", f'"{icon}",0')
    for extension in EXTENSIONS:
        value(classes + "\\" + extension + r"\OpenWithProgids", "WaweDock.MolecularFile", "")
        value(application + r"\SupportedTypes", extension, "")
        value(capabilities + r"\FileAssociations", extension, "WaweDock.MolecularFile")
    value(r"Software\RegisteredApplications", "WaweDock", capabilities)
    ctypes.windll.shell32.SHChangeNotify(0x08000000, 0, None, None)
    print("WaweDock is registered for Open with; existing defaults are unchanged.")


if __name__ == "__main__":
    register()
