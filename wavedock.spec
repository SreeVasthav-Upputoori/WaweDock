# PyInstaller build recipe for WaweDock.   Build:  pyinstaller --noconfirm wavedock.spec
# Produces dist/WaweDock/WaweDock.exe (one-folder build: fast start-up, QtWebEngine-friendly).
from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules

datas = [("mds/data", "mds/data"), ("mds/assets", "mds/assets")]
binaries = []
hiddenimports = collect_submodules("meeko") + collect_submodules("pdb2pqr") + collect_submodules("propka")
for pkg in ("meeko", "pdb2pqr", "propka", "gemmi"):
    d, b, h = collect_all(pkg)
    datas += d
    binaries += b
    hiddenimports += h
datas += collect_data_files("rdkit")

a = Analysis(
    ["mds/__main__.py"],
    pathex=["."],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports + ["openpyxl", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineCore"],
    excludes=["tkinter", "matplotlib", "IPython", "pytest", "openmm", "pdbfixer", "openbabel"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="WaweDock",
    icon="mds/assets/wavedock.ico",
    console=False,
    version=None,
)
coll = COLLECT(exe, a.binaries, a.datas, name="WaweDock")
