# PyInstaller build recipe for WaweDock.   Build:  pyinstaller --noconfirm wavedock.spec
# Produces dist/WaweDock/WaweDock.exe (one-folder build: fast start-up, QtWebEngine-friendly).
import json
import os
from importlib.util import find_spec
from pathlib import Path
from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules, copy_metadata

edition=os.environ.get('WAVEDOCK_BUILD_EDITION','Full')
if edition not in ('Full','Protein'): raise ValueError('WAVEDOCK_BUILD_EDITION must be Full or Protein')
metadata=Path('build/release-metadata/build-info.json')
metadata.parent.mkdir(parents=True,exist_ok=True)
metadata.write_text(json.dumps({'edition':edition}),encoding='utf-8')
datas = [("mds/data", "mds/data"), ("mds/assets", "mds/assets"),(str(metadata),"mds")]
binaries = []
hiddenimports=[]
packages=('meeko','pdb2pqr','propka','gemmi') if edition=='Full' else ('gemmi',)
for pkg in packages:
    d,b,h=collect_all(pkg)
    datas+=d; binaries+=b;hiddenimports+=h
if edition=='Full': datas+=collect_data_files('rdkit')
# Optional repair backends: templates, CPU plugin DLLs and metadata used by previews.
for package,distribution in (("pdbfixer","pdbfixer"),("openmm","openmm"),("Bio","biopython")):
    if find_spec(package) is not None:
        d,b,h=collect_all(package)
        datas+=d+copy_metadata(distribution)
        binaries+=b
        hiddenimports+=h
        if package=="openmm":
            # OpenMM wheels place native libraries beside the Python package.
            import openmm.version
            native_root=Path(openmm.version.openmm_library_path).parent
            if native_root.exists(): datas.append((str(native_root),"OpenMM.libs"))

a = Analysis(
    ["mds/__main__.py"],
    pathex=["."],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports + ["openpyxl", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineCore"],
    # openbabel is GPL-2.0 and kept out of the distribution deliberately (see pyproject.toml's
    # "interactions" extra and mds/chem/interactions.py) so this build is not GPL-bound; the
    # "Interactions (ProLIF)" feature is unavailable in this .exe unless it is installed
    # separately into a from-source environment.
    excludes=["tkinter", "matplotlib", "IPython", "pytest", "openbabel"]+
        (["rdkit","meeko","pdb2pqr","propka","prolif","MDAnalysis"] if edition=="Protein" else []),
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="WaweDock",
    icon="mds/assets/binding-pocket.ico",
    console=False,
    version=None,
)
coll = COLLECT(exe, a.binaries, a.datas, name="WaweDock")
