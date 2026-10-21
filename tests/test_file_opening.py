"""Explorer and command-line files use the same import path as drag and drop."""
from pathlib import Path
import shutil

from rdkit import Chem
from mds.gui.main_window import MainWindow
from mds.project import Project


def test_open_files_imports_structure_and_ligand_with_spaces(qapp, settings, tmp_path):
    window = MainWindow(settings, use_viewer=False)
    window.controller.set_project(Project.create(tmp_path / "project"))
    structure = tmp_path / "protein with spaces.PDB"
    shutil.copyfile(Path(__file__).parent / "data" / "3PTB.pdb", structure)
    ligand = tmp_path / "ligand with spaces.sdf"
    with Chem.SDWriter(str(ligand)) as writer:
        writer.write(Chem.MolFromSmiles("CCO"))
    try:
        window.open_files([structure, ligand])
        assert window.controller.source_pdb.name == structure.name
        assert window.controller.source_pdb.read_bytes() == structure.read_bytes()
        assert len(window.controller.entries) == 1
        assert window.controller.entries[0].mol.GetNumAtoms() == 3
    finally:
        window.close()


def test_open_files_reports_missing_and_unsupported_files(qapp, settings, tmp_path):
    window = MainWindow(settings, use_viewer=False)
    errors = []
    window.controller.error.disconnect(window._show_error)
    window.controller.error.connect(lambda title, message: errors.append((title, message)))
    unsupported = tmp_path / "document.pdf"
    unsupported.write_bytes(b"example")
    try:
        window.open_files([tmp_path / "missing.pdb", unsupported])
        assert len(errors) == 2
        assert "File not found" in errors[0][1]
        assert "Unsupported file type" in errors[1][1]
        assert window.controller.project is None
    finally:
        window.close()
