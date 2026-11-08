"""Read-only receptor review and preparation-option validation."""
from pathlib import Path
import gemmi
import pytest
from conftest import wait_until
from mds.chem.receptor import (ReceptorOptions, _apply_protonation_overrides,
                               _select_occupancy_conformers, _validate_protonation_overrides,
                               inspect_receptor, prepare_receptor, read_structure)
from mds.chem.report import PreparationError
from mds.chem.report import Report

DATA = Path(__file__).parent / "data"


@pytest.fixture
def incomplete_structure(tmp_path):
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for num, x in ((1,0), (4,10)):
        residue = gemmi.Residue()
        residue.name, residue.seqid, residue.het_flag = "ALA", gemmi.SeqId(num," "), "A"
        for name, offset, element in (("N",0,"N"),("CA",1,"C"),("C",2,"C"),("O",3,"O")):
            if num == 4 and name == "O": continue
            atom=gemmi.Atom()
            atom.name, atom.element, atom.pos = name, gemmi.Element(element), gemmi.Position(x+offset,0,0)
            if name == "CA" and num == 1: atom.altloc = "A"
            residue.add_atom(atom)
            if name == "CA" and num == 1:
                atom.altloc, atom.occ = "B", 0.8
                residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    other=gemmi.Chain("B")
    residue=gemmi.Residue()
    residue.name,residue.seqid,residue.het_flag="GLY",gemmi.SeqId(1," "),"A"
    for name, element in (("N","N"),("CA","C"),("C","C"),("O","O")):
        atom=gemmi.Atom(); atom.name=name; atom.element=gemmi.Element(element)
        residue.add_atom(atom)
    other.add_residue(residue); model.add_chain(other); structure.add_model(model)
    path=tmp_path/'incomplete.pdb'
    structure.write_pdb(str(path))
    return path


def test_checks_report_backbone_alternates_and_gaps_without_modifying_input(incomplete_structure):
    before=incomplete_structure.read_bytes()
    report=inspect_receptor(incomplete_structure, ReceptorOptions(chains=["A"]))
    assert not report.errors
    assert any("Incomplete backbone" in warning and "O" in warning for warning in report.warnings)
    assert any("Alternate conformations" in warning and "A, B" in warning for warning in report.warnings)
    assert any("Chain break" in warning for warning in report.warnings)
    assert incomplete_structure.read_bytes() == before
    assert not list(incomplete_structure.parent.glob('work'))
    assert not inspect_receptor(incomplete_structure, ReceptorOptions(chains=["B"])).warnings


def test_preparation_selects_highest_occupancy_conformer(incomplete_structure):
    structure = read_structure(incomplete_structure)
    assert _select_occupancy_conformers(structure) == 1
    ca_atoms = [atom for atom in structure[0][0][0] if atom.name == "CA"]
    assert len(ca_atoms) == 1
    assert ca_atoms[0].altloc in ("\0", " ")
    assert ca_atoms[0].occ == pytest.approx(1.0)


def test_protonation_override_validation_requires_matching_residue_and_disulfide():
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for number, name in ((10, "HIS"), (11, "ASP"), (12, "CYS")):
        residue = gemmi.Residue()
        residue.name, residue.seqid, residue.het_flag = name, gemmi.SeqId(number, " "), "A"
        atoms = [("CA", "C")] + ([("SG", "S")] if name == "CYS" else [])
        for atom_name, element in atoms:
            atom = gemmi.Atom()
            atom.name, atom.element = atom_name, gemmi.Element(element)
            atom.pos = gemmi.Position(float(number), 0, 0)
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    report = Report()
    valid = _validate_protonation_overrides(
        structure,
        ReceptorOptions(protonation_overrides={
            ("A", 10, ""): "HIE",
            ("A", 11, ""): "CYX",
            ("A", 99, ""): "ASH",
        }),
        report,
    )
    assert valid == {("A", 10, ""): "HIE"}
    assert len(report.errors) == 2


def test_protonation_override_rewrites_only_selected_residue(tmp_path):
    source = tmp_path / "input.pdb"
    target = tmp_path / "output.pdb"
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for number, name in ((10, "HIS"), (11, "ASP")):
        residue = gemmi.Residue()
        residue.name, residue.seqid, residue.het_flag = name, gemmi.SeqId(number, " "), "A"
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = "CA", gemmi.Element("C"), gemmi.Position(number, 0, 0)
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.write_minimal_pdb(str(source))
    _apply_protonation_overrides(source, target, {("A", 10, ""): "HIE"})
    changed = read_structure(target)
    assert changed[0]["A"][0].name == "HIE"
    assert changed[0]["A"][1].name == "ASP"


def test_checks_summarize_water_metal_cofactor_and_protonation():
    report=inspect_receptor(DATA/'3PTB.pdb',ReceptorOptions(keep_metals=False,ph=6.8))
    text=report.render()
    assert 'Waters: 62' in text
    assert 'Metal CA' in text and 'removed by the selected' in text
    assert 'BEN' in text and 'excluded from the receptor' in text
    assert 'pH 6.8' in text and 'does not predict pKa' in text
    assert inspect_receptor(DATA/'3PTB.pdb',ReceptorOptions(chains=['Z'])).errors


@pytest.mark.parametrize('options',[ReceptorOptions(ph=float('nan')),ReceptorOptions(ph=15),
                                    ReceptorOptions(drop_waters=False),
                                    ReceptorOptions(water_cutoff=1.0)])
def test_unsupported_options_fail_before_preparation(tmp_path,options):
    assert inspect_receptor(DATA/'3PTB.pdb',options).errors
    with pytest.raises(PreparationError,match='Invalid receptor'):
        prepare_receptor(DATA/'3PTB.pdb',tmp_path/'out',options)
    assert not (tmp_path/'out').exists()


def test_check_runs_in_background_and_updates_panel(qapp,tmp_path):
    from mds.gui.controller import Controller
    from mds.gui.panels.docking_panel import DockingPanel
    from mds.jobs import JobRunner
    from mds.project import Project
    runner=JobRunner()
    controller=Controller(runner)
    controller.set_project(Project.create(tmp_path/'project'))
    controller.load_structure(DATA/'3PTB.pdb')
    panel=DockingPanel(controller)
    reports=[]
    controller.report.connect(lambda title,text:reports.append((title,text)))
    try:
        panel.check_rec_btn.click()
        assert not panel.check_rec_btn.isEnabled()
        assert wait_until(lambda: panel.check_rec_btn.isEnabled() and bool(reports),10)
        assert reports[-1][0]=='Receptor check'
        assert 'warning(s)' in panel.rec_check.text()
        assert controller.receptor is None
        assert not list(controller.project.root.rglob('receptor.pdbqt'))
        panel.ph.setValue(6.8)
        assert not panel.rec_check.text()
    finally:
        runner.wait_for_done(10000)
        panel.close()
