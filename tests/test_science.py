"""Scientific-core tests: preparation, boxes, Vina adapter, redocking, reproducibility."""

from __future__ import annotations

import threading
import time
from pathlib import Path

import numpy as np
import pytest
Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend", exc_type=ImportError)

from mds.chem.box import DockingBox, box_from_coords, box_from_ligand_file, validate_box
from mds.chem.ligand import (entry_from_smiles, prepare_ligand, read_ligands)
from mds.chem.poses import reference_from_pdb_ligand, symmetric_rmsd
from mds.chem.receptor import ReceptorOptions, prepare_receptor
from mds.chem.report import PreparationError
from mds.chem.workflow import DockingParameters, run_docking
from mds.engines import DockingRequest, VinaEngine
from mds.engines.base import DockingCancelled
from mds.engines.vina import parse_poses

DATA = Path(__file__).parent / "data"
engine = VinaEngine()
needs_vina = pytest.mark.skipif(not engine.available(), reason="vina.exe not bundled")


@pytest.fixture(scope="module")
def trypsin(tmp_path_factory):
    out = tmp_path_factory.mktemp("trypsin")
    return prepare_receptor(DATA / "3PTB.pdb", out)


# ---------------------------------------------------------------- receptor
def test_receptor_preparation_reports_every_decision(trypsin):
    text = trypsin.report.render()
    assert "62 water molecule(s) removed" in text
    assert "Metal ion CA" in text and "retained" in text
    assert any(h.name == "BEN" and h.extracted_to and h.extracted_to.exists() for h in trypsin.hetero)
    assert "12xCYX" in text                       # six disulfide bridges detected
    assert trypsin.n_atoms > 1900
    atom_types = {l[77:79].strip() for l in trypsin.pdbqt.read_text().splitlines()
                  if l.startswith("ATOM")}
    assert {"C", "N", "OA", "HD"} <= atom_types


def test_receptor_chain_selection_error(tmp_path):
    with pytest.raises(PreparationError, match="not found"):
        prepare_receptor(DATA / "3PTB.pdb", tmp_path, ReceptorOptions(chains=["Z"]))


def test_receptor_rejects_garbage(tmp_path):
    bad = tmp_path / "bad.pdb"
    bad.write_text("this is not a structure\n")
    with pytest.raises(PreparationError):
        prepare_receptor(bad, tmp_path / "out")


# ------------------------------------------------------------------ ligand
def test_ligand_descriptors_and_pdbqt():
    p = prepare_ligand(entry_from_smiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin"), seed=3)
    assert p.properties["MW"] == pytest.approx(180.16, abs=0.01)
    assert p.properties["HBD"] == 1 and p.properties["Lipinski_pass"] == "yes"
    assert "ROOT" in p.pdbqt and "REMARK SMILES" in p.pdbqt
    assert any("Ionisable groups" in w for w in p.report.warnings)   # COOH flagged


def test_ligand_invalid_valence_is_error():
    entry = entry_from_smiles("C(C)(C)(C)(C)C", "pentavalent-carbon")
    assert entry.mol is None
    with pytest.raises(PreparationError):
        prepare_ligand(entry)


def test_ligand_undefined_stereo_warns_and_defined_is_kept():
    undefined = prepare_ligand(entry_from_smiles("CC(N)C(=O)O", "alanine?"))
    assert any("undefined stereocentre" in w for w in undefined.report.warnings)
    defined = prepare_ligand(entry_from_smiles("C[C@H](N)C(=O)O", "L-alanine"))
    assert not any("undefined stereocentre" in w for w in defined.report.warnings)
    assert "@" in defined.smiles


def test_ligand_unsupported_element_rejected():
    with pytest.raises(PreparationError, match="atom types"):
        prepare_ligand(entry_from_smiles("[Pt](Cl)(Cl)(N)N", "cisplatin"))


def test_ligand_file_formats(tmp_path):
    (tmp_path / "a.smi").write_text("CCO ethanol\nc1ccccc1 benzene\nnot_a_smiles bad\n")
    entries = read_ligands(tmp_path / "a.smi")
    assert [e.name for e in entries] == ["ethanol", "benzene", "bad"]
    assert entries[2].mol is None and entries[2].error
    (tmp_path / "b.csv").write_text("id,smiles\nm1,CCN\nm2,CCC\n")
    assert [e.name for e in read_ligands(tmp_path / "b.csv")] == ["m1", "m2"]
    with pytest.raises(PreparationError):
        read_ligands(tmp_path / "x.xyz")


def test_multiple_conformers_each_get_pdbqt():
    p = prepare_ligand(entry_from_smiles("CCCCO", "butanol"), n_conformers=3, seed=5)
    assert len(p.conformer_pdbqts) == 3


# --------------------------------------------------------------------- box
def test_box_from_coords_and_validation(trypsin):
    box = box_from_coords(np.array([[0, 0, 0], [10, 0, 0]]), padding=5)
    assert box.center == (5.0, 0.0, 0.0) and box.size[0] == 20.0
    huge = DockingBox((0, 0, 0), (60, 60, 60))
    assert any("Large search space" in w for w in validate_box(huge).warnings)
    far = DockingBox((500, 500, 500), (20, 20, 20))
    assert any("no receptor atoms" in e for e in validate_box(far, trypsin.pdbqt).errors)
    assert validate_box(DockingBox((0, 0, 0), (0, 5, 5))).errors


# -------------------------------------------------------------- vina adapter
def test_parse_poses():
    text = ("MODEL 1\nREMARK VINA RESULT:    -7.512      0.000      0.000\nATOM x\nENDMDL\n"
            "MODEL 2\nREMARK VINA RESULT:    -6.100      1.234      2.345\nATOM y\nENDMDL\n")
    poses = parse_poses(text)
    assert [p.score for p in poses] == [-7.512, -6.1]
    assert poses[1].rmsd_lb == 1.234 and poses[1].rank == 2


@needs_vina
def test_vina_version_and_validation(tmp_path):
    assert "AutoDock Vina" in engine.version()
    rep = engine.validate(DockingRequest(tmp_path / "no.pdbqt", tmp_path / "no2.pdbqt",
                                         DockingBox((0, 0, 0), (10, 10, 10)), exhaustiveness=0))
    assert len(rep.errors) >= 3


# ---------------------------------------------- redocking & reproducibility
@pytest.fixture(scope="module")
def benzamidine(trypsin):
    lig_file = next(h.extracted_to for h in trypsin.hetero if h.name == "BEN")
    reference = reference_from_pdb_ligand(lig_file, Chem.MolFromSmiles("N=C(N)c1ccccc1"))
    prepared = prepare_ligand(entry_from_smiles("N=C(N)c1ccccc1", "BEN"), seed=1, ionize="simple")
    return lig_file, reference, prepared


@needs_vina
def test_redocking_reproduces_crystal_pose(trypsin, benzamidine, tmp_path):
    lig_file, reference, prepared = benzamidine
    box = box_from_ligand_file(lig_file)
    out = run_docking(engine, trypsin.pdbqt, [prepared], box,
                      DockingParameters(exhaustiveness=8, seed=11, cpu=2), tmp_path / "run",
                      reference=reference)
    o = out[0]
    assert o.status == "ok", o.message
    assert o.reference_rmsd[0] is not None and o.reference_rmsd[0] < 2.0   # success criterion
    assert o.poses[0].score < -4.0 and o.result.engine == "AutoDock Vina"
    assert (tmp_path / "run" / "run.json").exists() and (tmp_path / "run" / "results.csv").exists()
    assert "BEN" in {p.name for p in (tmp_path / "run").iterdir()}


@needs_vina
def test_same_seed_gives_identical_results(trypsin, benzamidine, tmp_path):
    lig_file, _, prepared = benzamidine
    box = box_from_ligand_file(lig_file)
    scores = []
    for i in range(2):
        o = run_docking(engine, trypsin.pdbqt, [prepared], box,
                        DockingParameters(exhaustiveness=4, seed=99, cpu=1), tmp_path / f"r{i}")[0]
        scores.append([p.score for p in o.poses])
    assert scores[0] == scores[1]


def test_symmetric_rmsd_handles_atom_order_permutation():
    a = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1C(=O)[O-]"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(a, randomSeed=1)
    b = Chem.Mol(a)
    conf = b.GetConformer()
    # swap the two carboxylate oxygens' coordinates: chemically identical pose
    o_idx = [x.GetIdx() for x in b.GetAtoms() if x.GetSymbol() == "O"]
    p0, p1 = conf.GetAtomPosition(o_idx[0]), conf.GetAtomPosition(o_idx[1])
    conf.SetAtomPosition(o_idx[0], p1)
    conf.SetAtomPosition(o_idx[1], p0)
    assert symmetric_rmsd(a, b) == pytest.approx(0.0, abs=1e-6)


@needs_vina
def test_vina_run_can_be_cancelled(trypsin, benzamidine, tmp_path):
    lig_file, _, prepared = benzamidine
    box = box_from_ligand_file(lig_file)
    lig_pdbqt = tmp_path / "l.pdbqt"
    lig_pdbqt.write_text(prepared.pdbqt)
    flag = threading.Event()
    threading.Timer(0.5, flag.set).start()
    req = DockingRequest(trypsin.pdbqt, lig_pdbqt, box, exhaustiveness=512, cpu=1)
    t0 = time.monotonic()
    with pytest.raises(DockingCancelled):
        engine.run(req, tmp_path / "w", cancelled=flag.is_set)
    assert time.monotonic() - t0 < 10


# ---------------------------------------------------------- interaction fingerprints (ProLIF)
@pytest.fixture(scope="module")
def hiv_protease(tmp_path_factory):
    out = tmp_path_factory.mktemp("hiv_protease")
    return prepare_receptor(DATA / "1HSG.pdb", out)


@needs_vina
def test_prolif_interactions_classify_real_docked_pose(hiv_protease, tmp_path):
    from mds.chem.interactions import compute_interactions, protein_to_rdkit

    mk1 = next(h for h in hiv_protease.hetero if h.name == "MK1")
    prepared = prepare_ligand(entry_from_smiles(
        "CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H]1c2ccccc2C[C@@H]1O",
        "indinavir"), seed=1)
    box = box_from_ligand_file(mk1.extracted_to)
    outcome = run_docking(engine, hiv_protease.pdbqt, [prepared], box,
                          DockingParameters(exhaustiveness=8, seed=7, cpu=2), tmp_path / "run")[0]
    assert outcome.status == "ok"

    protein_rd = protein_to_rdkit(hiv_protease.prepared_pdb)
    assert protein_rd.GetNumAtoms() > 3000
    result = compute_interactions(protein_rd, outcome.mols[0])
    assert result.n_interactions > 0
    labels = {t for r in result.residues for t in r.interactions}
    assert any("H-bond" in l or "Hydrophobic" in l or "Van der Waals" in l for l in labels)
    assert result.network_html and "vis-network.min.js" in result.network_html
    assert "unpkg.com" not in result.network_html   # must not depend on a live CDN


def test_interaction_caveat_is_conservative():
    from mds.chem.interactions import CAVEAT
    assert "not" in CAVEAT.lower() and "verified" in CAVEAT.lower()


# ---------------------------------------------- real-world receptor robustness
@pytest.mark.parametrize("pdb_name,expect_atoms_above,expect_dropped", [
    ("1HVR", 1500, True),     # Meeko template mismatch on a residue
    ("1E66", 4500, True),     # Meeko RuntimeError on a residue with broken alt-conf geometry
])
def test_receptor_preparation_recovers_from_meeko_failures(tmp_path, pdb_name, expect_atoms_above,
                                                           expect_dropped):
    """Common real PDBs that fail the strict preparation path must now succeed, with the
    dropped residues reported."""
    result = prepare_receptor(DATA / f"{pdb_name}.pdb", tmp_path)
    assert result.n_atoms > expect_atoms_above
    text = result.report.render()
    if expect_dropped:
        assert "dropped by Meeko" in text or "removed:" in text
        assert "A from nearest extracted ligand" in text  # distance always reported
