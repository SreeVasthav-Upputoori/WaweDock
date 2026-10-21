"""Exercise the bundled 3Dmol viewer and its real Python edit bridge."""
from pathlib import Path
import json

import gemmi
import pytest
from conftest import wait_until
from PySide6.QtCore import QPoint, Qt
from PySide6.QtTest import QTest

from mds.chem.edit import atom_records, cleanup_selection, edit_atoms
from mds.gui.main_window import MainWindow
from mds.project import Project

DATA = Path(__file__).parent / "data"


def js(viewer, expression):
    results = []
    viewer.page().runJavaScript(expression, results.append)
    assert wait_until(lambda: bool(results), 10), "JavaScript did not return"
    return results[0]


@pytest.fixture
def live_window(qapp, settings, tmp_path):
    window = MainWindow(settings, use_viewer=True)
    assert window.viewer is not None
    window.controller.set_project(Project.create(tmp_path / "project"))
    window.show()
    assert wait_until(lambda: window.viewer._ready, 20), "Viewer failed to initialise"
    yield window
    window.close()


def test_label_toggle_clears_nearby_residue_labels(live_window):
    w = live_window
    w.controller.load_structure(DATA / "3PTB.pdb")
    assert wait_until(lambda: bool(js(w.viewer, "objects['3PTB'] && objects['3PTB'].editable")), 10)
    ligand = '\n'.join(line for line in w.controller.source_pdb.read_text().splitlines()
                       if line.startswith('HETATM') and line[17:20].strip() == 'BEN')
    w.viewer.add_object('probe', ligand, 'pdb', 'ligand')
    assert js(w.viewer, "labelResidues(true); showAround('probe',4.5,'3PTB'); viewer.labels.length > 0")
    assert js(w.viewer, "setObjectLabels('3PTB',false); viewer.labels.length") == 0
    assert js(w.viewer, "setRepresentation('3PTB','sticks'); showAround('probe',4.5,'3PTB'); viewer.labels.length") == 0
    assert js(w.viewer, "labelResidues(true); viewer.labels.length > 0")
    assert js(w.viewer, "labelResidues(false); viewer.labels.length") == 0
    assert js(w.viewer, "measureShapes.push(viewer.addLabel('distance', {position:{x:0,y:0,z:0}})); setObjectLabels('3PTB',false); viewer.labels.length") == 1
    assert js(w.viewer, "clearMeasurements(); viewer.labels.length") == 0


def test_residue_single_click_component_double_click_and_overlay_highlights(live_window, tmp_path):
    w = live_window
    records = [
        ('ATOM', 'N', 'ALA', 'A', 1, -1.2, 0, 'N'),
        ('ATOM', 'CA', 'ALA', 'A', 1, 0, 0, 'C'),
        ('ATOM', 'C', 'ALA', 'A', 1, 1.2, 0, 'C'),
        ('ATOM', 'CA', 'GLY', 'B', 2, 8, 0, 'C'),
        ('HETATM', 'O', 'HOH', 'W', 1, 0, 8, 'O'),
        ('HETATM', 'H1', 'HOH', 'W', 1, .8, 8, 'H'),
        ('HETATM', 'O', 'HOH', 'W', 2, 8, 8, 'O'),
        ('HETATM', 'ZN', 'ZN', 'Z', 3, 16, 8, 'ZN'),
        ('HETATM', 'N', 'MSE', 'C', 3, 14.8, 0, 'N'),
        ('HETATM', 'CA', 'MSE', 'C', 3, 16, 0, 'C'),
        ('HETATM', 'C', 'MSE', 'C', 3, 17.2, 0, 'C'),
    ]
    path = tmp_path / 'clicks.pdb'
    path.write_text(''.join(
        f'{record:<6}{i:5d} {atom:^4} {resn:>3} {chain}{resi:4d}    '
        f'{x:8.3f}{y:8.3f}{0:8.3f}  1.00 20.00          {element:>2}  \n'
        for i, (record, atom, resn, chain, resi, x, y, element) in enumerate(records, 1)) + 'END\n')
    w.controller.load_structure(path)
    assert wait_until(lambda: bool(js(w.viewer, "objects.clicks && objects.clicks.editable")), 10)
    js(w.viewer, "setRepresentation('clicks','spheres'); setSelectionMode('chain'); true")
    # Click selects exactly one residue even if a broader tool mode was chosen.
    point = json.loads(js(w.viewer, "JSON.stringify(viewer.modelToScreen(objects.clicks.model.selectedAtoms({atom:'CA',chain:'A'})[0]))"))
    QTest.mouseClick(w.viewer.focusProxy() or w.viewer, Qt.MouseButton.LeftButton,
                     pos=QPoint(round(point['x']), round(point['y'])))
    assert wait_until(lambda: len(w.controller.selection.get('clicks', [])) == 3)
    # Dispatch the real DOM double-click listener against rendered geometry.
    assert js(w.viewer, """(() => {
      const p = viewer.modelToScreen(objects.clicks.model.selectedAtoms({atom:'CA',chain:'A'})[0]);
      container.dispatchEvent(new MouseEvent('dblclick', {button:0,clientX:p.x,clientY:p.y,bubbles:true,detail:2}));
      return selection.clicks.size;
    })()""") == 7
    js(w.viewer, "pickMolecularSelection('clicks',objects.clicks.model.selectedAtoms({resn:'HOH'})[0],{button:0,detail:1}); true")
    assert wait_until(lambda: len(w.controller.selection.get('clicks', [])) == 2)
    js(w.viewer, "pickMolecularSelection('clicks',objects.clicks.model.selectedAtoms({resn:'HOH'})[0],{button:0},true); true")
    assert wait_until(lambda: len(w.controller.selection.get('clicks', [])) == 3)
    assert js(w.viewer, "selectionShapes.length") == 1
    assert js(w.viewer, "selectionShapes[0].stylespec.color") == '#ffff00'
    assert js(w.viewer, "selectionShapes[0].boundingSphere.radius > 1.5")
    # Prepared overlays highlight visually without gaining editing identities.
    w.viewer.add_object('prepared-overlay', path.read_text(), 'pdb', 'protein')
    js(w.viewer, "pickMolecularSelection('prepared-overlay',objects['prepared-overlay'].model.selectedAtoms({resn:'HOH'})[0],{button:0}); true")
    assert js(w.viewer, "visualSelectionCount()") == 2
    assert not w.controller.selection
    assert js(w.viewer, "objects['prepared-overlay'].editable") is False
    assert js(w.viewer, "clearSelection(); visualSelectionCount()") == 0
    assert js(w.viewer, "selectionShapes.length") == 0


def test_protein_styles_do_not_keep_all_atom_overlay(live_window):
    w = live_window
    w.controller.load_structure(DATA / "3PTB.pdb")
    name = w.controller.source_pdb.stem
    key = json.dumps(name)
    assert wait_until(lambda: bool(js(w.viewer, f"objects[{key}] && objects[{key}].editable")), 10)
    for rep in ("ribbon", "backbone", "cartoon", "ribbon"):
        w.viewer.set_representation(name, rep)
        styles = json.loads(js(w.viewer, f"""JSON.stringify((() => {{
            const m = objects[{key}].model;
            return {{rep: objects[{key}].rep,
                ca: m.selectedAtoms({{hetflag: false, atom: 'CA'}})[0].style,
                side: m.selectedAtoms({{hetflag: false, atom: 'CB'}})[0].style,
                water: m.selectedAtoms({{resn: 'HOH'}})[0].style,
                ligand: m.selectedAtoms({{resn: 'BEN'}})[0].style}};
        }})())"""))
        assert styles["rep"] == rep
        assert "line" not in styles["side"]
        assert "sphere" not in styles["side"]
        if rep == "backbone":
            assert "stick" in styles["ca"]
            assert "stick" not in styles["side"]
        else:
            assert "cartoon" in styles["ca"]
            assert styles["ca"]["cartoon"]["style"] == ("trace" if rep == "ribbon" else "rectangle")
        assert "sphere" in styles["water"]
        assert "stick" in styles["ligand"]


def test_context_hit_testing_and_molecular_controls(live_window):
    w = live_window
    w.controller.load_structure(DATA / "3PTB.pdb")
    assert wait_until(lambda: bool(js(w.viewer, "objects['3PTB'] && objects['3PTB'].editable")), 10)
    def run(code):
        return js(w.viewer, code)

    # Use real rendered spheres to exercise the bundled raycaster.
    run("setRepresentation('3PTB', 'spheres'); zoomTo('3PTB'); true")
    payload = json.loads(run("""JSON.stringify((() => {
      const a = objects['3PTB'].model.selectedAtoms({hetflag:false, atom:'CA'})[0];
      const p = viewer.modelToScreen(a);
      return contextPayload({pageX:p.x,pageY:p.y,clientX:p.x,clientY:p.y});
    })())"""))
    assert payload["name"] == "3PTB"
    assert payload["editable"]
    assert payload["atomIndex"] is not None
    assert payload["selectionCount"] == 0  # Right-click never changes selection.
    assert run("contextPayload({pageX:-10000,pageY:-10000,clientX:0,clientY:0}).name === null")
    run(f"selectContextTarget('3PTB', {payload['atomIndex']}, 'atom'); true")
    assert wait_until(lambda: len(w.controller.selection.get('3PTB', [])) == 1)

    count = run("objects['3PTB'].model.selectedAtoms({}).length")
    run("setComponentVisible('3PTB','waters',false); setCustomColor('3PTB','#123456'); true")
    assert run("objects['3PTB'].model.selectedAtoms({resn:'HOH'}).every(a => !Object.keys(a.style).length)")
    assert run("objects['3PTB'].model.selectedAtoms({hetflag:false,atom:'CA'})[0].style.sphere.color") == '#123456'
    run("setRepresentation('3PTB','sticks'); true")
    assert run("objects['3PTB'].model.selectedAtoms({resn:'HOH'}).every(a => !Object.keys(a.style).length)")
    run("setComponentVisible('3PTB','waters',true); setComponentVisible('3PTB','ligands',false); true")
    assert run("objects['3PTB'].model.selectedAtoms({resn:'BEN'}).every(a => !Object.keys(a.style).length)")
    assert run("objects['3PTB'].model.selectedAtoms({resn:'HOH'}).every(a => !!a.style.sphere)")
    assert run("objects['3PTB'].model.selectedAtoms({}).length") == count
    run("setColorScheme('3PTB','residue'); setObjectLabels('3PTB',true); true")
    assert run("objects['3PTB'].customColor === null")
    assert run("objects['3PTB'].objectLabels.length > 0")
    run("labelResidues(false); true")
    assert run("objects['3PTB'].objectLabels.length") == 0
    run("labelResidues(true); true")
    assert run("objects['3PTB'].labels")
    w.viewer.set_visible('3PTB', False)
    assert not w.project_explorer.is_visible('3PTB')
    assert run("objects['3PTB'].objectLabels.length") == 0
    w.viewer.set_visible('3PTB', True)
    assert w.project_explorer.is_visible('3PTB')
    run("setObjectLabels('3PTB',false); setBackgroundColor('#123456'); true")
    assert run("backgroundColor") == '#123456'


def test_surface_switching_uses_numeric_surface_ids(live_window):
    w = live_window
    w.controller.load_structure(DATA / "3PTB.pdb")
    assert wait_until(lambda: bool(js(w.viewer, "objects['3PTB'] && objects['3PTB'].editable")), 10)
    assert js(w.viewer, "setRepresentation('3PTB','surface'); typeof objects['3PTB'].surf") == 'number'
    assert js(w.viewer, "setColorScheme('3PTB','residue'); typeof objects['3PTB'].surf") == 'number'
    assert js(w.viewer, "setRepresentation('3PTB','cartoon'); objects['3PTB'].surf === undefined")
    assert js(w.viewer, "Object.keys(viewer.surfaces).length") == 0


def test_context_event_opens_native_menu_and_hydrogens_stay_hidden(live_window, monkeypatch):
    from PySide6.QtWidgets import QMenu
    w = live_window
    w.controller.load_structure(DATA / "3PTB.pdb")
    assert wait_until(lambda: bool(js(w.viewer, "objects['3PTB'] && objects['3PTB'].editable")), 10)
    opened = []
    monkeypatch.setattr(QMenu, "popup", lambda menu, point: opened.append(
        [action.text() for action in menu.actions()]))
    js(w.viewer, "container.dispatchEvent(new MouseEvent('contextmenu', {clientX:0,clientY:0,bubbles:true})); true")
    assert wait_until(lambda: bool(opened))
    assert "Background Settings" in opened[0]
    # A separate small molecule ensures the H toggle exercises actual hydrogens.
    pdb = (
        "HETATM    1  C1  LIG B   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "HETATM    2  H1  LIG B   1       1.000   0.000   0.000  1.00 20.00           H  \nEND\n"
    )
    w.viewer.add_object('hydrogen-test', pdb, 'pdb', 'ligand')
    assert js(w.viewer, "objects['hydrogen-test'].model.selectedAtoms({elem:'H'}).length") == 1
    assert js(w.viewer, """(() => {
      setComponentVisible('hydrogen-test','hydrogens',false);
      setRepresentation('hydrogen-test','ball-and-stick');
      const m = objects['hydrogen-test'].model;
      return m.selectedAtoms({elem:'H'}).every(a=>!Object.keys(a.style).length)
        && m.selectedAtoms({elem:'C'}).every(a=>!!a.style.stick);
    })()""")


def test_results_workspace_displays_selected_pose_and_keeps_viewer_state(live_window):
    from test_results_workspace import sample_outcomes
    w = live_window
    w.controller.load_structure(DATA / "3PTB.pdb")
    assert wait_until(lambda: bool(js(w.viewer, "objects['3PTB'] && objects['3PTB'].editable")), 10)
    js(w.viewer, "setSelectionMode('atom'); pickSelection('3PTB', objects['3PTB'].model.selectedAtoms({atom:'CA'})[0], {}); true")
    assert wait_until(lambda: len(w.controller.selection.get('3PTB', [])) == 1)
    original_viewer = w.viewer
    w.controller.outcomes = sample_outcomes()
    w.controller.docking_finished.emit(w.controller.outcomes)
    assert w.workspace.tabs.currentIndex() == 2
    assert not w.results_panel.isHidden()
    assert wait_until(lambda: bool(js(w.viewer, "!!objects['pose: ligand B #1']")), 10)
    w.results_panel.select_pose(w.controller.outcomes[0], 1)
    assert wait_until(lambda: bool(js(w.viewer, "!!objects['pose: ligand A #2']")), 10)
    assert js(w.viewer, "objects['pose: ligand B #1'] === undefined")
    camera = js(w.viewer, "JSON.stringify(viewer.getView())")
    w.workspace.select(0)
    w.workspace.select(1)
    w.workspace.select(2)
    assert w.viewer is original_viewer
    assert js(w.viewer, "JSON.stringify(viewer.getView())") == camera
    assert len(w.controller.selection.get('3PTB', [])) == 1
    w.results_panel.search.setText('no matching pose')
    assert wait_until(lambda: bool(js(w.viewer, "objects['pose: ligand A #2'] === undefined")), 10)


def test_real_viewer_atom_identity_delete_and_undo(live_window):
    w = live_window
    c = w.controller
    c.load_structure(DATA / "3PTB.pdb")
    name = c.source_pdb.stem
    key = json.dumps(name)
    assert wait_until(lambda: bool(js(w.viewer, f"objects[{key}] && objects[{key}].editable")), 10)
    original = c.source_pdb.read_bytes()
    original_count = len(atom_records(gemmi.read_structure(str(c.source_pdb))))
    # Every atom, including waters and HETATM records, survives viewer import.
    assert js(w.viewer, f"objects[{key}].model.selectedAtoms({{}}).length") == original_count
    js(w.viewer, f"setSelectionMode('atom'); pickSelection({key}, objects[{key}].model.selectedAtoms({{}})[0], {{}}); true")
    assert wait_until(lambda: sum(map(len, c.selection.values())) == 1)
    js(w.viewer, "document.dispatchEvent(new KeyboardEvent('keydown', {key:'Delete'})); true")
    assert wait_until(lambda: c.source_pdb.stem != name)
    edited_name = c.source_pdb.stem
    assert len(atom_records(gemmi.read_structure(str(c.source_pdb)))) == original_count - 1
    assert w.project_explorer.find(edited_name) is not None
    assert w.project_explorer.find(name) is None
    js(w.viewer, "document.dispatchEvent(new KeyboardEvent('keydown', {key:'z', ctrlKey:true})); true")
    assert wait_until(lambda: c.source_pdb.stem == name)
    assert c.source_pdb.read_bytes() == original
    c.redo_structure_edit()
    assert c.source_pdb.stem == edited_name


def test_mouse_pick_freehand_lasso_and_radius(live_window, tmp_path):
    w, c = live_window, live_window.controller
    # Three isolated atoms make hit-testing deterministic. The third is well
    # outside the radius and polygon around the first two atoms.
    pdb = tmp_path / "points.pdb"
    pdb.write_text(
        "HETATM    1 ZN    ZN A   1       0.000   0.000   0.000  1.00 20.00          ZN  \n"
        "HETATM    2 ZN    ZN A   2       3.000   0.000   0.000  1.00 20.00          ZN  \n"
        "HETATM    3 ZN    ZN B   3      20.000  10.000   0.000  1.00 20.00          ZN  \nEND\n"
    )
    c.load_structure(pdb)
    assert wait_until(lambda: bool(js(w.viewer, "objects.points && objects.points.editable")), 10)
    js(w.viewer, "setSelectionMode('atom'); setRepresentation('points','spheres'); zoomTo('points'); true")
    point = json.loads(js(w.viewer, "JSON.stringify(viewer.modelToScreen(objects.points.model.selectedAtoms({})[0]))"))
    target = w.viewer.focusProxy() or w.viewer
    QTest.mouseClick(target, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(round(point["x"]), round(point["y"])))
    assert wait_until(lambda: len(c.selection.get("points", [])) == 1), "3Dmol mouse picking failed"
    assert c.selection["points"][0]["atom_index"] == 0
    js(w.viewer, "selectWithinRadius(3.1); true")
    assert wait_until(lambda: len(c.selection.get("points", [])) == 2)
    assert {r["atom_index"] for r in c.selection["points"]} == {0, 1}
    # Real event listeners draw a nonrectangular polygon using page coordinates.
    assert js(w.viewer, """(() => {
      clearSelection();
      const p = viewer.modelToScreen(objects.points.model.selectedAtoms({})[0]);
      const points = [[p.x-8,p.y-10],[p.x+12,p.y-10],[p.x+9,p.y+10],[p.x-12,p.y+8]];
      container.dispatchEvent(new MouseEvent('mousedown', {clientX:points[0][0],clientY:points[0][1],altKey:true,button:0,bubbles:true}));
      points.slice(1).forEach(q => document.dispatchEvent(new MouseEvent('mousemove', {clientX:q[0],clientY:q[1],altKey:true,bubbles:true})));
      document.dispatchEvent(new MouseEvent('mouseup', {altKey:true,button:0,bubbles:true}));
      return pointInPolygon(p, points.map(q=>({x:q[0],y:q[1]})));
    })()""")
    assert wait_until(lambda: len(c.selection.get("points", [])) == 1)
    assert c.selection["points"][0]["atom_index"] == 0


def test_exact_atom_deletion_keeps_other_altlocs_and_models(tmp_path):
    st = gemmi.Structure()
    for mi in range(2):
        model = gemmi.Model(mi + 1)
        chain = gemmi.Chain("A")
        for icode in (" ", "A"):
            residue = gemmi.Residue()
            residue.name, residue.seqid = "ALA", gemmi.SeqId(25, icode)
            for alt in ("A", "B"):
                atom = gemmi.Atom()
                atom.name, atom.element, atom.altloc = "CA", gemmi.Element("C"), alt
                atom.pos = gemmi.Position(mi * 10, ord(icode), ord(alt))
                residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
        st.add_model(model)
    src, dst = tmp_path / "source.pdb", tmp_path / "edited.pdb"
    st.write_minimal_pdb(str(src))
    records = atom_records(gemmi.read_structure(str(src)))
    assert edit_atoms(src, dst, [records[3]]) == 1
    remaining = atom_records(gemmi.read_structure(str(dst)))
    identity = lambda r: (r["model"], r["chain"], r["resi"], r["icode"], r["altloc"])
    assert {identity(r) for r in remaining} == {identity(r) for i, r in enumerate(records) if i != 3}


def test_cleanup_retains_selected_water_and_ligand(tmp_path):
    src, dst = DATA / "3PTB.pdb", tmp_path / "cleaned.pdb"
    st = gemmi.read_structure(str(src))
    records = atom_records(st)
    water = next(r for r in records if r["resn"] == "HOH")
    retained = [r for r in records if r["resn"] == "BEN" or
                (r["resn"] == "HOH" and r["resi"] == water["resi"])]
    victims = cleanup_selection(src, waters=True, hetero=True, retain=retained)
    edit_atoms(src, dst, victims)
    remaining = atom_records(gemmi.read_structure(str(dst)))
    assert sum(r["resn"] == "HOH" for r in remaining) == 1
    assert sum(r["resn"] == "BEN" for r in remaining) == sum(r["resn"] == "BEN" for r in records)


def test_history_branch_and_source_model_names(qapp, settings, tmp_path):
    w = MainWindow(settings, use_viewer=False)
    c = w.controller
    c.set_project(Project.create(tmp_path / "project"))
    try:
        c.load_structure(DATA / "3PTB.pdb")
        original = c.source_pdb
        for _ in range(2):
            name = c.source_pdb.stem
            records = atom_records(gemmi.read_structure(str(c.source_pdb)))
            c.set_selection({name: [records[0]]}, "")
            c.delete_selected_residues()
        c.undo_structure_edit()
        c.undo_structure_edit()
        assert c.source_pdb == original
        c.redo_structure_edit()
        c.set_selection({c.source_pdb.stem: [{"chain":"A", "resi":57}]}, "")
        c.delete_selected_residues()
        assert not c._redo_structures
        assert w.project_explorer.find(c.source_pdb.stem) is not None
    finally:
        w.close()


@pytest.mark.parametrize("suffix", [".pdb", ".cif"])
def test_viewer_keeps_all_models_hydrogens_and_altlocs(live_window, tmp_path, suffix):
    st = gemmi.Structure()
    for mi in range(2):
        model, chain = gemmi.Model(mi + 1), gemmi.Chain("A")
        for code in (" ", "A"):
            res = gemmi.Residue()
            res.name, res.seqid = "ALA", gemmi.SeqId(25, code)
            for atom_name, element, alt in (("CA", "C", "A"), ("CA", "C", "B"), ("H", "H", "\0")):
                atom = gemmi.Atom()
                atom.name, atom.element, atom.altloc = atom_name, gemmi.Element(element), alt
                atom.pos = gemmi.Position(mi * 10, len(res), ord(code))
                res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
        st.add_model(model)
    src = tmp_path / ("ensemble" + suffix)
    if suffix == ".cif":
        st.make_mmcif_document().write_file(str(src))
    else:
        st.write_minimal_pdb(str(src))
    c, viewer = live_window.controller, live_window.viewer
    c.load_structure(src)
    assert wait_until(lambda: bool(js(viewer, "objects.ensemble && objects.ensemble.editable")), 10)
    assert js(viewer, "objects.ensemble.model.selectedAtoms({}).length") == 12
    js(viewer, "setSelectionMode('residue'); pickSelection('ensemble', objects.ensemble.model.selectedAtoms({})[0], {}); true")
    assert wait_until(lambda: len(c.selection.get("ensemble", [])) == 3)
    assert {r["model"] for r in c.selection["ensemble"]} == {0}
    assert {r["icode"] for r in c.selection["ensemble"]} == {""}
    c.delete_selected_residues()
    assert len(atom_records(gemmi.read_structure(str(c.source_pdb)))) == 9
    assert c.source_pdb.suffix == suffix


def test_delete_everything_can_be_undone_and_reopened(qapp, settings, tmp_path):
    w = MainWindow(settings, use_viewer=False)
    c = w.controller
    c.set_project(Project.create(tmp_path / "project"))
    try:
        c.load_structure(DATA / "3PTB.pdb")
        original = c.source_pdb
        records = atom_records(gemmi.read_structure(str(original)))
        c.set_selection({original.stem: records}, "all atoms")
        c.delete_selected_residues()
        assert not atom_records(gemmi.read_structure(str(c.source_pdb)))
        c.undo_structure_edit()
        assert c.source_pdb == original
        c.redo_structure_edit()
        c.load_project_state()
        assert not atom_records(gemmi.read_structure(str(c.source_pdb)))
    finally:
        w.close()
