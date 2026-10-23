"""Results workspace navigation, pose identity, comparison, charts and exports."""

import pytest
from PySide6.QtCore import QItemSelectionModel, QSettings, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QDockWidget
Chem = pytest.importorskip("rdkit.Chem", reason="Requires the optional RDKit backend", exc_type=ImportError)

from mds.chem.workflow import LigandOutcome
from mds.engines.base import Pose
from mds.gui.controller import Controller
from mds.gui.main_window import MainWindow
from mds.gui.panels.results_panel import ResultsPanel
from mds.jobs import JobRunner
from mds.project import Project


def sample_outcomes():
    molecule = Chem.MolFromSmiles('CO')
    conformer = Chem.Conformer(2)
    conformer.SetAtomPosition(0, (0, 0, 0))
    conformer.SetAtomPosition(1, (1.4, 0, 0))
    molecule.AddConformer(conformer)
    a = LigandOutcome('ligand A', 'ok', poses=[Pose(1, -7, 0, 0, ''), Pose(2, -6, 1, 2, '')],
                      mols=[Chem.Mol(molecule), Chem.Mol(molecule)], reference_rmsd=[None, 1.5])
    b = LigandOutcome('ligand B', 'ok', poses=[Pose(1, -9, 0, 0, '')],
                      mols=[Chem.Mol(molecule)], reference_rmsd=[None])
    return [a, b]


@pytest.fixture
def panel(qapp, tmp_path, monkeypatch):
    controller = Controller(JobRunner())
    controller.set_project(Project.create(tmp_path / 'project'))
    monkeypatch.setattr(controller, 'contacts_for', lambda *_args: [])
    widget = ResultsPanel(controller)
    controller.outcomes = sample_outcomes()
    widget.show_outcomes(controller.outcomes)
    yield widget
    widget.close()


def selected_pose(panel):
    return panel.interactions_panel._current


def test_search_and_numeric_sort_keep_actual_pose_identity(panel):
    assert panel.table.rowCount() == 3
    assert selected_pose(panel)[0].name == 'ligand B'
    assert float(panel.table.item(0, 2).text()) == -9
    panel.search.setText('LIGAND A')
    assert panel.match_count.text() == '2 / 3 poses'
    assert selected_pose(panel)[0].name == 'ligand A'
    panel.table.sortItems(2, Qt.SortOrder.DescendingOrder)
    panel.select_pose(panel.c.outcomes[0], 1)
    assert selected_pose(panel) == (panel.c.outcomes[0], 1)
    panel.search.setText('no matches')
    assert selected_pose(panel) is None
    assert panel.comparison.rowCount() == 0
    panel.search.clear()
    assert selected_pose(panel) is not None


def test_missing_reference_rmsd_does_not_drop_poses(panel):
    outcomes = sample_outcomes()
    outcomes[0].reference_rmsd = []
    panel.show_outcomes(outcomes)
    assert panel.table.rowCount() == 3
    assert len(panel._rows) == 3


def test_multiple_selected_poses_have_comparable_metrics(panel):
    panel.table.selectRow(0)
    panel.table.selectionModel().select(panel.table.model().index(1, 0),
        QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows)
    assert panel.comparison.rowCount() == 2
    assert {panel.comparison.item(i, 0).text() for i in range(2)} == {'ligand A', 'ligand B'}
    assert {panel.comparison.item(i, 2).text() for i in range(2)} == {'-9.000', '-7.000'}


def test_clicking_score_chart_selects_best_pose_and_clears_filter(panel, qapp):
    panel.resize(750, 850)
    panel.show()
    panel.analysis_tabs.setCurrentWidget(panel.score_chart)
    panel.search.setText('ligand B')
    chart = panel.score_chart
    assert len(chart.bars) == 2
    bar = chart.bars[-1]
    chart.ensureVisible(bar)
    qapp.processEvents()
    point = chart.mapFromScene(bar.sceneBoundingRect().center())
    QTest.mouseClick(chart.viewport(), Qt.MouseButton.LeftButton, pos=point)
    assert panel.search.text() == ''
    assert selected_pose(panel) == (bar.outcome, bar.pose_index)


def test_report_exports_all_results_and_escapes_names(panel, tmp_path):
    panel.c.outcomes[0].name = '<ligand & A>'
    panel.c.outcomes.append(LigandOutcome('<failed>', 'failed', '<message>'))
    panel.show_outcomes(panel.c.outcomes)
    panel.search.setText('ligand B')
    path = tmp_path / 'report.html'
    panel.export_report(path)
    report = path.read_text(encoding='utf-8')
    assert '&lt;ligand &amp; A&gt;' in report
    assert '&lt;failed&gt;' in report and '&lt;message&gt;' in report
    assert '-7' in report and '-9' in report
    assert '<ligand & A>' not in report


def test_csv_excel_and_sdf_export_real_pose_data(panel, tmp_path):
    import csv
    import openpyxl
    csv_path, xlsx_path, sdf_path = [tmp_path / f'results.{suffix}' for suffix in ('csv', 'xlsx', 'sdf')]
    panel.c.export_csv(csv_path)
    panel.c.export_excel(xlsx_path)
    panel.c.export_sdf(sdf_path)
    with csv_path.open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    assert {float(row['Score (kcal/mol)']) for row in rows} == {-9, -7, -6}
    workbook = openpyxl.load_workbook(xlsx_path, read_only=True)
    assert workbook['Docking results'].max_row == 4
    workbook.close()
    molecules = list(Chem.SDMolSupplier(str(sdf_path)))
    assert len(molecules) == 3
    assert {float(m.GetProp('Score_kcal_per_mol')) for m in molecules} == {-9, -7, -6}


@pytest.mark.parametrize('fmt,method', [('csv', 'export_csv'), ('xlsx', 'export_excel'), ('sdf', 'export_sdf')])
def test_export_menu_preserves_existing_exporters(panel, tmp_path, monkeypatch, fmt, method):
    path = tmp_path / f'results.{fmt}'
    calls = []
    monkeypatch.setattr(panel.c, method, lambda target: calls.append(target))
    monkeypatch.setattr('mds.gui.panels.results_panel.QFileDialog.getSaveFileName', lambda *_args: (str(path), ''))
    panel._export(fmt)
    assert calls == [path]


def test_interaction_job_failure_allows_retry_and_selection_reset(panel):
    interactions = panel.interactions_panel
    interactions._started()
    assert not interactions.analyse_btn.isEnabled()
    panel.c.runner._names[5] = 'Interaction analysis (ProLIF)'
    panel.c.runner.job_failed.emit(5, 'DLL unavailable')
    assert interactions.analyse_btn.isEnabled()
    assert 'failed' in interactions.status.text()
    panel.show_outcomes([])
    assert interactions._current is None
    assert panel.contacts.rowCount() == 0
    assert not panel.score_chart.bars


def test_workspace_replaces_results_dock_and_restores_setup(qapp, settings):
    window = MainWindow(settings, use_viewer=False)
    try:
        assert not window.findChildren(QDockWidget, 'dockResults')
        viewport = window.viewport
        window.dock_selection.show()
        window.workspace.select(2)
        assert window.centralWidget() is window.workspace
        assert not window.results_panel.isHidden()
        assert window.dock_docking.isHidden()
        assert window.dock_selection.isHidden()
        window.workspace.select(0)
        assert window.results_panel.isHidden()
        assert not window.dock_selection.isHidden()
        assert not window.dock_docking.isHidden()
        assert window.workspace.splitter.widget(1) is viewport
        window.workspace.select(1)
        assert window.docking_panel.current_step == 3
        window.actions_['analyse.results'].trigger()
        assert window.workspace.tabs.currentIndex() == 2
        window.actions_['view.fullscreen'].setChecked(True)
        assert window.results_panel.isHidden() and window.workspace.tabs.isHidden()
        window.actions_['view.exit_fullscreen'].trigger()
        assert not window.results_panel.isHidden()
        assert not window.workspace.tabs.isHidden()
    finally:
        window.close()


def test_session_and_named_layout_restore_workspace(qapp, tmp_path):
    path = str(tmp_path / 'workspace.ini')
    first = MainWindow(QSettings(path, QSettings.Format.IniFormat), use_viewer=False)
    first.dock_selection.show()
    first.workspace.select(2)
    first.layouts.save_named('analysis')
    first.workspace.select(0)
    first.layouts.restore_named('analysis')
    assert first.workspace.tabs.currentIndex() == 2
    first.close()
    second = MainWindow(QSettings(path, QSettings.Format.IniFormat), use_viewer=False)
    try:
        assert second.workspace.tabs.currentIndex() == 2
        assert not second.results_panel.isHidden()
        assert second.dock_docking.isHidden()
        second.workspace.select(0)
        assert not second.dock_selection.isHidden()
    finally:
        second.close()


def test_review_navigation_follows_filtered_sorted_pose_rows(panel):
    bar=panel.review_bar
    assert bar.position.text()=='1 / 3'
    bar.next.click()
    assert panel._current_pose[0].name=='ligand A'
    assert bar.position.text()=='2 / 3'
    panel.search.setText('ligand A')
    panel.table.sortItems(2,Qt.SortOrder.DescendingOrder)
    panel.select_pose(panel.c.outcomes[0],1)
    assert bar.position.text().endswith('/ 3')  # select_pose clears the search
    panel.search.setText('ligand A')
    panel.table.selectRow(next(r for r in range(panel.table.rowCount()) if not panel.table.isRowHidden(r)))
    assert bar.position.text()=='1 / 2'
    bar.next.click()
    assert bar.position.text()=='2 / 2'
    assert not bar.next.isEnabled()
    bar.previous.click()
    assert bar.position.text()=='1 / 2'
