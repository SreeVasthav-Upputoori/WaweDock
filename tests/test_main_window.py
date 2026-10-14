from __future__ import annotations

import logging

import pytest
from conftest import wait_until
from PySide6.QtCore import QSettings

from mds.gui.layouts import LayoutManager
from mds.gui.main_window import MainWindow
from mds.gui.theme import Theme
from mds.jobs import JobState
from mds.jobs.diagnostics import JobSystemSelfTest


@pytest.fixture
def window(qapp, settings):
    w = MainWindow(settings, use_viewer=False)   # headless tests skip the WebEngine viewer
    yield w
    w.close()


def test_menus_and_docks_exist(window):
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    assert titles == ["File", "Edit", "View", "Build", "Prepare", "Docking",
                      "Analyse", "Tools", "Window", "Help"]
    for dock in (window.dock_project, window.dock_docking, window.dock_properties, window.dock_results,
                 window.dock_report, window.dock_log, window.dock_jobs):
        assert dock.objectName()
    assert window.windowTitle() == "WaweDock"


def test_working_actions_enabled_and_future_ones_labelled(window):
    assert window.actions_["docking.run"].isEnabled()
    assert window.actions_["file.open_structure"].isEnabled()
    screen = window.actions_["docking.screen"]
    assert not screen.isEnabled() and "Phase 9" in screen.toolTip()


def test_log_panel_receives_records(window):
    logging.getLogger("mds.test").warning("hello from test")
    assert wait_until(lambda: "hello from test" in window.log_panel._view.toPlainText())


def test_theme_switch_persists(window, settings):
    window.set_theme(Theme.DARK)
    assert settings.value("appearance/theme") == "dark"
    assert window.actions_["view.theme.dark"].isChecked()
    window.set_theme(Theme.LIGHT)
    assert window.actions_["view.theme.light"].isChecked()


def test_selection_updates_properties_and_status(window):
    window.project_explorer.view.setCurrentIndex(window.project_explorer.model.index(0, 0))
    assert window.properties._title.text() == "<b>Proteins</b>"


def test_named_layout_round_trip(window):
    window.dock_jobs.hide()
    window.layouts.save_named("minimal")
    window.dock_jobs.show()
    assert window.layouts.names() == ["minimal"]
    assert window.layouts.restore_named("minimal")
    assert window.dock_jobs.isHidden()
    window.layouts.delete_named("minimal")
    assert window.layouts.names() == []


def test_session_layout_restored_in_new_window(qapp, tmp_path):
    path = str(tmp_path / "s.ini")
    first = MainWindow(QSettings(path, QSettings.Format.IniFormat), use_viewer=False)
    first.dock_report.hide()
    first.close()
    second = MainWindow(QSettings(path, QSettings.Format.IniFormat), use_viewer=False)
    try:
        assert second.dock_report.isHidden()
    finally:
        second.close()


def test_self_test_job_appears_in_jobs_panel(window):
    job_id = window.job_runner.submit(JobSystemSelfTest(steps=3, interval_s=0.01))
    assert window.jobs_panel.row_count() == 1
    assert wait_until(lambda: window.job_runner.state(job_id) is JobState.FINISHED)


def test_layout_version_rejects_stale_state(window, settings):
    settings.setValue("session/state", b"garbage")
    assert LayoutManager(window, settings).restore_session() is False
