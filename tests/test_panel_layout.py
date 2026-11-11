"""The analysis panels stay compact, and the caveats stay readable while doing so.

Height budgets are the point of this file. The controls used to take so much of the results
panel that the metrics table was left about one row tall, so these bounds are deliberately
generous against the current layout and still far below what the old one needed.
"""

import pytest

from mds.gui.panels.common import CaveatBar, ElidedLabel
from mds.gui.panels.pose_review import NO_REFERENCE_HINT, PoseReviewBar

LONG = ("A pharmacophore derived here describes one docked or crystal pose, which is a "
        "hypothesis about what a ligand needs at those positions and not a validated model.")


def test_elided_label_shortens_what_it_paints_but_not_what_it_reports(qapp):
    label = ElidedLabel(LONG)
    label.resize(120, 16)
    try:
        assert label.text() == LONG, "text() must round-trip, whatever the width"
        assert label.toolTip() == LONG
        assert len(label.displayedText()) < len(LONG)
        assert label.sizeHint().height() < 40, "a status line must never wrap"
    finally:
        label.close()


def test_caveat_stays_on_screen_in_one_line_with_the_full_wording_reachable(qapp):
    bar = CaveatBar(LONG)
    bar.resize(300, 20)
    bar.show()
    try:
        assert bar.isVisible(), "the caveat is shortened, never hidden"
        assert bar.text() == LONG
        assert bar.toolTip() == LONG and bar.label.toolTip() == LONG
        assert bar.sizeHint().height() <= 40
        assert len(bar.label.displayedText()) < len(LONG)
    finally:
        bar.close()


def test_review_bar_is_one_row_once_a_reference_exists(qapp):
    bar = PoseReviewBar()
    bar.show()
    try:
        empty = bar.sizeHint().height()
        assert bar.hint.isVisible() and bar.hint.text() == NO_REFERENCE_HINT
        bar.set_references(["STI (crystal)"])
        bar.layout().activate()           # sizeHint is cached until the layout is re-run
        assert not bar.hint.isVisible()
        assert bar.sizeHint().height() < empty, "the hint's row collapses once it is unneeded"
        assert bar.sizeHint().height() <= 40, "the review controls are a single row"
    finally:
        bar.close()


def test_display_toggles_are_in_the_menu_and_still_drive_the_viewer(qapp):
    """Moving them off the face of the bar must not change what they do."""
    bar = PoseReviewBar()
    changes = []
    bar.settings_changed.connect(lambda: changes.append(True))
    try:
        assert bar.display.menu().actions() == [bar.keep_camera, bar.hide_protein]
        bar.set_references(["ref"])
        changes.clear()
        bar.hide_protein.setChecked(True)
        assert changes, "hiding the protein must still ask the viewer to update"
        changes.clear()
        bar.lock_reference.setChecked(True)
        assert changes and bar.keep_camera.isChecked() and not bar.keep_camera.isEnabled()
    finally:
        bar.close()


def test_results_focus_review_hides_analysis_chrome(results):
    results.focus_review.setChecked(True)
    assert results.analysis_tabs.isHidden()
    assert results.focus_review.text() == "Show analysis"
    results.focus_review.setChecked(False)
    assert not results.analysis_tabs.isHidden()
    assert results.focus_review.text() == "Focus review"


@pytest.fixture
def results(qapp, tmp_path):
    from mds.gui.controller import Controller
    from mds.gui.panels.results_panel import ResultsPanel
    from mds.jobs import JobRunner
    from mds.project import Project

    controller = Controller(JobRunner())
    controller.set_project(Project.create(tmp_path / "project"))
    panel = ResultsPanel(controller)
    yield panel
    panel.close()


def test_the_results_panel_leaves_room_for_its_tables(results):
    """Chrome above and below the tables has a budget; the tables get the rest."""
    results.resize(980, 760)
    results.show()
    assert results.minimumSizeHint().height() <= 480
    assert results.consistency_panel.sizeHint().height() <= 320
    # Every tab fits on the bar, so none of them hides behind a scroll arrow.
    assert results.analysis_tabs.count() == 6
    assert all(len(results.analysis_tabs.tabText(i)) <= 14
               for i in range(results.analysis_tabs.count()))


def test_the_consistency_detail_pane_appears_only_when_it_has_something_to_say(results):
    tab = results.consistency_panel
    # isHidden, not isVisible: the panel itself has no parent window in this test.
    assert tab.detail.isHidden()
    tab._show_detail(_stub_result())
    assert not tab.detail.isHidden() and "pose A" in tab.detail.toPlainText()
    tab.clear()
    assert tab.detail.isHidden()


def _stub_result():
    from mds.chem.pose_consistency import PoseConsistency

    return PoseConsistency(label="pose A", consensus=0.8, verdict="Conserved",
                           notes=("a note",))


def test_the_axes_button_names_what_is_switched_off(results):
    """A disabled axis must be visible from the closed button, not only inside the menu."""
    tab = results.consistency_panel
    assert tab.axes_btn.text() == "All axes"
    tab.shape.setChecked(False)
    assert tab.axes_btn.text() == "Axes (2/3)"
    assert "Shape and volume overlap" in tab.axes_btn.toolTip()
    tab.shape.setChecked(True)
    assert tab.axes_btn.text() == "All axes"
