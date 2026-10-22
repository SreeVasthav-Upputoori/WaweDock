from unittest.mock import Mock

import pytest
from PySide6.QtCore import QObject, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QAbstractButton

from mds.gui.panels.selection_panel import SelectionPanel


@pytest.fixture(autouse=True)
def sidebar_font(qapp):
    original = qapp.font()
    qapp.setFont(QFont("Segoe UI", 9))
    yield
    qapp.setFont(original)


class Controller(QObject):
    selection_changed = Signal(dict, str)
    structure_loaded = Signal(str)
    edit_history_changed = Signal(bool, bool)
    editing_view_requested = Signal()
    selection_mode_requested = Signal(str)
    selection_radius_requested = Signal(float)

    def __init__(self):
        super().__init__()
        self.selection = {}
        for name in ('undo_structure_edit', 'redo_structure_edit', 'cleanup_structure',
                     'expand_selection_to_chain', 'invert_selection', 'clear_selection',
                     'keep_only_selected_residues', 'delete_selected_residues'):
            setattr(self, name, Mock())


@pytest.mark.parametrize('width', [200, 240, 320])
def test_sidebar_fits_without_horizontal_clipping(qapp, width):
    panel = SelectionPanel(Controller())
    panel.resize(width, 740)
    panel.show()
    for _ in range(3):
        qapp.processEvents()
    try:
        assert panel.scroll.horizontalScrollBar().maximum() == 0
        viewport = panel.scroll.viewport()
        for button in panel.findChildren(QAbstractButton):
            if button.isVisible():
                left = button.mapTo(viewport, button.rect().topLeft()).x()
                assert left >= 0
                assert left + button.width() <= viewport.width()
                assert button.fontMetrics().horizontalAdvance(button.text()) + 24 <= button.width()
    finally:
        panel.close()


def test_selection_actions_and_cleanup_still_work(qapp):
    controller = Controller()
    panel = SelectionPanel(controller)
    assert not panel.delete_btn.isEnabled()
    controller.selection = {'protein': [1]}
    controller.selection_changed.emit(controller.selection, '1 atom selected')
    panel.delete_btn.click()
    controller.delete_selected_residues.assert_called_once()
    assert panel.radius_btn.isEnabled()
    panel.remove_water.setChecked(True)
    panel.remove_chains.setText(' B, C ')
    panel.cleanup_btn.click()
    controller.cleanup_structure.assert_called_once_with(
        waters=True, hetero=False, ligands=False, remove_chains=['B', 'C'], retain_selected=True)
    controller.edit_history_changed.emit(True, False)
    assert panel.undo_btn.isEnabled()
    assert not panel.redo_btn.isEnabled()
    panel.close()
