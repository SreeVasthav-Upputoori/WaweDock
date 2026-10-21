"""Workspace navigation around one persistent molecular viewer."""
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QSplitter, QTabBar, QVBoxLayout, QWidget


class WorkspaceHost(QWidget):
    workspace_changed = Signal(int)

    def __init__(self, viewport, parent=None):
        super().__init__(parent)
        self.tabs = QTabBar()
        self.tabs.setExpanding(False)
        for title in ('Prepare', 'Dock', 'Results & Analysis'):
            self.tabs.addTab(title)
        self.splitter = QSplitter(Qt.Orientation.Horizontal)
        self.splitter.setObjectName('resultsWorkspaceSplit')
        self.splitter.setChildrenCollapsible(False)
        self.splitter.addWidget(viewport)
        self._results = None
        self._focus = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addWidget(self.tabs)
        layout.addWidget(self.splitter, 1)
        self.tabs.currentChanged.connect(self._changed)

    def set_results(self, results):
        self._results = results
        self.splitter.insertWidget(0, results)
        self.splitter.setStretchFactor(0, 2)
        self.splitter.setStretchFactor(1, 3)
        self.splitter.setSizes([550, 850])
        results.hide()

    def _changed(self, index):
        self._update_visibility()
        self.workspace_changed.emit(index)

    def _update_visibility(self):
        self.tabs.setVisible(not self._focus)
        if self._results is not None:
            self._results.setVisible(self.tabs.currentIndex() == 2 and not self._focus)

    def set_focus_mode(self, enabled):
        self._focus = enabled
        self._update_visibility()

    def select(self, index):
        self.tabs.setCurrentIndex(index)
