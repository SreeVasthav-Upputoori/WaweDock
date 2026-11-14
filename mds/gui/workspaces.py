"""Workspace navigation around one persistent molecular viewer."""
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QLabel, QSplitter, QTabBar, QVBoxLayout, QWidget

STAGES = (
    ("1  Prepare", "Load a protein, prepare the receptor, and add ligands."),
    ("2  Dock", "Choose the binding site, review the box, and run docking."),
    ("3  Review", "Browse poses with a locked reference and inspect their analysis."),
)


class WorkspaceHost(QWidget):
    workspace_changed = Signal(int)

    def __init__(self, viewport, parent=None):
        super().__init__(parent)
        self.tabs = QTabBar()
        self.tabs.setExpanding(False)
        self.tabs.setDocumentMode(True)
        self.tabs.setDrawBase(False)
        for title, _description in STAGES:
            self.tabs.addTab(title)
        self.stage_title = QLabel()
        self.stage_title.setObjectName("workspaceStageTitle")
        self.stage_title.setStyleSheet("font-size: 15px; font-weight: 600;")
        self.stage_description = QLabel()
        self.stage_description.setObjectName("workspaceStageDescription")
        self.stage_description.setStyleSheet("color: palette(mid);")
        self.stage_description.setWordWrap(True)
        intro = QVBoxLayout()
        intro.setContentsMargins(8, 4, 8, 0)
        intro.setSpacing(1)
        intro.addWidget(self.stage_title)
        intro.addWidget(self.stage_description)
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
        layout.addLayout(intro)
        layout.addWidget(self.splitter, 1)
        self.tabs.currentChanged.connect(self._changed)
        self._update_stage_text(0)

    def set_results(self, results):
        self._results = results
        self.splitter.insertWidget(0, results)
        self.splitter.setStretchFactor(0, 2)
        self.splitter.setStretchFactor(1, 3)
        self.splitter.setSizes([550, 850])
        results.hide()

    def _changed(self, index):
        self._update_stage_text(index)
        self._update_visibility()
        self.workspace_changed.emit(index)

    def _update_stage_text(self, index):
        title, description = STAGES[index]
        self.stage_title.setText(title)
        self.stage_description.setText(description)

    def _update_visibility(self):
        self.tabs.setVisible(not self._focus)
        if self._results is not None:
            self._results.setVisible(self.tabs.currentIndex() == 2 and not self._focus)

    def set_focus_mode(self, enabled):
        self._focus = enabled
        self._update_visibility()

    def select(self, index):
        self.tabs.setCurrentIndex(index)
