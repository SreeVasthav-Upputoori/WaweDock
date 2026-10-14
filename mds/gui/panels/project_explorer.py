"""Project Explorer: every molecular object in the workspace, with visibility checkboxes.

Model/view architecture (``QTreeView`` + ``QStandardItemModel``).  Object items carry
``ObjectRole`` = the viewer object name; category folders carry none.
"""

from __future__ import annotations

from PySide6.QtCore import QModelIndex, QPoint, Qt, Signal
from PySide6.QtGui import QStandardItem, QStandardItemModel
from PySide6.QtWidgets import QMenu, QStyle, QTreeView, QVBoxLayout, QWidget

from mds.render.base import COLOR_SCHEMES, REPRESENTATIONS

CATEGORIES = ("Proteins", "Ligands", "Docking Boxes", "Docking Poses")
ObjectRole = Qt.ItemDataRole.UserRole + 2
PropertiesRole = Qt.ItemDataRole.UserRole + 1
KindRole = Qt.ItemDataRole.UserRole + 3


class ProjectExplorer(QWidget):
    item_selected = Signal(str, dict)            # title, properties
    visibility_changed = Signal(str, bool)       # object name, visible
    representation_requested = Signal(str, str)
    color_requested = Signal(str, str)
    zoom_requested = Signal(str)
    remove_requested = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.model = QStandardItemModel(self)
        self.model.setHorizontalHeaderLabels(["Workspace"])
        icon = self.style().standardIcon(QStyle.StandardPixmap.SP_DirIcon)
        self._categories: dict[str, QStandardItem] = {}
        for name in CATEGORIES:
            item = QStandardItem(icon, name)
            item.setEditable(False)
            item.setData({"Type": "Category"}, PropertiesRole)
            self.model.appendRow(item)
            self._categories[name] = item
        self._updating = False
        self.model.itemChanged.connect(self._item_changed)

        self.view = QTreeView()
        self.view.setModel(self.model)
        self.view.setHeaderHidden(True)
        self.view.expandAll()
        self.view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.view.customContextMenuRequested.connect(self._context_menu)
        self.view.selectionModel().currentChanged.connect(self._current_changed)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.view)

    # ---------------------------------------------------------------- content
    def category(self, name: str) -> QStandardItem:
        return self._categories[name]

    def find(self, obj: str) -> QStandardItem | None:
        for cat in self._categories.values():
            for r in range(cat.rowCount()):
                if cat.child(r).data(ObjectRole) == obj:
                    return cat.child(r)
        return None

    def upsert_object(self, category: str, name: str, props: dict, *, visible: bool = True,
                      kind: str = "protein") -> None:
        item = self.find(name)
        self._updating = True
        if item is None:
            item = QStandardItem(name)
            item.setEditable(False)
            item.setCheckable(True)
            item.setData(name, ObjectRole)
            item.setData(kind, KindRole)
            self._categories[category].appendRow(item)
        item.setData(props, PropertiesRole)
        item.setCheckState(Qt.CheckState.Checked if visible else Qt.CheckState.Unchecked)
        self._updating = False
        self.view.expandAll()

    def remove_object(self, name: str) -> None:
        item = self.find(name)
        if item is not None:
            item.parent().removeRow(item.row())

    def clear_category(self, category: str) -> None:
        self._categories[category].removeRows(0, self._categories[category].rowCount())

    def set_visible(self, name: str, visible: bool) -> None:
        item = self.find(name)
        if item is not None:
            self._updating = True
            item.setCheckState(Qt.CheckState.Checked if visible else Qt.CheckState.Unchecked)
            self._updating = False

    def object_names(self, category: str) -> list[str]:
        cat = self._categories[category]
        return [cat.child(r).data(ObjectRole) for r in range(cat.rowCount())]

    # ----------------------------------------------------------------- events
    def _item_changed(self, item: QStandardItem) -> None:
        if self._updating or item.data(ObjectRole) is None:
            return
        self.visibility_changed.emit(item.data(ObjectRole), item.checkState() == Qt.CheckState.Checked)

    def _current_changed(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if not current.isValid():
            self.item_selected.emit("", {})
            return
        item = self.model.itemFromIndex(current)
        props = dict(item.data(PropertiesRole) or {})
        if props.get("Type") == "Category":
            props["Items"] = item.rowCount()
        self.item_selected.emit(item.text(), props)

    def _context_menu(self, pos: QPoint) -> None:
        index = self.view.indexAt(pos)
        menu = QMenu(self)
        item = self.model.itemFromIndex(index) if index.isValid() else None
        obj = item.data(ObjectRole) if item is not None else None
        if obj is not None:
            menu.addAction("Zoom to", lambda: self.zoom_requested.emit(obj))
            rep = menu.addMenu("Representation")
            for r in REPRESENTATIONS:
                rep.addAction(r.title(), lambda r=r: self.representation_requested.emit(obj, r))
            col = menu.addMenu("Colour by")
            for c in COLOR_SCHEMES:
                col.addAction(c.title(), lambda c=c: self.color_requested.emit(obj, c))
            menu.addSeparator()
            menu.addAction("Remove from workspace", lambda: self.remove_requested.emit(obj))
            menu.addSeparator()
        menu.addAction("Expand All", self.view.expandAll)
        menu.addAction("Collapse All", self.view.collapseAll)
        menu.exec(self.view.viewport().mapToGlobal(pos))
