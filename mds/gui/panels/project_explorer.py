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
    reference_requested = Signal(str)   # use this ligand as the pose-review reference
    objects_changed = Signal()                   # object was added or removed

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

    def kind_of(self, obj: str) -> str:
        """Return the stored ``kind`` of ``obj`` (``"protein"``, ``"ligand"``, …) or ``""``."""
        item = self.find(obj)
        return (item.data(KindRole) or "") if item is not None else ""

    def is_visible(self, obj: str) -> bool:
        item = self.find(obj)
        return bool(item is not None and item.checkState() == Qt.CheckState.Checked)

    def all_object_names(self, kinds: tuple[str, ...] | None = None) -> list[str]:
        """Every object name in the workspace, optionally restricted to given kinds."""
        out: list[str] = []
        for cat in self._categories.values():
            for r in range(cat.rowCount()):
                item = cat.child(r)
                if kinds is None or (item.data(KindRole) or "") in kinds:
                    out.append(item.data(ObjectRole))
        return out

    def upsert_object(self, category: str, name: str, props: dict, *, visible: bool = True,
                      kind: str = "protein") -> None:
        item = self.find(name)
        self._updating = True
        added = item is None
        if added:
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
        if added:
            self.objects_changed.emit()

    def remove_object(self, name: str) -> None:
        item = self.find(name)
        if item is not None:
            item.parent().removeRow(item.row())
            self.objects_changed.emit()

    def clear_category(self, category: str) -> None:
        cat = self._categories[category]
        if cat.rowCount():
            cat.removeRows(0, cat.rowCount())
            self.objects_changed.emit()

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
        item = self.model.itemFromIndex(index) if index.isValid() else None
        self.build_context_menu(item).exec(self.view.viewport().mapToGlobal(pos))

    def build_context_menu(self, item: QStandardItem | None) -> QMenu:
        """The menu for ``item`` (``None`` for empty space); separated so it can be inspected."""
        menu = QMenu(self)
        obj = item.data(ObjectRole) if item is not None else None
        if obj is not None:
            menu.addAction("Zoom to", lambda: self.zoom_requested.emit(obj))
            rep = menu.addMenu("Representation")
            for r in REPRESENTATIONS:
                rep.addAction(r.title(), lambda r=r: self.representation_requested.emit(obj, r))
            col = menu.addMenu("Colour by")
            for c in COLOR_SCHEMES:
                col.addAction(c.title(), lambda c=c: self.color_requested.emit(obj, c))
            if (item.data(KindRole) or "") == "ligand":
                menu.addSeparator()
                action = menu.addAction("Use as pose reference",
                                        lambda: self.reference_requested.emit(obj))
                action.setToolTip("Keep this ligand visible while browsing docking poses. "
                                  "Its coordinates are not changed.")
            menu.addSeparator()
            menu.addAction("Remove from workspace", lambda: self.remove_requested.emit(obj))
            menu.addSeparator()
        menu.addAction("Expand All", self.view.expandAll)
        menu.addAction("Collapse All", self.view.collapseAll)
        return menu
