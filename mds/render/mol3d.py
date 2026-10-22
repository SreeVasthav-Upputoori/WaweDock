"""3Dmol.js viewer embedded through QtWebEngine (implements ``ViewerBackend``).

Python -> JavaScript: ``runJavaScript("fn(args)")`` with JSON-encoded arguments.
JavaScript -> Python: the page writes ``console.log("WD:" + json)`` and
``_Page.javaScriptConsoleMessage`` parses it (no QWebChannel needed).
Calls made before the page reports ``ready`` are queued and replayed in order.
"""

from __future__ import annotations

import base64
import json
import logging
from collections.abc import Callable

from PySide6.QtCore import QPoint, QTimer, QUrl, Signal
from PySide6.QtGui import QColor
from PySide6.QtWebEngineCore import QWebEnginePage
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QWidget

from mds.render.context_menu import build_context_menu

from mds.core.paths import assets_dir

log = logging.getLogger("mds.render")


class _Page(QWebEnginePage):
    message = Signal(str, object)

    def javaScriptConsoleMessage(self, level, message, line, source) -> None:  # noqa: N802
        if message.startswith("WD:"):
            try:
                payload = json.loads(message[3:])
                self.message.emit(payload["type"], payload.get("data"))
            except (ValueError, KeyError):
                log.debug("Bad viewer message: %s", message)
        elif level != QWebEnginePage.JavaScriptConsoleMessageLevel.InfoMessageLevel:
            log.warning("Viewer JS: %s (line %d)", message, line)


class Mol3DViewer(QWebEngineView):
    atom_info = Signal(str)
    measurement = Signal(str)
    selection_changed = Signal(dict, str)       # detail {object: [{chain, resi}, …]}, summary
    ready = Signal()
    edit_requested = Signal(str)
    health_inspect_requested = Signal()
    visibility_changed = Signal(str, bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._context_menu = None
        self._ready = False
        self._queue: list[str] = []
        self._resize_timer = QTimer(self)
        self._resize_timer.setSingleShot(True)
        self._resize_timer.setInterval(80)
        self._resize_timer.timeout.connect(self._sync_viewport)
        self._page = _Page(self)
        self._page.setBackgroundColor(QColor("black"))
        self.setPage(self._page)
        self._page.message.connect(self._on_message)
        # Protocol attributes (callbacks) map onto the Qt signals.
        self.on_atom_info: Callable[[str], None] = self.atom_info.emit
        self.on_measurement: Callable[[str], None] = self.measurement.emit
        self.on_selection: Callable[[dict, str], None] = self.selection_changed.emit
        self.load(QUrl.fromLocalFile(str(assets_dir() / "viewer" / "viewer.html")))

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if hasattr(self, "_resize_timer"):
            self._resize_timer.start()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._resize_timer.start()

    def _sync_viewport(self) -> None:
        if self._ready and self.isVisible() and self.width() > 0 and self.height() > 0:
            self._page.runJavaScript("requestViewportSync();")

    # -------------------------------------------------------------- plumbing
    def _on_message(self, kind: str, data: object) -> None:
        if kind == "ready":
            self._ready = True
            for js in self._queue:
                self._page.runJavaScript(js)
            self._queue.clear()
            self._resize_timer.start()
            self.ready.emit()
        elif kind == "background":
            color = QColor(str(data))
            if color.isValid():
                # Leave the Chromium callback before changing its native surface.
                QTimer.singleShot(0, lambda value=color: self._page.setBackgroundColor(value))
        elif kind == "jserror":
            log.error("Viewer script error: %s", data)
        elif kind == "atom":
            self.on_atom_info(str(data))
        elif kind == "measure":
            self.on_measurement(str(data))
        elif kind == "selection":
            payload = data if isinstance(data, dict) else {}
            self.on_selection(payload.get("detail", {}) or {}, str(payload.get("summary", "")))
        elif kind == "edit" and data in ("delete", "undo", "redo"):
            self.edit_requested.emit(str(data))
        elif kind == "context" and isinstance(data, dict):
            # Leave the Chromium console callback before opening native UI.
            # A nested menu.exec() here can block the JavaScript bridge.
            QTimer.singleShot(0, lambda state=dict(data): self._show_context_menu(state))

    def _show_context_menu(self, state: dict) -> None:
        if self._context_menu is not None:
            self._context_menu.close()
        menu = build_context_menu(self, state)
        self._context_menu = menu
        menu.aboutToHide.connect(lambda: setattr(self, "_context_menu", None))
        menu.aboutToHide.connect(menu.deleteLater)
        menu.popup(self.mapToGlobal(QPoint(int(state.get("x", 0)), int(state.get("y", 0)))))

    def _call(self, fn: str, *args: object) -> None:
        js = f"{fn}({', '.join(json.dumps(a) for a in args)});"
        if self._ready:
            self._page.runJavaScript(js)
        else:
            self._queue.append(js)

    # ------------------------------------------------------------ ViewerBackend
    def add_object(self, name: str, data: str, fmt: str, kind: str = "protein", identities=None,
                   preserve_view: bool = False, preserve_style: bool = False) -> None:
        self._call("addObject", name, data, fmt, kind, identities, preserve_view, preserve_style)

    def copy_appearance(self, source: str, target: str) -> None:
        self._call("copyObjectAppearance", source, target)

    def replace_structure(self, old_name: str, name: str, data: str, fmt: str, identities=None) -> None:
        self._call("replaceStructure", old_name, name, data, fmt, identities)

    def remove_object(self, name: str) -> None:
        self._call("removeObject", name)

    def set_visible(self, name: str, visible: bool) -> None:
        self._call("setVisible", name, visible)
        self.visibility_changed.emit(name, visible)

    def set_representation(self, name: str, representation: str) -> None:
        self._call("setRepresentation", name, representation)

    def set_color_scheme(self, name: str, scheme: str) -> None:
        self._call("setColorScheme", name, scheme)

    def zoom_to(self, name: str | None = None) -> None:
        self._call("zoomTo", name)

    def show_residues_around(self, ligand: str, cutoff: float, protein: str) -> None:
        self._call("showAround", ligand, cutoff, protein)

    def clear_residue_selection(self) -> None:
        self._call("clearAround")

    def show_box(self, center, size, status="ok") -> None:
        self._call("showBox", list(center), list(size), status)

    def clear_box(self) -> None:
        self._call("clearBox")

    def set_measure_mode(self, mode: str) -> None:
        self._call("setMeasureMode", mode)

    def clear_measurements(self) -> None:
        self._call("clearMeasurements")

    def label_residues(self, enabled: bool) -> None:
        self._call("labelResidues", enabled)

    def set_background(self, dark: bool) -> None:
        self._call("setBackground", dark)

    # Selection ---------------------------------------------------------------
    def clear_selection(self) -> None:
        self._call("clearSelection")

    def set_health_issues(self, name: str, issues) -> None:
        colors = {"Error": "#ff4444", "Warning": "#ff9933", "Note": "#33aaff"}
        markers, seen = [], set()
        # Bound geometry independently of report size. Keep exact editing styles intact.
        for issue in sorted(issues, key=lambda i: {"Error":0,"Warning":1,"Note":2}.get(i.severity,3)):
            for atom in issue.atoms:
                key = (atom["model"], atom["chain_index"], atom["residue_index"])
                if key in seen or len(markers) >= 250: continue
                seen.add(key)
                markers.append({"atom_index": atom["atom_index"], "color": colors.get(issue.severity,"#33aaff")})
        self._call("setHealthIssues", name, markers)

    def set_selection_mode(self, mode: str) -> None:
        self._call("setSelectionMode", mode)

    def select_within_radius(self, radius: float) -> None:
        self._call("selectWithinRadius", radius)

    def expand_selection_to_chain(self) -> None:
        self._call("expandSelectionToChain")

    def invert_selection(self) -> None:
        self._call("invertSelection")

    def set_selection(self, detail: dict) -> None:
        """Overwrite the viewer's selection with ``detail`` (object -> list of residues)."""
        self._call("setSelectionFromPython", detail)

    def export_image(self, path: str, scale: int = 3) -> None:
        """Render at ``scale`` x the on-screen resolution and save a PNG."""
        def _save(uri: object) -> None:
            if not isinstance(uri, str) or "," not in uri:
                log.error("Image export failed: viewer returned no image")
                return
            with open(path, "wb") as fh:
                fh.write(base64.b64decode(uri.split(",", 1)[1]))
            log.info("Image saved: %s", path)

        self._page.runJavaScript(f"exportPng({int(scale)})", _save)
