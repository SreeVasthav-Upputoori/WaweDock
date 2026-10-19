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

from PySide6.QtCore import QUrl, Signal
from PySide6.QtWebEngineCore import QWebEnginePage
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QWidget

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

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ready = False
        self._queue: list[str] = []
        self._page = _Page(self)
        self.setPage(self._page)
        self._page.message.connect(self._on_message)
        # Protocol attributes (callbacks) map onto the Qt signals.
        self.on_atom_info: Callable[[str], None] = self.atom_info.emit
        self.on_measurement: Callable[[str], None] = self.measurement.emit
        self.on_selection: Callable[[dict, str], None] = self.selection_changed.emit
        self.load(QUrl.fromLocalFile(str(assets_dir() / "viewer" / "viewer.html")))

    # -------------------------------------------------------------- plumbing
    def _on_message(self, kind: str, data: object) -> None:
        if kind == "ready":
            self._ready = True
            for js in self._queue:
                self._page.runJavaScript(js)
            self._queue.clear()
            self.ready.emit()
        elif kind == "jserror":
            log.error("Viewer script error: %s", data)
        elif kind == "atom":
            self.on_atom_info(str(data))
        elif kind == "measure":
            self.on_measurement(str(data))
        elif kind == "selection":
            payload = data if isinstance(data, dict) else {}
            self.on_selection(payload.get("detail", {}) or {}, str(payload.get("summary", "")))

    def _call(self, fn: str, *args: object) -> None:
        js = f"{fn}({', '.join(json.dumps(a) for a in args)});"
        if self._ready:
            self._page.runJavaScript(js)
        else:
            self._queue.append(js)

    # ------------------------------------------------------------ ViewerBackend
    def add_object(self, name: str, data: str, fmt: str, kind: str = "protein") -> None:
        self._call("addObject", name, data, fmt, kind)

    def remove_object(self, name: str) -> None:
        self._call("removeObject", name)

    def set_visible(self, name: str, visible: bool) -> None:
        self._call("setVisible", name, visible)

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

    def show_box(self, center, size) -> None:
        self._call("showBox", list(center), list(size))

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
