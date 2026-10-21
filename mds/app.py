"""Application bootstrap: QApplication, logging, settings, main window."""

from __future__ import annotations

import os
import sys

# A windowed .exe has no console: give libraries that print a harmless sink.
for _name in ("stdout", "stderr"):
    if getattr(sys, _name) is None:
        setattr(sys, _name, open(os.devnull, "w", encoding="utf-8"))  # noqa: SIM115

from PySide6.QtCore import QCoreApplication, QSettings, QTimer, Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen

try:  # must be imported before QApplication exists (Qt requirement for QtWebEngine)
    from PySide6 import QtWebEngineWidgets  # noqa: F401
except ImportError:  # the app still runs, without the 3D viewer
    QtWebEngineWidgets = None

from mds import APP_NAME, ORG_NAME, __version__
from mds.core.logging_setup import configure_logging, install_exception_hook
from mds.core.paths import assets_dir
from mds.gui.main_window import MainWindow, app_icon, default_log_dir


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv[1:]
    if args in (["--register-file-associations"], ["--unregister-file-associations"]):
        from mds.windows_integration import register, unregister
        action = unregister if args[0].startswith("--unregister") else register
        action()
        return 0
    if args[:1] == ["--validate"]:
        return _validate(args[1:])
    # Give source launches their own Windows taskbar identity instead of Python's.
    if sys.platform == "win32":
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("WaweDock.Desktop")
    # Organisation/application names determine where QSettings and
    # QStandardPaths store data, so they must be set before either is used.
    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    QApplication.setOrganizationName(ORG_NAME)
    QApplication.setApplicationName(APP_NAME)
    QApplication.setApplicationVersion(__version__)
    app = QApplication(sys.argv if argv is None else argv)
    app.setWindowIcon(app_icon())

    log_dir = default_log_dir()
    configure_logging(log_dir)
    install_exception_hook()
    if sys.platform == "win32":
        from mds.windows_integration import register
        try:
            register()
        except OSError:
            import logging
            logging.getLogger("mds").warning("Windows Open with registration failed", exc_info=True)

    splash = QSplashScreen(QPixmap(str(assets_dir() / "wavedock-splash.png")).scaledToWidth(
        280, Qt.TransformationMode.SmoothTransformation))
    splash.setWindowIcon(app_icon())
    splash.show()
    app.processEvents()
    window = MainWindow(QSettings(), log_dir=log_dir)
    window.show()
    splash.finish(window)
    if args:
        from pathlib import Path
        QTimer.singleShot(0, lambda: window.open_files([Path(arg) for arg in args]))
    return app.exec()


def _validate(args: list[str]) -> int:
    """Headless redocking benchmark (see mds.validation); writes validation.json."""
    import json
    from pathlib import Path

    from mds.validation import redock

    if len(args) != 4:
        print("usage: --validate <structure.pdb> <RESNAME> <SMILES> <output_dir>")
        return 2
    out = Path(args[3])
    out.mkdir(parents=True, exist_ok=True)
    try:
        result = redock(Path(args[0]), args[1], args[2], out)
    except Exception as exc:  # noqa: BLE001
        (out / "validation_error.txt").write_text(f"{type(exc).__name__}: {exc}")
        print(f"FAILED: {exc}")
        return 1
    print(json.dumps(result, indent=2))
    return 0 if result["success_top_pose"] else 3


if __name__ == "__main__":
    sys.exit(main())
