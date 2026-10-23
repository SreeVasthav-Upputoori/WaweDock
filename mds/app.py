"""Guarded application bootstrap; importing this module needs no native libraries."""
from __future__ import annotations
import logging
from pathlib import Path
import sys
from mds.core.startup import prepare_startup, report_startup_failure


def main(argv: list[str] | None = None) -> int:
    arguments=list(sys.argv if argv is None else argv)
    args=arguments[1:]
    check_startup=args in (['--startup-check'],['--startup-check-full'])
    require_full=args==['--startup-check-full']
    log_dir=None
    try:
        # This must happen before Qt, RDKit, or MainWindow imports. A pythonw
        # import failure used to exit silently before logging was configured.
        log_dir=prepare_startup()
        if len(args)==2 and args[0]=='--probe-chemistry':
            from mds.core.capabilities import write_probe
            return write_probe(args[1])
        if args in (["--register-file-associations"], ["--unregister-file-associations"]):
            from mds.windows_integration import register, unregister
            action=unregister if args[0].startswith('--unregister') else register
            action()
            return 0
        if args[:1]==['--validate']:
            return _validate(args[1:])
        options={'check_startup':check_startup}
        if require_full: options['require_full']=True
        return _run_gui(arguments,[] if check_startup else args,log_dir,**options)
    except Exception as exc:
        report_startup_failure(exc,log_dir,show_dialog=not check_startup and args[:1] not in (['--probe-chemistry'],['--validate']))
        return 1


def _run_gui(arguments, args, log_dir, *, check_startup=False, require_full=False):
    import tempfile
    import time
    from PySide6.QtCore import QCoreApplication, QSettings, QTimer, Qt
    from PySide6.QtGui import QIcon, QPixmap
    from PySide6.QtWidgets import QApplication, QSplashScreen
    try:  # QtWebEngine must load before QApplication exists.
        from PySide6 import QtWebEngineWidgets  # noqa: F401
    except ImportError:
        QtWebEngineWidgets=None
    from mds import APP_NAME, ORG_NAME, __version__
    from mds.core.paths import assets_dir
    # Keep the native dependency import inside the startup error boundary.
    from mds.gui.main_window import MainWindow

    if sys.platform=='win32':
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('WaweDock.Desktop')
    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    QApplication.setOrganizationName(ORG_NAME)
    QApplication.setApplicationName(APP_NAME)
    QApplication.setApplicationVersion(__version__)
    app=QApplication(arguments[:1]+args)
    icon=QIcon(str(assets_dir()/'binding-pocket.ico'))
    app.setWindowIcon(icon)
    if sys.platform=='win32' and not check_startup:
        from mds.windows_integration import register
        try: register()
        except OSError: logging.getLogger('mds').warning('Windows Open with registration failed',exc_info=True)
    splash=QSplashScreen(QPixmap(str(assets_dir()/'wavedock-splash.png')).scaledToWidth(280,Qt.TransformationMode.SmoothTransformation))
    splash.setWindowIcon(icon);splash.show();app.processEvents()
    # A release smoke check must not modify the user's stored layout or projects.
    temporary=tempfile.TemporaryDirectory(prefix='wavedock-startup-') if check_startup else None
    settings=QSettings(str(Path(temporary.name)/'settings.ini'),QSettings.Format.IniFormat) if temporary else QSettings()
    window=None
    try:
        window=MainWindow(settings,log_dir=log_dir)
        window.show();splash.finish(window)
        if check_startup:
            deadline=time.monotonic()+20
            timer=QTimer(window);timer.setInterval(100)
            def check():
                if window.isVisible() and window.viewer is not None and window.viewer._ready and window.controller.chemistry is not None:
                    import json
                    from dataclasses import asdict
                    status=window.controller.chemistry
                    (log_dir/'startup-check.json').write_text(json.dumps({'core_ready':True,'chemistry':asdict(status),'full_required':require_full},indent=2),encoding='utf-8')
                    code=0 if status.available or not require_full else 1
                    logging.getLogger('mds').log(logging.INFO if code==0 else logging.ERROR,
                        'Startup check: core ready; chemistry available=%s; full required=%s',status.available,require_full)
                    app.exit(code)
                elif time.monotonic()>deadline:
                    logging.getLogger('mds').error('Startup check failed: main window / 3D viewer not ready within 20 seconds')
                    app.exit(1)
            timer.timeout.connect(check);timer.start()
        elif args:
            QTimer.singleShot(0,lambda:window.open_files([Path(arg) for arg in args]))
        return app.exec()
    finally:
        splash.close()
        if window is not None and window.isVisible():
            window.close()
        if temporary: temporary.cleanup()


def _validate(args: list[str]) -> int:
    """Headless redocking benchmark (see mds.validation); writes validation.json."""
    import json
    from mds.validation import redock
    if len(args)!=4:
        print('usage: --validate <structure.pdb> <RESNAME> <SMILES> <output_dir>')
        return 2
    out=Path(args[3]);out.mkdir(parents=True,exist_ok=True)
    try: result=redock(Path(args[0]),args[1],args[2],out)
    except Exception as exc:
        (out/'validation_error.txt').write_text(f'{type(exc).__name__}: {exc}')
        print(f'FAILED: {exc}');return 1
    print(json.dumps(result,indent=2))
    return 0 if result['success_top_pose'] else 3


if __name__=='__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    sys.exit(main())
