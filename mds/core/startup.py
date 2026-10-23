"""Dependency-free startup diagnostics, including failures before Qt can load."""
from __future__ import annotations
import faulthandler
import logging
import os
from pathlib import Path
import sys
import tempfile
import traceback
from datetime import datetime, timezone
from mds import APP_NAME, ORG_NAME, __version__
from mds.core.logging_setup import configure_logging, install_exception_hook

_crash_stream = None


def default_log_dir() -> Path:
    """Match the application-data location without importing Qt or chemistry DLLs."""
    if sys.platform=='win32':
        base=Path(os.environ.get('LOCALAPPDATA') or Path.home()/'AppData'/'Local')/ORG_NAME/APP_NAME
    elif sys.platform=='darwin':
        base=Path.home()/'Library'/'Application Support'/ORG_NAME/APP_NAME
    else:
        base=Path(os.environ.get('XDG_DATA_HOME') or Path.home()/'.local'/'share')/ORG_NAME/APP_NAME
    return base/'logs'


def prepare_startup() -> Path:
    global _crash_stream
    candidates=(default_log_dir(),Path(tempfile.gettempdir())/APP_NAME/'logs')
    last_error=None
    for directory in candidates:
        try:
            directory.mkdir(parents=True,exist_ok=True)
            # pythonw and frozen GUI executables have no stderr/stdout. Keep
            # warnings and backend error messages instead of discarding them.
            for name in ('stdout','stderr'):
                if getattr(sys,name) is None:
                    setattr(sys,name,(directory/'startup-console.log').open('a',encoding='utf-8',buffering=1))
            configure_logging(directory)
            install_exception_hook()
            try:
                if _crash_stream is None:
                    _crash_stream=(directory/'native-crash.log').open('a',encoding='utf-8',buffering=1)
                faulthandler.enable(file=_crash_stream,all_threads=True)
            except (OSError,RuntimeError):
                logging.getLogger('mds').warning('Native crash tracing could not be enabled',exc_info=True)
            logging.getLogger('mds').info('Starting %s %s | Python %s | executable %s | frozen=%s',
                APP_NAME,__version__,sys.version.split()[0],sys.executable,bool(getattr(sys,'frozen',False)))
            return directory
        except OSError as exc:
            last_error=exc
    raise OSError('Could not initialise startup logging in application data or the temporary folder.') from last_error


def failure_message(exc: Exception, report: Path | None) -> str:
    reason=f'{type(exc).__name__}: {exc}'
    blocked=any(word in str(exc).lower() for word in ('application control','code integrity','signing level'))
    if blocked:
        explanation=('Windows Application Control blocked a required dependency. '
                     'WaweDock cannot load that library while the block is active. '
                     'A trusted, properly signed build is needed; security settings were not changed.')
    else:
        explanation='A startup error prevented the application from opening. The details below identify the failure.'
    return f'{APP_NAME} could not start.\n\n{explanation}\n\n{reason[:1200]}\n\n'+(
        f'Diagnostic report:\n{report}' if report else 'The diagnostic report could not be written.')


def show_startup_error(message: str) -> None:
    """Use a native Windows dialog even if PySide6 itself failed to import."""
    if sys.platform=='win32':
        import ctypes
        ctypes.windll.user32.MessageBoxW(None,message,f'{APP_NAME} - startup failed',0x10)
    else:
        stream=sys.stderr or sys.__stderr__
        if stream: print(message,file=stream)


def report_startup_failure(exc: Exception, directory: Path | None, *, show_dialog=True) -> Path | None:
    report=None
    text=(f'{APP_NAME} {__version__} startup failure\nUTC: {datetime.now(timezone.utc).isoformat()}\n'
          f'Python: {sys.version}\nExecutable: {sys.executable}\nFrozen: {bool(getattr(sys,"frozen",False))}\n\n'+
          ''.join(traceback.format_exception(type(exc),exc,exc.__traceback__)))
    for folder in (directory,Path(tempfile.gettempdir())/APP_NAME/'logs'):
        if folder is None: continue
        try:
            folder.mkdir(parents=True,exist_ok=True)
            target=folder/'startup-error.txt'
            target.write_text(text,encoding='utf-8')
            report=target
            break
        except OSError: continue
    logging.getLogger('mds').critical('Startup failed; diagnostic report: %s',report,
                                     exc_info=(type(exc),exc,exc.__traceback__))
    message=failure_message(exc,report)
    if show_dialog:
        try: show_startup_error(message)
        except Exception:
            logging.getLogger('mds').error('Could not display the startup error dialog',exc_info=True)
    else:
        stream=sys.stderr or sys.__stderr__
        if stream: print(message,file=stream)
    return report
