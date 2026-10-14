"""Structured application logging.

All application modules log through children of the ``mds`` logger
(``logging.getLogger("mds.jobs")`` etc.).  Records go to a rotating file
in the user's application-data folder and to the console; the GUI log
panel attaches its own handler on top of this.
"""

from __future__ import annotations

import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path
from types import TracebackType

LOGGER_NAME = "mds"
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

log = logging.getLogger(LOGGER_NAME)


def configure_logging(log_dir: Path, console_level: int = logging.INFO) -> Path:
    """Attach file and console handlers to the ``mds`` logger.

    Safe to call more than once; handlers are only added the first time.
    Returns the path of the active log file.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "mds.log"

    log.setLevel(logging.DEBUG)
    log.propagate = False
    if any(isinstance(h, RotatingFileHandler) for h in log.handlers):
        return log_file

    formatter = logging.Formatter(LOG_FORMAT)

    file_handler = RotatingFileHandler(
        log_file, maxBytes=2_000_000, backupCount=5, encoding="utf-8"
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    log.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(console_level)
    console_handler.setFormatter(formatter)
    log.addHandler(console_handler)

    return log_file


def install_exception_hook() -> None:
    """Log uncaught exceptions instead of letting them vanish inside Qt slots."""

    def _hook(
        exc_type: type[BaseException],
        exc: BaseException,
        tb: TracebackType | None,
    ) -> None:
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc, tb)
            return
        log.critical("Uncaught exception", exc_info=(exc_type, exc, tb))

    sys.excepthook = _hook
