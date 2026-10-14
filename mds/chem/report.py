"""Preparation report: an ordered list of every operation and warning."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Level(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True)
class ReportEntry:
    level: Level
    text: str


@dataclass
class Report:
    entries: list[ReportEntry] = field(default_factory=list)

    def info(self, text: str) -> None:
        self.entries.append(ReportEntry(Level.INFO, text))

    def warn(self, text: str) -> None:
        self.entries.append(ReportEntry(Level.WARNING, text))

    def error(self, text: str) -> None:
        self.entries.append(ReportEntry(Level.ERROR, text))

    @property
    def warnings(self) -> list[str]:
        return [e.text for e in self.entries if e.level is Level.WARNING]

    @property
    def errors(self) -> list[str]:
        return [e.text for e in self.entries if e.level is Level.ERROR]

    def render(self) -> str:
        tag = {Level.INFO: "INFO   ", Level.WARNING: "WARNING", Level.ERROR: "ERROR  "}
        return "\n".join(f"[{tag[e.level]}] {e.text}" for e in self.entries)


class PreparationError(Exception):
    """Raised when a molecule cannot be prepared; ``report`` explains why."""

    def __init__(self, message: str, report: Report | None = None) -> None:
        super().__init__(message)
        self.report = report or Report()
