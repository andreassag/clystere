"""Shared I/O helpers for antiSMASH JSON processing scripts."""

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from typing import TypeVar

T = TypeVar("T")


def _iter_json_files(directory: Path) -> list[Path]:
    found: list[Path] = []
    for entry in directory.iterdir():
        if entry.name.startswith("."):
            continue
        if entry.is_dir():
            found.extend(_iter_json_files(entry))
        elif entry.suffix == ".json":
            found.append(entry)
    return found


def list_antismash_json_files(directory: Path) -> list[Path]:
    """Return sorted antiSMASH JSON files in flat or <sample>/<sample>.json layout."""
    return sorted(_iter_json_files(directory))


def load_antismash_json(json_path: Path) -> dict:
    """Load and decode one antiSMASH JSON file."""
    with json_path.open() as json_file:
        return json.load(json_file)


def run_threaded(paths: list[Path], worker: Callable[[Path], T], threads: int | None = None) -> list[T]:
    """Execute a path-based worker function in parallel and preserve input order."""
    if not paths:
        return []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        return list(executor.map(worker, paths))
