"""Shared I/O helpers for antiSMASH JSON processing scripts."""

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
from typing import Callable, TypeVar

T = TypeVar("T")


def list_antismash_json_files(directory: Path) -> list[Path]:
    """Return sorted antiSMASH JSON files under <sample>/<sample>.json layout."""
    return sorted(directory.glob("*/*.json"))


def load_antismash_json(json_path: Path) -> dict:
    """Load and decode one antiSMASH JSON file."""
    with json_path.open() as json_file:
        return json.load(json_file)


def run_threaded(paths: list[Path], worker: Callable[[Path], T], threads: int | None = None) -> list[T]:
    """Execute a path-based worker function in parallel and collect results."""
    if not paths:
        return []

    results: list[T] = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(worker, path): path for path in paths}
        for future in as_completed(futures):
            results.append(future.result())

    return results
