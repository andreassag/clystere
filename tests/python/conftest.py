"""Add the pipeline's bin/ directory to sys.path so tests can import the scripts directly."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "bin"))
