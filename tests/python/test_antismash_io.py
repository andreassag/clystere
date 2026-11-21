"""Unit tests for bin/antismash_io.py."""

import json
import time

from antismash_io import list_antismash_json_files, load_antismash_json, run_threaded


def test_list_antismash_json_files_nested_and_flat(tmp_path):
    nested_dir = tmp_path / "sample_b"
    nested_dir.mkdir()
    nested_json = nested_dir / "sample_b.json"
    nested_json.write_text("{}")

    flat_json = tmp_path / "sample_a.json"
    flat_json.write_text("{}")

    hidden_dir = tmp_path / ".hidden"
    hidden_dir.mkdir()
    (hidden_dir / "ignored.json").write_text("{}")
    (tmp_path / ".dotfile.json").write_text("{}")

    found = list_antismash_json_files(tmp_path)
    assert found == [flat_json, nested_json]


def test_load_antismash_json(tmp_path):
    json_path = tmp_path / "sample.json"
    payload = {"input_file": "sample.gbff", "records": []}
    json_path.write_text(json.dumps(payload))

    assert load_antismash_json(json_path) == payload


def test_run_threaded_empty_returns_empty():
    assert run_threaded([], lambda p: p.name) == []


def test_run_threaded_preserves_input_order(tmp_path):
    paths = [tmp_path / f"item_{i}.json" for i in range(4)]

    def slow_first(path):
        idx = int(path.stem.split("_")[1])
        time.sleep((3 - idx) * 0.01)
        return path.stem

    assert run_threaded(paths, slow_first, threads=4) == ["item_0", "item_1", "item_2", "item_3"]
