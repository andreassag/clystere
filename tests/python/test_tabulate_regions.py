"""Unit tests for bin/tabulate_regions.py.

Tests cover the pure-Python functions:
    - extract_top_known_cluster_hit
    - extract_known_cluster_blast_results
    - parse_antismash_json
"""

import json

import pytest
from tabulate_regions import (
    extract_known_cluster_blast_results,
    extract_top_known_cluster_hit,
    parse_antismash_json,
)

# ── Helpers ──────────────────────────────────────────────────────────────────

NO_HIT = {"knownclusterblast_hit": "", "knownclusterblast_accession": "", "knownclusterblast_similarity": ""}


def _make_ranking_entry(similarity: int) -> list:
    """Build a minimal KnownClusterBlast ranking list for one region."""
    return [
        [
            {"description": "Epothilone D", "accession": "BGC0000168"},
            {"similarity": similarity},
        ]
    ]


def _make_kcb_results(similarity: int, region_index: int = 0) -> list[dict]:
    """Build a KnownClusterBlast results list with a hit at ``region_index``."""
    results = [{"ranking": []}] * (region_index + 1)
    results[region_index] = {"ranking": _make_ranking_entry(similarity)}
    return results


# ── extract_top_known_cluster_hit ─────────────────────────────────────────────


def test_no_hit_when_results_are_none():
    assert extract_top_known_cluster_hit(None, 0) == NO_HIT


def test_no_hit_when_results_list_is_empty():
    assert extract_top_known_cluster_hit([], 0) == NO_HIT


def test_no_hit_when_ranking_is_empty():
    assert extract_top_known_cluster_hit([{"ranking": []}], 0) == NO_HIT


def test_no_hit_at_similarity_threshold_boundary():
    """Similarity of exactly 15 % is below the reporting cutoff."""
    assert extract_top_known_cluster_hit(_make_kcb_results(15), 0) == NO_HIT


def test_no_hit_below_threshold():
    assert extract_top_known_cluster_hit(_make_kcb_results(10), 0) == NO_HIT


@pytest.mark.parametrize("similarity,expected_category", [(16, "low"), (51, "medium"), (76, "high")])
def test_similarity_categories(similarity: int, expected_category: str):
    result = extract_top_known_cluster_hit(_make_kcb_results(similarity), 0)
    assert result["knownclusterblast_similarity"] == expected_category


def test_hit_fields_are_populated():
    result = extract_top_known_cluster_hit(_make_kcb_results(80), 0)
    assert result["knownclusterblast_hit"] == "Epothilone D"
    assert result["knownclusterblast_accession"] == "BGC0000168"


def test_correct_region_selected_by_index():
    """Results for region 0 are a no-hit; region 1 has a real hit."""
    results = [{"ranking": []}, {"ranking": _make_ranking_entry(80)}]
    assert extract_top_known_cluster_hit(results, 0) == NO_HIT
    assert extract_top_known_cluster_hit(results, 1)["knownclusterblast_similarity"] == "high"


# ── extract_known_cluster_blast_results ──────────────────────────────────────


def test_extract_returns_none_for_empty_modules():
    assert extract_known_cluster_blast_results({"modules": {}}) is None


def test_extract_returns_none_when_modules_key_missing():
    assert extract_known_cluster_blast_results({}) is None


def test_extract_returns_none_for_missing_clusterblast_module():
    record = {"modules": {"antismash.modules.nrps_pks": {}}}
    assert extract_known_cluster_blast_results(record) is None


def test_extract_returns_results_list():
    sentinel = [{"ranking": []}]
    record = {"modules": {"antismash.modules.clusterblast": {"knowncluster": {"results": sentinel}}}}
    assert extract_known_cluster_blast_results(record) is sentinel


# ── parse_antismash_json ──────────────────────────────────────────────────────


def _minimal_antismash_json(tmp_path, areas: list, features: list, modules: dict | None = None) -> list[dict]:
    """Write a minimal antiSMASH JSON fixture and parse it; return the region rows."""
    data = {
        "input_file": "GCF_fixture.gbff",
        "records": [
            {
                "name": "CM001234.1",
                "description": "Streptomyces sp. chromosome",
                "areas": areas,
                "features": features,
                "modules": modules or {},
            }
        ],
    }
    json_path = tmp_path / "GCF_fixture" / "GCF_fixture.json"
    json_path.parent.mkdir()
    json_path.write_text(json.dumps(data))
    return parse_antismash_json(json_path)


def test_parse_produces_one_row_per_region(tmp_path):
    features = [
        {
            "type": "region",
            "location": "[0:45000](+)",
            "qualifiers": {"region_number": ["1"], "contig_edge": ["False"], "product": ["T1PKS"]},
        }
    ]
    rows = _minimal_antismash_json(tmp_path, areas=[{"products": ["T1PKS"]}], features=features)
    assert len(rows) == 1


def test_parse_row_fields(tmp_path):
    features = [
        {
            "type": "region",
            "location": "[0:45000](+)",
            "qualifiers": {"region_number": ["1"], "contig_edge": ["False"], "product": ["T1PKS"]},
        }
    ]
    rows = _minimal_antismash_json(tmp_path, areas=[{"products": ["T1PKS"]}], features=features)
    row = rows[0]
    assert row["file"] == "GCF_fixture"
    assert row["record_id"] == "CM001234.1"
    assert row["region"] == "1"
    assert row["start"] == "0"
    assert row["end"] == "45000"
    assert row["contig_edge"] == "False"
    assert row["product"] == "T1PKS"


def test_parse_hybrid_product_joined_with_slash(tmp_path):
    features = [
        {
            "type": "region",
            "location": "[0:45000](+)",
            "qualifiers": {"region_number": ["1"], "contig_edge": ["False"], "product": ["NRPS", "T1PKS"]},
        }
    ]
    rows = _minimal_antismash_json(tmp_path, areas=[{"products": ["NRPS", "T1PKS"]}], features=features)
    assert rows[0]["product"] == "NRPS / T1PKS"


def test_parse_skips_record_with_no_areas(tmp_path):
    data = {
        "input_file": "GCF_empty.gbff",
        "records": [
            {
                "name": "CM001234.1",
                "description": "Empty contig",
                "areas": [],
                "features": [],
                "modules": {},
            }
        ],
    }
    json_path = tmp_path / "GCF_empty" / "GCF_empty.json"
    json_path.parent.mkdir()
    json_path.write_text(json.dumps(data))
    rows = parse_antismash_json(json_path)
    assert rows == []
