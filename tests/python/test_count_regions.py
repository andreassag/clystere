"""Unit tests for bin/count_regions.py.

Tests cover the pure-Python functions:
  - tabulate_bgc_counts
  - parse_antismash_json
"""

import json

from count_regions import parse_antismash_json, tabulate_bgc_counts

# ── Shared fixtures ───────────────────────────────────────────────────────────

# Two contigs: contig_1 has a solo T1PKS and a hybrid NRPS+T1PKS region;
# contig_2 has a solo RiPP-like region.
GENOME_PRODUCTS: dict[str, dict[str, list[list[str]]]] = {
    "GCF_A": {
        "contig_1": [["T1PKS"], ["NRPS", "T1PKS"]],
        "contig_2": [["RiPP-like"]],
    }
}

CONTIG_DESCRIPTIONS: dict[str, dict[str, str]] = {
    "GCF_A": {
        "contig_1": "Streptomyces sp. chromosome 1",
        "contig_2": "Streptomyces sp. chromosome 2",
    }
}


# ── tabulate_bgc_counts – per-genome (default) ────────────────────────────────


def test_per_genome_single_row():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    assert len(rows) == 1


def test_per_genome_total_count():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    assert rows[0]["total_count"] == 3  # T1PKS + hybrid(NRPS/T1PKS) + RiPP-like


def test_per_genome_hybrid_counted_once():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    assert rows[0]["hybrid"] == 1


def test_per_genome_solo_types_counted():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    assert rows[0]["T1PKS"] == 1  # solo T1PKS only; the hybrid is counted as "hybrid"
    assert rows[0]["RiPP-like"] == 1


def test_per_genome_record_is_genome_name():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    assert rows[0]["record"] == "GCF_A"


def test_per_genome_description_includes_contig_count():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    assert "2 total records" in str(rows[0]["description"])


# ── tabulate_bgc_counts – per-contig ─────────────────────────────────────────


def test_per_contig_two_rows():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS, per_contig=True)
    assert len(rows) == 2


def test_per_contig_record_names():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS, per_contig=True)
    records = {str(row["record"]) for row in rows}
    assert "GCF_A|contig_1" in records
    assert "GCF_A|contig_2" in records


def test_per_contig_counts():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS, per_contig=True)
    contig1 = next(row for row in rows if "contig_1" in str(row["record"]))
    assert contig1["total_count"] == 2
    assert contig1["hybrid"] == 1


# ── tabulate_bgc_counts – split_hybrids ──────────────────────────────────────


def test_split_hybrids_nrps_and_t1pks_both_counted():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS, split_hybrids=True)
    assert rows[0]["NRPS"] == 1
    assert rows[0]["T1PKS"] == 2  # solo + one from split hybrid


def test_split_hybrids_no_hybrid_column():
    rows = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS, split_hybrids=True)
    assert rows[0].get("hybrid", 0) == 0


def test_total_count_unaffected_by_split():
    per_genome = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS)
    split = tabulate_bgc_counts(GENOME_PRODUCTS, CONTIG_DESCRIPTIONS, split_hybrids=True)
    assert per_genome[0]["total_count"] == split[0]["total_count"]


# ── tabulate_bgc_counts – edge cases ─────────────────────────────────────────


def test_empty_genome_products_returns_empty():
    assert tabulate_bgc_counts({}, {}) == []


def test_genome_with_no_bgcs_produces_row():
    """A genome with contigs that have zero areas still appears with total_count=0."""
    genome_products = {"GCF_empty": {"contig_1": []}}
    contig_descriptions = {"GCF_empty": {"contig_1": "Empty contig"}}
    rows = tabulate_bgc_counts(genome_products, contig_descriptions)
    assert rows[0]["total_count"] == 0


# ── parse_antismash_json ──────────────────────────────────────────────────────


def test_parse_returns_correct_input_file(tmp_path):
    data = {
        "input_file": "GCF_fixture.gbff",
        "records": [{"name": "CM001234.1", "description": "Test contig", "areas": []}],
    }
    json_path = tmp_path / "GCF_fixture.json"
    json_path.write_text(json.dumps(data))

    input_file, _, _ = parse_antismash_json(json_path)
    assert input_file == "GCF_fixture.gbff"


def test_parse_products_by_contig(tmp_path):
    data = {
        "input_file": "GCF_fixture.gbff",
        "records": [
            {
                "name": "CM001234.1",
                "description": "Chromosome",
                "areas": [{"products": ["T1PKS"]}, {"products": ["NRPS", "T1PKS"]}],
            }
        ],
    }
    json_path = tmp_path / "GCF_fixture.json"
    json_path.write_text(json.dumps(data))

    _, products_by_contig, _ = parse_antismash_json(json_path)
    assert "CM001234.1" in products_by_contig
    assert products_by_contig["CM001234.1"] == [["T1PKS"], ["NRPS", "T1PKS"]]


def test_parse_contig_descriptions(tmp_path):
    data = {
        "input_file": "GCF_fixture.gbff",
        "records": [{"name": "CM001234.1", "description": "Streptomyces sp.", "areas": []}],
    }
    json_path = tmp_path / "GCF_fixture.json"
    json_path.write_text(json.dumps(data))

    _, _, descriptions = parse_antismash_json(json_path)
    assert descriptions["CM001234.1"] == "Streptomyces sp."


def test_parse_missing_description_defaults_to_empty(tmp_path):
    data = {
        "input_file": "GCF_fixture.gbff",
        "records": [{"name": "CM001234.1", "areas": []}],
    }
    json_path = tmp_path / "GCF_fixture.json"
    json_path.write_text(json.dumps(data))

    _, _, descriptions = parse_antismash_json(json_path)
    assert descriptions["CM001234.1"] == ""
