"""Unit tests for bin/extract_ripp_candidates.py."""

import csv
import json

import extract_ripp_candidates


def test_extract_ripp_modules_and_rrefinder(tmp_path):
    sample_json = tmp_path / "genomeA.json"
    payload = {
        "input_file": "genomeA.gbff",
        "records": [
            {
                "name": "ctg1",
                "features": [
                    {
                        "type": "region",
                        "location": "[1000:15000](+)",
                        "qualifiers": {"region_number": ["1"], "product": ["lanthipeptide"]},
                    },
                    {
                        "type": "CDS",
                        "location": "[2000:2200](+)",
                        "qualifiers": {"locus_tag": ["LOC_0001"]},
                    },
                    {
                        "type": "CDS",
                        "location": "[5000:5400](+)",
                        "qualifiers": {"locus_tag": ["LOC_0002"]},
                    },
                    {
                        "type": "CDS_motif",
                        "location": "[3000:3150](+)",
                        "qualifiers": {
                            "locus_tag": ["LOC_0003"],
                            "prepeptide": ["lassopeptide"],
                            "leader_sequence": ["MKKQ"],
                            "core_sequence": ["GGAGHVPEYFVGIGTPISFYG"],
                            "note": ["macrolactam"],
                        },
                    },
                ],
                "modules": {
                    "antismash.modules.rrefinder": {
                        "hits_by_cds": {
                            "LOC_0001": [{"domain": "Lan_RRE", "bitscore": 42.6}],
                            "LOC_0002": [{"domain": "Lasso_F_RRE", "bitscore": 31.2}],
                        }
                    },
                    "antismash.modules.lanthipeptides": {
                        "motifs": {
                            "LOC_0001": [
                                {
                                    "locus_tag": "LOC_0001",
                                    "location": "[2000:2180](+)",
                                    "leader": "MSTKDFNLDLV",
                                    "core": "ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK",
                                    "detailed_info": {"predicted_class": ["Class I"], "bridges": ["5"]},
                                }
                            ]
                        }
                    },
                },
            }
        ],
    }
    sample_json.write_text(json.dumps(payload))

    out_tsv = tmp_path / "ripp_candidates.tsv"
    rc = extract_ripp_candidates.main(tmp_path, out_tsv)
    assert rc == 0

    with out_tsv.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert len(rows) == 3
    by_locus = {row["locus_tag"]: row for row in rows}
    assert by_locus["LOC_0001"]["ripp_subclass"] == "lanthipeptide"
    assert by_locus["LOC_0001"]["region"] == "1"
    assert by_locus["LOC_0001"]["rre_bitscore"] == "42.6"
    assert by_locus["LOC_0001"]["core_sequence"].startswith("ITSISL")
    assert float(by_locus["LOC_0001"]["predicted_molecular_weight"]) > 3000.0
    assert "predicted_class=Class I" in by_locus["LOC_0001"]["tailoring_modifications"]

    assert by_locus["LOC_0003"]["ripp_subclass"] == "lassopeptide"
    assert by_locus["LOC_0003"]["core_sequence"] == "GGAGHVPEYFVGIGTPISFYG"

    assert by_locus["LOC_0002"]["ripp_subclass"] == "RRE:Lasso_F_RRE"
    assert by_locus["LOC_0002"]["rre_bitscore"] == "31.2"


def test_extract_ripp_empty_directory_writes_header_only(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    out_tsv = tmp_path / "ripp_candidates.tsv"
    assert extract_ripp_candidates.main(empty_dir, out_tsv) == 0
    lines = out_tsv.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t") == extract_ripp_candidates.OUTPUT_COLS


def test_extract_ripp_zero_default_molecular_weight_falls_back_to_estimate(tmp_path):
    """antiSMASH 8 Prepeptide.to_json() emits 0.0 when weight is unset; must estimate from core sequence."""
    sample_json = tmp_path / "genomeB.json"
    payload = {
        "records": [
            {
                "name": "ctg1",
                "features": [],
                "modules": {
                    "antismash.modules.lassopeptides": {
                        "motifs": {
                            "LOC_0010": [
                                {
                                    "locus_tag": "LOC_0010",
                                    "leader": "MKKQ",
                                    "core": "GGAGHVPEYFVGIGTPISFYG",
                                    "molecular_weight": 0.0,
                                    "monoisotopic_mass": 0.0,
                                }
                            ]
                        }
                    }
                },
            }
        ]
    }
    sample_json.write_text(json.dumps(payload))
    out_tsv = tmp_path / "ripp.tsv"
    assert extract_ripp_candidates.main(tmp_path, out_tsv) == 0
    with out_tsv.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert float(rows[0]["predicted_molecular_weight"]) > 1500.0
