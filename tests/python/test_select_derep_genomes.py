"""Unit tests for bin/select_derep_genomes.py."""

import csv

import checkm2_to_drep_info as ctd
import pytest
import select_derep_genomes as sdg


def _write_fasta(path, seq_lens):
    lines = []
    for idx, length in enumerate(seq_lens, start=1):
        lines.append(f">ctg{idx}\n" + ("A" * length) + "\n")
    path.write_text("".join(lines))


def test_skani_sparse_clustering_and_n50_selection(tmp_path):
    g1 = tmp_path / "mag1.fna"
    g2 = tmp_path / "mag2.fna"
    g3 = tmp_path / "mag3.fna"
    _write_fasta(g1, [5000, 5000])  # N50 = 5000
    _write_fasta(g2, [10000])  # N50 = 10000 (should win over mag1)
    _write_fasta(g3, [4000])  # singleton (not present in sparse skani edges)

    skani_tsv = tmp_path / "skani.tsv"
    skani_tsv.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "mag1.fna\tmag2.fna\t99.5\t95.0\t96.0\tctg1\tctg1\n"
    )

    out_ids = tmp_path / "ids.txt"
    out_clusters = tmp_path / "clusters.tsv"
    rc = sdg.main(
        [
            "--skani",
            str(skani_tsv),
            "--genomes",
            str(g1),
            str(g2),
            str(g3),
            "--ani",
            "99.0",
            "--out-ids",
            str(out_ids),
            "--out-clusters",
            str(out_clusters),
        ]
    )
    assert rc == 0
    assert out_ids.read_text().splitlines() == ["mag2", "mag3"]

    with out_clusters.open() as fh:
        rows = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["mag2"]["is_representative"] == "true"
    assert rows["mag1"]["is_representative"] == "false"
    assert rows["mag3"]["is_representative"] == "true"
    assert rows["mag1"]["cluster_id"] == rows["mag2"]["cluster_id"]


def test_skani_dense_matrix_and_checkm2_weighting(tmp_path):
    g1 = tmp_path / "magA.fna"
    g2 = tmp_path / "magB.fna"
    _write_fasta(g1, [10000])
    _write_fasta(g2, [5000])

    # magB has lower N50 but much higher completeness & lower contamination -> must win
    checkm2 = tmp_path / "checkm2.tsv"
    checkm2.write_text("Name\tCompleteness\tContamination\nmagA\t75.0\t8.0\nmagB\t99.0\t0.5\n")

    skani_dense = tmp_path / "skani_dense.tsv"
    skani_dense.write_text("2\nmagA.fna\nmagB.fna\t99.4\n")

    out_ids = tmp_path / "ids.txt"
    out_clusters = tmp_path / "clusters.tsv"
    rc = sdg.main(
        [
            "--skani",
            str(skani_dense),
            "--genomes",
            str(g1),
            str(g2),
            "--checkm2",
            str(checkm2),
            "--ani",
            "99.0",
            "--out-ids",
            str(out_ids),
            "--out-clusters",
            str(out_clusters),
        ]
    )
    assert rc == 0
    assert out_ids.read_text().splitlines() == ["magB"]


def test_checkm2_to_drep_info_matches_staged_fastas(tmp_path):
    g1 = tmp_path / "mag1.fna"
    g2 = tmp_path / "mag2.fasta"
    _write_fasta(g1, [1000])
    _write_fasta(g2, [1000])

    checkm2_tsv = tmp_path / "checkm2.tsv"
    checkm2_tsv.write_text("Name\tCompleteness\tContamination\nmag1\t98.5\t1.2\nmag2.fasta\t84.0\t3.5\n")

    out_csv = tmp_path / "drep_genome_info.csv"
    rc = ctd.main(["--checkm2", str(checkm2_tsv), "--genomes", str(g1), str(g2), "--out", str(out_csv)])
    assert rc == 0
    raw = out_csv.read_bytes()
    assert b"\r\n" not in raw

    with out_csv.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert rows == [
        {"genome": "mag1.fna", "completeness": "98.50", "contamination": "1.20"},
        {"genome": "mag2.fasta", "completeness": "84.00", "contamination": "3.50"},
    ]


def test_checkm2_to_drep_info_missing_genome_raises(tmp_path):
    g1 = tmp_path / "mag1.fna"
    g2 = tmp_path / "mag_missing.fna"
    _write_fasta(g1, [1000])
    _write_fasta(g2, [1000])

    checkm2_tsv = tmp_path / "checkm2.tsv"
    checkm2_tsv.write_text("Name\tCompleteness\tContamination\nmag1\t98.5\t1.2\n")

    out_csv = tmp_path / "drep_genome_info.csv"
    with pytest.raises(ValueError, match="mag_missing.fna"):
        ctd.main(["--checkm2", str(checkm2_tsv), "--genomes", str(g1), str(g2), "--out", str(out_csv)])
