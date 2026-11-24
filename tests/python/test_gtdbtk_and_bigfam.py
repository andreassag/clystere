"""Unit tests for bin/gtdbtk_to_bigslice_taxonomy.py and bin/extract_bigfam_scores.py."""

import csv
import sqlite3

import extract_bigfam_scores
import gtdbtk_to_bigslice_taxonomy


def test_gtdbtk_to_bigslice_taxonomy_parses_seven_ranks(tmp_path):
    bac_tsv = tmp_path / "gtdbtk.bac120.summary.tsv"
    bac_tsv.write_text(
        "user_genome\tclassification\n"
        "GCA_000005845.fna\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
        "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n"
        "MAG_02\td__Bacteria;p__Bacillota;c__Bacilli;o__;f__;g__;s__\n"
    )
    out_tax = tmp_path / "taxonomy.tsv"
    rc = gtdbtk_to_bigslice_taxonomy.main([str(bac_tsv), "-o", str(out_tax)])
    assert rc == 0

    lines = out_tax.read_text().splitlines()
    assert lines[0].startswith("# path_startswith\tKingdom\tPhylum")
    assert len(lines) == 3

    ecoli_cols = lines[1].split("\t")
    assert ecoli_cols == [
        "GCA_000005845/",
        "Bacteria",
        "Pseudomonadota",
        "Gammaproteobacteria",
        "Enterobacterales",
        "Enterobacteriaceae",
        "Escherichia",
        "Escherichia coli",
        "Escherichia coli",
    ]
    mag_cols = lines[2].split("\t")
    assert mag_cols[0] == "MAG_02/"
    assert mag_cols[2] == "Bacillota"
    assert mag_cols[8] == "MAG_02"


def test_gtdbtk_taxonomy_prevents_sample1_sample10_collision(tmp_path):
    bac_tsv = tmp_path / "gtdbtk.bac120.summary.tsv"
    bac_tsv.write_text(
        "user_genome\tclassification\n"
        "sample1.fna\td__Bacteria;p__Pseudomonadota;c__;o__;f__;g__Escherichia;s__\n"
        "sample10.fna\td__Bacteria;p__Bacillota;c__;o__;f__;g__Bacillus;s__\n"
    )
    out_tax = tmp_path / "taxonomy.tsv"
    assert gtdbtk_to_bigslice_taxonomy.main([str(bac_tsv), "-o", str(out_tax)]) == 0
    rows = [line.split("\t") for line in out_tax.read_text().splitlines()[1:]]
    prefixes = [r[0] for r in rows]
    assert prefixes == ["sample1/", "sample10/"]
    assert not prefixes[1].startswith(prefixes[0])


def test_extract_bigfam_scores_from_sqlite_report(tmp_path):
    report_dir = tmp_path / "reports" / "1"
    report_dir.mkdir(parents=True)
    db_path = report_dir / "data.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE bgc (id INTEGER PRIMARY KEY, name TEXT, orig_filename TEXT)")
    conn.execute("CREATE TABLE gcf_membership (gcf_id INTEGER, bgc_id INTEGER, membership_value REAL, rank INTEGER)")
    conn.execute("INSERT INTO bgc VALUES (1, '/work/query_gbks/sample1.region001', 'sample1.region001.gbk')")
    conn.execute("INSERT INTO bgc VALUES (2, 'bgc_novel', 'sample2.region001.gbk')")
    conn.execute("INSERT INTO bgc VALUES (3, 'bgc_unassigned', 'sample3.region001.gbk')")
    conn.execute("INSERT INTO gcf_membership VALUES (42, 1, 320.5, 0)")
    conn.execute("INSERT INTO gcf_membership VALUES (99, 2, 1125.0, 0)")
    conn.commit()
    conn.close()

    out_tsv = tmp_path / "bigfam_novelty_scores.tsv"
    rc = extract_bigfam_scores.main([str(tmp_path / "reports"), str(out_tsv), "--novelty-threshold", "900"])
    assert rc == 0
    assert b"\r\n" not in out_tsv.read_bytes()

    with out_tsv.open() as fh:
        rows = {row["bgc_id"]: row for row in csv.DictReader(fh, delimiter="\t")}

    assert rows["sample1.region001"]["membership_type"] == "core"
    assert rows["sample1.region001"]["is_putative_novel"] == "false"
    assert rows["sample1.region001"]["gcf_id"] == "GCF_42"

    assert rows["bgc_novel"]["membership_type"] == "orphan_novel"
    assert rows["bgc_novel"]["is_putative_novel"] == "true"

    assert rows["bgc_unassigned"]["membership_type"] == "unassigned"
    assert rows["bgc_unassigned"]["is_putative_novel"] == "true"
