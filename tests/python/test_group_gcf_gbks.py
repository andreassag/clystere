"""Unit tests for bin/group_gcf_gbks.py."""

import csv
import sqlite3

import group_gcf_gbks


def _write_dummy_gbk(path, name="stub"):
    path.write_text(
        f"LOCUS       {name:<16} 100 bp    DNA     linear   BCT 01-JAN-2026\nORIGIN\n        1 atgcgtacgt\n//\n"
    )


def test_group_gcfs_from_bigscape_db_and_bigslice_tsv(tmp_path):
    gbk_dir = tmp_path / "gbks"
    gbk_dir.mkdir()
    g1 = gbk_dir / "sample1.region001.gbk"
    g2 = gbk_dir / "sample2.region001.gbk"
    g3 = gbk_dir / "sample3.region001.gbk"
    for g in (g1, g2, g3):
        _write_dummy_gbk(g, g.stem)

    # Create mock BiG-SCAPE v2 SQLite DB with 1 multi-member GCF (g1, g2) and 1 singleton (g3)
    db_path = tmp_path / "bigscape.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE gbk (id INTEGER PRIMARY KEY, path TEXT)")
    conn.execute("CREATE TABLE bgc_record (id INTEGER PRIMARY KEY, gbk_id INTEGER)")
    conn.execute("CREATE TABLE family (id INTEGER PRIMARY KEY, center_id INTEGER, cutoff REAL)")
    conn.execute("CREATE TABLE bgc_record_family (record_id INTEGER, family_id INTEGER)")
    conn.executemany(
        "INSERT INTO gbk VALUES (?, ?)",
        [
            (1, "bgc_inputs/sample1.region001.gbk"),
            (2, "bgc_inputs/sample2.region001.gbk"),
            (3, "bgc_inputs/sample3.region001.gbk"),
        ],
    )
    conn.executemany("INSERT INTO bgc_record VALUES (?, ?)", [(10, 1), (20, 2), (30, 3)])
    conn.executemany("INSERT INTO family VALUES (?, ?, ?)", [(1, 10, 0.3), (2, 30, 0.3), (3, 10, 0.7)])
    conn.executemany(
        "INSERT INTO bgc_record_family VALUES (?, ?)",
        [(10, 1), (20, 1), (30, 2), (10, 3), (20, 3), (30, 3)],
    )
    conn.commit()
    conn.close()

    # Also create mock BiG-SLiCE TSV export grouping g1 and g2
    sl_dir = tmp_path / "tsv_export"
    sl_dir.mkdir()
    (sl_dir / "bgc_metadata.tsv").write_text(
        "bgc_id\torig_filename\n1\tsample1.region001.gbk\n2\tsample2.region001.gbk\n"
    )
    (sl_dir / "gcf_membership.tsv").write_text("gcf_id\tbgc_id\tmembership_value\trank\n7\t1\t0.1\t0\n7\t2\t0.25\t0\n")

    outdir = tmp_path / "out"
    rc = group_gcf_gbks.main(
        [
            "--gbks",
            str(gbk_dir),
            "--bigscape-db",
            str(db_path),
            "--bigslice-tsv",
            str(sl_dir),
            "--outdir",
            str(outdir),
        ]
    )
    assert rc == 0

    bs_gcf = outdir / "gcf_groups" / "bigscape_GCF_0001"
    sl_gcf = outdir / "gcf_groups" / "bigslice_GCF_0007"
    assert bs_gcf.is_dir()
    assert sorted(p.name for p in bs_gcf.glob("*.gbk")) == ["sample1.region001.gbk", "sample2.region001.gbk"]
    assert sl_gcf.is_dir()
    assert sorted(p.name for p in sl_gcf.glob("*.gbk")) == ["sample1.region001.gbk", "sample2.region001.gbk"]
    # Singleton GCF_0002 should NOT get a clinker folder (< 2 members)
    assert not (outdir / "gcf_groups" / "bigscape_GCF_0002").exists()

    # Representatives: g1 (center of GCF 1) and g3 (center of singleton GCF 2)
    queries = sorted(p.name for p in (outdir / "cblaster_queries").glob("*.gbk"))
    assert queries == ["sample1.region001.gbk", "sample3.region001.gbk"]

    with (outdir / "gcf_membership_summary.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 5
    assert {r["tool"] for r in rows} == {"BiG-SCAPE", "BiG-SLiCE"}


def test_group_gcfs_unclustered_falls_back_to_all_queries(tmp_path):
    gbk_dir = tmp_path / "gbks"
    gbk_dir.mkdir()
    _write_dummy_gbk(gbk_dir / "s1.region001.gbk")
    outdir = tmp_path / "out"
    assert group_gcf_gbks.main(["--gbks", str(gbk_dir), "--outdir", str(outdir)]) == 0
    assert list((outdir / "gcf_groups").iterdir()) == []
    assert [p.name for p in (outdir / "cblaster_queries").glob("*.gbk")] == ["s1.region001.gbk"]


def test_group_gcfs_parses_real_bigslice_v2_export_tsv_wide_format(tmp_path):
    """Real BiG-SLiCE 2.0.2 `--export-tsv` writes unnamed index + `bgc` col and `bgc_id\trun#1` pivot table."""
    gbk_dir = tmp_path / "gbks"
    gbk_dir.mkdir()
    _write_dummy_gbk(gbk_dir / "sample1.region001.gbk", "s1")
    _write_dummy_gbk(gbk_dir / "sample2.region001.gbk", "s2")
    _write_dummy_gbk(gbk_dir / "sample3.region001.gbk", "s3")

    sl_dir = tmp_path / "tsv_export"
    sl_dir.mkdir()
    (sl_dir / "bgc_metadata.tsv").write_text(
        "\tdataset\tgenome\tbgc\tcontig_edge\tlen_nt\n"
        "1\tdataset\tsample1\tsample1.region001.gbk\tFalse\t6000\n"
        "2\tdataset\tsample2\tsample2.region001.gbk\tFalse\t6000\n"
        "3\tdataset\tsample3\tsample3.region001.gbk\tTrue\t5000\n"
    )
    (sl_dir / "gcf_membership.tsv").write_text("bgc_id\trun#1\trun#2\n1\t1\t5\n2\t1\t5\n3\t2\t-1\n")

    outdir = tmp_path / "out_v2"
    rc = group_gcf_gbks.main(["--gbks", str(gbk_dir), "--bigslice-tsv", str(sl_dir), "--outdir", str(outdir)])
    assert rc == 0

    gcf_dir = outdir / "gcf_groups" / "bigslice_GCF_0005"
    assert gcf_dir.is_dir()
    assert sorted(p.name for p in gcf_dir.glob("*.gbk")) == ["sample1.region001.gbk", "sample2.region001.gbk"]
    with (outdir / "gcf_membership_summary.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
    assert {r["bgc_file"] for r in rows} == {"sample1.region001.gbk", "sample2.region001.gbk"}
    assert all(r["gcf_id"] == "bigslice_GCF_0005" for r in rows)
