"""Unit tests for bin/comBGC.py."""

from pathlib import Path
import shutil

import comBGC
import pandas as pd

FIXTURES = Path(__file__).parents[1] / "data" / "combgc"


def test_parse_knownclusterblast():
    ids = comBGC.parse_knownclusterblast(FIXTURES / "knownclusterblast" / "ctg1.1_c1.txt")
    assert ids == ["BGC0000349"]


def test_antismash_workflow_no_duplicate_and_order_independent():
    gbk = str(FIXTURES / "sample1.gbk")
    kcb = str(FIXTURES / "knownclusterblast")

    df_a = comBGC.antismash_workflow([gbk, kcb])
    df_b = comBGC.antismash_workflow([kcb, gbk])

    assert len(df_a) == 1
    assert df_a.iloc[0]["MIBiG_ID"] == "BGC0000349"
    assert df_a.iloc[0]["Product_class"] == "Siderophore"
    assert df_a.iloc[0]["PFAM_domains"] == "IucA_IucC"
    assert df_a.equals(df_b)


def test_antismash_workflow_matches_non_c_contig_kcb(tmp_path):
    gbk_text = (FIXTURES / "sample1.gbk").read_text().replace("ctg1", "NC_000913")
    gbk_path = tmp_path / "NC_sample.gbk"
    gbk_path.write_text(gbk_text)

    kcb_dir = tmp_path / "knownclusterblast"
    kcb_dir.mkdir()
    shutil.copy(FIXTURES / "knownclusterblast" / "ctg1.1_c1.txt", kcb_dir / "NC_000913.1_c1.txt")

    df = comBGC.antismash_workflow([str(gbk_path), str(kcb_dir)])
    assert len(df) == 1
    assert df.iloc[0]["Contig_ID"] == "NC_000913.1"
    assert df.iloc[0]["MIBiG_ID"] == "BGC0000349"


def test_deepbgc_workflow_parses_and_handles_empty(tmp_path):
    df = comBGC.deepbgc_workflow(FIXTURES / "sample1.bgc.tsv")
    assert len(df) == 1
    assert df.iloc[0]["Prediction_tool"] == "deepBGC"
    assert df.iloc[0]["Sample_ID"] == "sample1"

    empty_tsv = tmp_path / "empty.bgc.tsv"
    header = (FIXTURES / "sample1.bgc.tsv").read_text().splitlines()[0] + "\n"
    empty_tsv.write_text(header)
    empty_df = comBGC.deepbgc_workflow(empty_tsv)
    assert empty_df.empty
    assert list(empty_df.columns) == comBGC.INTERNAL_COLS


def test_gecco_workflow_parses_interpro_and_handles_empty(tmp_path):
    df = comBGC.gecco_workflow([str(FIXTURES / "sample1.clusters.tsv"), str(FIXTURES / "sample1_cluster_1.gbk")])
    assert len(df) == 1
    assert df.iloc[0]["Prediction_tool"] == "GECCO"
    assert df.iloc[0]["InterPro_ID"] == "IPR007310"

    empty_tsv = tmp_path / "empty.clusters.tsv"
    header = (FIXTURES / "sample1.clusters.tsv").read_text().splitlines()[0] + "\n"
    empty_tsv.write_text(header)
    empty_df = comBGC.gecco_workflow([str(empty_tsv)])
    assert empty_df.empty
    assert list(empty_df.columns) == comBGC.INTERNAL_COLS


def test_gecco_converted_region_gbk_extracts_interpro(tmp_path):
    converted_gbk = tmp_path / "sample1.region001.gbk"
    shutil.copy(FIXTURES / "sample1_cluster_1.gbk", converted_gbk)

    as_in, _, ge_in = comBGC.classify_inputs([str(FIXTURES / "sample1.clusters.tsv"), str(converted_gbk)])
    assert as_in == []
    assert len(ge_in) == 2

    df = comBGC.gecco_workflow(ge_in)
    assert len(df) == 1
    assert df.iloc[0]["InterPro_ID"] == "IPR007310"


def test_classify_inputs_filters_regions_and_recognizes_all_tools():
    as_in, db_in, ge_in = comBGC.classify_inputs(
        [
            str(FIXTURES / "sample1.gbk"),
            str(FIXTURES / "sample1.region001.gbk"),
            str(FIXTURES / "knownclusterblast"),
            str(FIXTURES / "sample1.clusters.tsv"),
            str(FIXTURES / "sample1_cluster_1.gbk"),
            str(FIXTURES / "sample1.bgc.tsv"),
        ]
    )
    assert len(as_in) == 2
    assert db_in.endswith("sample1.bgc.tsv")
    assert len(ge_in) == 2


def test_main_writes_combined_summary_tsv(tmp_path):
    outdir = tmp_path / "out"
    rc = comBGC.main(
        [
            "--outdir",
            str(outdir),
            "--input",
            str(FIXTURES / "sample1.gbk"),
            str(FIXTURES / "knownclusterblast"),
            str(FIXTURES / "sample1.clusters.tsv"),
            str(FIXTURES / "sample1_cluster_1.gbk"),
            str(FIXTURES / "sample1.bgc.tsv"),
        ]
    )
    assert rc == 0
    summary = pd.read_csv(outdir / "combgc_summary.tsv", sep="\t")
    assert list(summary.columns) == comBGC.OUTPUT_COLS
    assert set(summary["Prediction_tool"]) == {"antiSMASH", "GECCO", "deepBGC"}


def test_main_empty_inputs_writes_header_only_tsv(tmp_path):
    work = tmp_path / "empty_in"
    work.mkdir()
    empty_tsv = work / "empty.bgc.tsv"
    header = (FIXTURES / "sample1.bgc.tsv").read_text().splitlines()[0] + "\n"
    empty_tsv.write_text(header)

    outdir = tmp_path / "out_empty"
    rc = comBGC.main(["--outdir", str(outdir), "--input", str(empty_tsv)])
    assert rc == 0
    lines = (outdir / "combgc_summary.tsv").read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t") == comBGC.OUTPUT_COLS


def test_prepare_multisample_input_antismash(tmp_path):
    sdir = tmp_path / "sampleA"
    (sdir / "knownclusterblast").mkdir(parents=True)
    (sdir / "index.html").write_text("<html></html>")
    shutil.copy(FIXTURES / "sample1.gbk", sdir / "sampleA.gbk")

    paths = comBGC.prepare_multisample_input_antismash(tmp_path)
    assert len(paths) == 1
    assert paths[0][0].endswith("sampleA.gbk")
    assert paths[0][1].endswith("knownclusterblast")
