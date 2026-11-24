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


def test_sanntis_workflow_parses_gff_and_handles_empty(tmp_path):
    gff = tmp_path / "sample1.sanntis.full.gff"
    gff.write_text(
        "##gff-version 3\n"
        "ctg1\tSanntiS\tcluster\t100\t1200\t0.91\t+\t.\t"
        "ID=SanntiS_0001;nearest_MiBIG=BGC0000349;nearest_MiBIG_class=NRP;"
        "nearest_MiBIG_diceScore=0.87;partial=false;cds_count=4\n"
        "ctg1\tSanntiS\tcluster\t2000\t2500\t0.65\t-\t.\t"
        "ID=SanntiS_0002;nearest_MiBIG_class=Polyketide;partial=true\n"
    )
    df = comBGC.sanntis_workflow(gff)
    assert len(df) == 2
    assert list(df["Prediction_tool"]) == ["SanntiS", "SanntiS"]
    assert df.iloc[0]["Sample_ID"] == "sample1"
    assert df.iloc[0]["Product_class"] == "NRP"
    assert df.iloc[0]["MIBiG_ID"] == "BGC0000349"
    assert df.iloc[0]["BGC_probability"] == "0.87"
    assert df.iloc[0]["BGC_complete"] == "Yes"
    assert df.iloc[0]["CDS_count"] == 4
    assert df.iloc[1]["BGC_complete"] == "No"
    assert df.iloc[1]["MIBiG_ID"] == "NA"

    empty_gff = tmp_path / "empty.sanntis.full.gff"
    empty_gff.write_text("##gff-version 3\n")
    empty_df = comBGC.sanntis_workflow(empty_gff)
    assert empty_df.empty
    assert list(empty_df.columns) == comBGC.INTERNAL_COLS


def test_sanntis_workflow_parses_real_v0935_attributes(tmp_path):
    """Real SanntiS 0.9.3.5 writes CLUSTER with col 6 '.', score=..., nearest_MiBIG_diceDistance=..., partial=00/10."""
    gff = tmp_path / "bsub.sanntis.full.gff"
    gff.write_text(
        "##gff-version 3\n"
        "NC_000964.3\tSanntiSv0.9.3.5\tCLUSTER\t1000\t15000\t.\t.\t.\t"
        "ID=NC_000964.3_sanntis_1;nearest_MiBIG=BGC0000433;nearest_MiBIG_class=NRP;"
        "nearest_MiBIG_diceDistance=0.125;score=0.941;partial=00\n"
        "NC_000964.3\tSanntiSv0.9.3.5\tCLUSTER\t50000\t60000\t.\t.\t.\t"
        "ID=NC_000964.3_sanntis_2;nearest_MiBIG=BGC0000309;nearest_MiBIG_class=Polyketide;"
        "nearest_MiBIG_diceDistance=0.250;partial=10\n"
        "NC_000964.3\tSanntiSv0.9.3.5\tCLUSTER_border\t1200\t14800\t.\t.\t.\t"
        "ID=NC_000964.3_sanntis_1_1;nearest_MiBIG=BGC0000433;partial=00\n"
    )
    df = comBGC.sanntis_workflow(gff)
    assert len(df) == 2
    assert df.iloc[0]["BGC_probability"] == "0.941"
    assert df.iloc[0]["BGC_complete"] == "Yes"
    assert df.iloc[0]["MIBiG_ID"] == "BGC0000433"
    assert df.iloc[1]["BGC_probability"] == "0.750"
    assert df.iloc[1]["BGC_complete"] == "No"


def test_main_combines_all_four_predictors_including_sanntis(tmp_path):
    sanntis_gff = tmp_path / "sample1.sanntis.full.gff"
    sanntis_gff.write_text(
        "##gff-version 3\n"
        "ctg1.1\tSanntiS\tcluster\t10\t200\t0.88\t+\t.\t"
        "ID=SanntiS_1;nearest_MiBIG=BGC0000349;nearest_MiBIG_class=Siderophore;"
        "nearest_MiBIG_diceScore=0.79;partial=false\n"
    )
    classified = comBGC.classify_all_inputs([str(FIXTURES / "sample1.gbk"), str(sanntis_gff)])
    assert classified["SanntiS"] == [str(sanntis_gff)]

    outdir = tmp_path / "out_4tools"
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
            str(sanntis_gff),
        ]
    )
    assert rc == 0
    summary = pd.read_csv(outdir / "combgc_summary.tsv", sep="\t")
    assert set(summary["Prediction_tool"]) == {"antiSMASH", "GECCO", "deepBGC", "SanntiS"}
