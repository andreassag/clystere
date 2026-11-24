"""Unit tests for bin/prepare_bgc_reference.py and bin/aggregate_bgc_abundance.py."""

import csv
from pathlib import Path

import aggregate_bgc_abundance as aba
from Bio import SeqIO
import prepare_bgc_reference as pbr

FIXTURES = Path(__file__).parents[1] / "data" / "combgc"


def test_prepare_bgc_reference_deduplicates_sequences(tmp_path):
    in_dir = tmp_path / "regions"
    in_dir.mkdir()
    src = FIXTURES / "sample1.region001.gbk"
    (in_dir / "sample1.region001.gbk").write_text(src.read_text())
    # Exact duplicate sequence under another name should be deduplicated
    (in_dir / "sample1_dup.region001.gbk").write_text(src.read_text())

    out_fa = tmp_path / "bgc_catalog.fna"
    out_meta = tmp_path / "bgc_catalog_metadata.tsv"
    rc = pbr.main(["--gbks", str(in_dir), "--out-fasta", str(out_fa), "--out-metadata", str(out_meta)])
    assert rc == 0

    records = list(SeqIO.parse(out_fa, "fasta"))
    assert len(records) == 1
    with out_meta.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert int(rows[0]["length_bp"]) == len(records[0].seq)


def test_aggregate_bgc_abundance_dna_and_rna_matrices(tmp_path):
    meta_tsv = tmp_path / "bgc_catalog_metadata.tsv"
    meta_tsv.write_text(
        "bgc_id\tsample_id\tsource_file\tlength_bp\tproduct\tgcf_id\n"
        "bgc_1\ts1\ts1.region001.gbk\t10000\tNRPS\tGCF_0001\n"
        "bgc_2\ts2\ts2.region001.gbk\t12000\tNRPS\tGCF_0001\n"
        "bgc_3\ts2\ts2.region002.gbk\t8000\tTerpene\tNA\n"
    )
    manifest = tmp_path / "assay_manifest.tsv"
    manifest.write_text("sample_id\tassay_type\nmg_1\tmetagenome\nmt_1\tmetatranscriptome\n")

    mg_cov = tmp_path / "mg_1.depth.tsv"
    mg_cov.write_text(
        "Contig\tmg_1.sorted TPM\tmg_1.sorted RPKM\tmg_1.sorted Mean\t"
        "mg_1.sorted Covered Fraction\tmg_1.sorted Read Count\n"
        "bgc_1\t100.0\t50.0\t10.0\t0.95\t500\n"
        "bgc_2\t200.0\t80.0\t20.0\t0.85\t900\n"
        "bgc_3\t0.0\t0.0\t0.0\t0.0\t0\n"
    )
    mt_cov = tmp_path / "mt_1.depth.tsv"
    mt_cov.write_text(
        "Contig\tmt_1.sorted TPM\tmt_1.sorted RPKM\tmt_1.sorted Mean\t"
        "mt_1.sorted Covered Fraction\tmt_1.sorted Read Count\n"
        "bgc_1\t300.0\t120.0\t30.0\t0.90\t1200\n"
        "bgc_2\t50.0\t20.0\t5.0\t0.40\t200\n"
        "bgc_3\t0.0\t0.0\t0.0\t0.0\t0\n"
    )

    outdir = tmp_path / "abund_out"
    rc = aba.main(
        [
            "--coverm",
            str(mg_cov),
            str(mt_cov),
            "--metadata",
            str(meta_tsv),
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
        ]
    )
    assert rc == 0

    with (outdir / "gcf_abundance_tpm.tsv").open() as fh:
        gcf_rows = {r["gcf_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    # GCF_0001 sums bgc_1 + bgc_2
    assert float(gcf_rows["GCF_0001"]["mg_1"]) == 300.0
    assert float(gcf_rows["GCF_0001"]["mt_1"]) == 350.0
    assert "singleton_bgc_3" in gcf_rows

    with (outdir / "bgc_activity_summary.tsv").open() as fh:
        act = {r["bgc_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert float(act["bgc_1"]["rna_dna_tpm_ratio"]) == 3.0
    assert act["bgc_1"]["prevalence_samples"] == "2"
    assert act["bgc_2"]["prevalence_samples"] == "1"  # mt_1 cov=0.40 < 0.50
    assert act["bgc_3"]["rna_dna_tpm_ratio"] == "NA"

    with (outdir / "bgc_abundance_mqc.tsv").open() as fh:
        mqc = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert mqc["mg_1"]["assay_type"] == "metagenome"
    assert mqc["mg_1"]["detected_bgcs"] == "2"
    assert mqc["mt_1"]["assay_type"] == "metatranscriptome"
    assert mqc["mt_1"]["detected_bgcs"] == "1"
