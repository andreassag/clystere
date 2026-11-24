"""Unit tests for bin/consensus_novelty_summary.py."""

import csv

import consensus_novelty_summary as cns


def test_consensus_merges_four_predictors_and_joins_qc_taxonomy_novelty(tmp_path):
    combgc_tsv = tmp_path / "sample1_combgc.tsv"
    combgc_tsv.write_text(
        "sample_id\tcontig_id\tPrediction_tool\tProduct_class\tBGC_probability\tBGC_complete\t"
        "BGC_start\tBGC_end\tBGC_length\tBGC_region_contig_ids\tCDS_count\tPFAM_domains\tMIBiG_ID\tInterPro_ID\n"
        "sample1\tctg1\tantiSMASH\tNRPS\tNA\tYes\t1000\t20000\t19001\tcds1\t10\tAMP-binding\tBGC0000433\tNA\n"
        "sample1\tctg1\tGECCO\tNRP\t0.92\tNA\t1500\t19500\t18000\tcds1\t9\tPF00501\tNA\tIPR000873\n"
        "sample1\tctg1\tdeepBGC\tNRP\t0.89\tNA\t1200\t20500\t19300\tcds1\t10\tPF00501\tNA\tNA\n"
        "sample1\tctg1\tSanntiS\tNRP\t0.84\tYes\t1100\t19800\t18701\tSanntiS_1\t9\tNA\tBGC0000433\tNA\n"
        "sample1\tctg1\tantiSMASH\tTerpene\tNA\tNo\t50000\t60000\t10001\tcds2\t4\tNA\tNA\tNA\n"
        "sample1\tctg1\tGECCO\tTerpene\t0.81\tNA\t51000\t59000\t8000\tcds2\t4\tNA\tNA\tNA\n"
        "sample1\tctg2\tdeepBGC\tPolyketide\t0.75\tNA\t500\t5500\t5000\tcds3\t3\tNA\tNA\tNA\n"
    )

    all_regions = tmp_path / "all_regions.tsv"
    all_regions.write_text(
        "file\trecord_id\tregion\tstart\tend\tcontig_edge\tproduct\t"
        "knownclusterblast_hit\tknownclusterblast_accession\tknownclusterblast_similarity\trecord_desc\n"
        "sample1\tctg1\t1\t1000\t20000\tFalse\tNRPS\tsurfactin\tBGC0000433\thigh\ttest genome\n"
        "sample1\tctg1\t2\t50000\t60000\tTrue\tterpene\t\t\t\ttest genome\n"
    )

    quast_tsv = tmp_path / "quast.tsv"
    quast_tsv.write_text("Assembly\tsample1.fna\n# contigs\t2\nTotal length\t2000000\nN50\t1500000\n")

    checkm2_tsv = tmp_path / "checkm2.tsv"
    checkm2_tsv.write_text(
        "Name\tCompleteness\tContamination\tGenome_Size\tTotal_Contigs\tContig_N50\n"
        "sample1\t98.5\t0.8\t2000000\t2\t1500000\n"
    )

    gtdbtk_tsv = tmp_path / "gtdbtk.bac120.summary.tsv"
    gtdbtk_tsv.write_text(
        "user_genome\tclassification\n"
        "sample1.fna\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;"
        "f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\n"
    )

    bigfam_tsv = tmp_path / "bigfam_novelty_scores.tsv"
    bigfam_tsv.write_text(
        "bgc_id\tbgc_filename\tgcf_id\tdistance_d\trank\tmembership_type\tis_putative_novel\n"
        "sample1.region001\tsample1.region001.gbk\tGCF_10\t1100.0\t0\torphan_novel\ttrue\n"
    )

    amp_tsv = tmp_path / "Ampcombi_summary.tsv"
    amp_tsv.write_text("sample\tcontig\tamp_seq\nsample1\tctg1\tMKKLLVVA\nsample1\tctg2\tMKKAAVV\n")

    outdir = tmp_path / "summary"
    rc = cns.main(
        [
            "--combgc",
            str(combgc_tsv),
            "--all-regions",
            str(all_regions),
            "--quast",
            str(quast_tsv),
            "--checkm2",
            str(checkm2_tsv),
            "--gtdbtk",
            str(gtdbtk_tsv),
            "--bigfam-scores",
            str(bigfam_tsv),
            "--amp-summary",
            str(amp_tsv),
            "--outdir",
            str(outdir),
        ]
    )
    assert rc == 0

    with (outdir / "consensus_loci.tsv").open() as fh:
        loci = list(csv.DictReader(fh, delimiter="\t"))
    assert len(loci) == 3

    tiers = {row["contig_id"] + ":" + row["start"]: row for row in loci}
    high_locus = tiers["ctg1:1000"]
    assert high_locus["tool_count"] == "4"
    assert high_locus["support_tier"] == "High (3-4 tools)"
    assert high_locus["knownclusterblast_hit"] == "surfactin"
    assert high_locus["complete_status"] == "Complete"
    assert high_locus["is_putative_novel"] == "true"

    med_locus = tiers["ctg1:50000"]
    assert med_locus["tool_count"] == "2"
    assert med_locus["support_tier"] == "Medium (2 tools)"
    assert med_locus["complete_status"] == "Partial"

    single_locus = tiers["ctg2:500"]
    assert single_locus["tool_count"] == "1"
    assert single_locus["support_tier"] == "Single-tool (1 tool)"

    with (outdir / "genome_qc_bgc_summary.tsv").open() as fh:
        genomes = list(csv.DictReader(fh, delimiter="\t"))
    assert len(genomes) == 1
    g = genomes[0]
    assert g["sample_id"] == "sample1"
    assert g["completeness"] == "98.5"
    assert g["n50"] == "1500000"
    assert g["gtdb_classification"].endswith("Bacillus subtilis")
    assert g["antismash_bgcs"] == "2"
    assert g["gecco_bgcs"] == "2"
    assert g["deepbgc_bgcs"] == "2"
    assert g["sanntis_bgcs"] == "1"
    assert g["consensus_loci"] == "3"
    assert g["high_support_loci"] == "1"
    assert g["medium_support_loci"] == "1"
    assert g["single_tool_loci"] == "1"
    assert g["bgcs_per_mbp"] == "1.5"
    assert g["complete_bgc_fraction"] == "0.5"
    assert g["contig_edge_rate"] == "0.5"
    assert g["amp_count"] == "2"

    assert (outdir / "consensus_support_mqc.tsv").is_file()
    assert (outdir / "genome_overview_mqc.tsv").is_file()
    assert (outdir / "amp_summary_mqc.tsv").is_file()


def test_consensus_empty_inputs_emit_headers(tmp_path):
    outdir = tmp_path / "empty_out"
    assert cns.main(["--outdir", str(outdir)]) == 0
    assert (outdir / "consensus_loci.tsv").read_text().splitlines()[0].split("\t") == cns.CONSENSUS_COLS
    assert (outdir / "genome_qc_bgc_summary.tsv").read_text().splitlines()[0].split("\t") == cns.GENOME_SUMMARY_COLS


def test_consensus_joins_gcf_and_bigfam_at_locus_level_in_multi_bgc_genome(tmp_path):
    """Multi-BGC genomes must match each locus to its own exact GCF and BiG-FAM distance_d, not sample aggregate."""
    combgc_tsv = tmp_path / "sample1_combgc.tsv"
    combgc_tsv.write_text(
        "sample_id\tcontig_id\tPrediction_tool\tProduct_class\tBGC_probability\tBGC_complete\t"
        "BGC_start\tBGC_end\tBGC_length\tBGC_region_contig_ids\tCDS_count\tPFAM_domains\tMIBiG_ID\tInterPro_ID\n"
        "sample1\tctg1\tantiSMASH\tNRPS\tNA\tYes\t1000\t15000\t14001\tcds1\t8\tAMP-binding\tNA\tNA\n"
        "sample1\tctg1\tantiSMASH\tTerpene\tNA\tYes\t50000\t65000\t15001\tcds2\t6\tTerpene\tNA\tNA\n"
    )
    all_regions = tmp_path / "all_regions.tsv"
    all_regions.write_text(
        "file\trecord_id\tregion\tstart\tend\tcontig_edge\tproduct\trecord_desc\n"
        "sample1\tctg1\t1\t1000\t15000\tFalse\tNRPS\ttest\n"
        "sample1\tctg1\t2\t50000\t65000\tFalse\tterpene\ttest\n"
    )
    gcf_tsv = tmp_path / "gcf_membership_summary.tsv"
    gcf_tsv.write_text(
        "bgc_file\tsample_id\ttool\tgcf_id\tgcf_size\tis_representative\tcutoff\n"
        "sample1.region001.gbk\tsample1\tBiG-SCAPE\tbigscape_GCF_0001\t3\ttrue\t0.3\n"
        "sample1.region002.gbk\tsample1\tBiG-SCAPE\tbigscape_GCF_0002\t2\tfalse\t0.3\n"
    )
    bigfam_tsv = tmp_path / "bigfam_novelty_scores.tsv"
    bigfam_tsv.write_text(
        "bgc_id\tbgc_filename\tgcf_id\tdistance_d\trank\tmembership_type\tis_putative_novel\n"
        "sample1.region001\tsample1.region001.gbk\tGCF_10\t120.0\t0\tcore\tfalse\n"
        "sample1.region002\tsample1.region002.gbk\tGCF_99\t1150.0\t0\torphan_novel\ttrue\n"
    )

    outdir = tmp_path / "multi_out"
    rc = cns.main(
        [
            "--combgc",
            str(combgc_tsv),
            "--all-regions",
            str(all_regions),
            "--gcf-membership",
            str(gcf_tsv),
            "--bigfam-scores",
            str(bigfam_tsv),
            "--outdir",
            str(outdir),
        ]
    )
    assert rc == 0
    with (outdir / "consensus_loci.tsv").open() as fh:
        by_start = {row["start"]: row for row in csv.DictReader(fh, delimiter="\t")}

    assert by_start["1000"]["bigscape_gcf"] == "bigscape_GCF_0001"
    assert by_start["1000"]["bigfam_gcf"] == "GCF_10"
    assert by_start["1000"]["bigfam_distance_d"] == "120.00"
    assert by_start["1000"]["is_putative_novel"] == "false"

    assert by_start["50000"]["bigscape_gcf"] == "bigscape_GCF_0002"
    assert by_start["50000"]["bigfam_gcf"] == "GCF_99"
    assert by_start["50000"]["bigfam_distance_d"] == "1150.00"
    assert by_start["50000"]["is_putative_novel"] == "true"
