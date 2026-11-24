/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSENSUS_SUMMARY
    Merges cross-tool BGC calls across up to 4 predictors (antiSMASH, GECCO, deepBGC,
    SanntiS) and joins GCF/BiG-FAM novelty scores, QUAST, CheckM2, GTDB-Tk, and AMPs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CONSENSUS_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    path combgc_tsvs, stageAs: 'combgc_inputs/?/*'
    path all_regions_tsv
    path gcf_membership_tsv
    path bigfam_scores_tsv
    path quast_tsvs, stageAs: 'quast_inputs/?/*'
    path checkm2_tsvs, stageAs: 'checkm2_inputs/?/*'
    path gtdbtk_tsvs, stageAs: 'gtdbtk_inputs/?/*'
    path amp_summary_tsv
    val min_overlap

    output:
    path 'consensus_loci.tsv', emit: consensus_loci_tsv
    path 'genome_qc_bgc_summary.tsv', emit: genome_summary_tsv
    path 'consensus_support_mqc.tsv', emit: consensus_support_mqc
    path 'genome_overview_mqc.tsv', emit: genome_overview_mqc
    path 'amp_summary_mqc.tsv', optional: true, emit: amp_summary_mqc
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_consensus_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def all_regions_arg = all_regions_tsv ? "--all-regions ${all_regions_tsv}" : ''
    def gcf_arg = gcf_membership_tsv ? "--gcf-membership ${gcf_membership_tsv}" : ''
    def bigfam_arg = bigfam_scores_tsv ? "--bigfam-scores ${bigfam_scores_tsv}" : ''
    def amp_arg = amp_summary_tsv ? "--amp-summary ${amp_summary_tsv}" : ''
    """
    shopt -s nullglob
    combgc_files=(combgc_inputs/*/*.tsv)
    quast_files=(quast_inputs/*/*.tsv)
    checkm2_files=(checkm2_inputs/*/*.tsv)
    gtdbtk_files=(gtdbtk_inputs/*/*.tsv)

    cmd=(consensus_novelty_summary.py --min-overlap "${min_overlap}" --outdir .)
    if [[ \${#combgc_files[@]} -gt 0 ]]; then
        cmd+=(--combgc "\${combgc_files[@]}")
    fi
    if [[ \${#quast_files[@]} -gt 0 ]]; then
        cmd+=(--quast "\${quast_files[@]}")
    fi
    if [[ \${#checkm2_files[@]} -gt 0 ]]; then
        cmd+=(--checkm2 "\${checkm2_files[@]}")
    fi
    if [[ \${#gtdbtk_files[@]} -gt 0 ]]; then
        cmd+=(--gtdbtk "\${gtdbtk_files[@]}")
    fi

    "\${cmd[@]}" \\
        ${all_regions_arg} \\
        ${gcf_arg} \\
        ${bigfam_arg} \\
        ${amp_arg}
    """

    stub:
    def amp_stub = amp_summary_tsv
        ? "printf 'sample_id\\tamp_candidates\\n' > amp_summary_mqc.tsv"
        : ''
    """
    printf 'locus_id\\tsample_id\\tcontig_id\\tstart\\tend\\tlength\\ttools\\ttool_count\\tsupport_tier\\tmean_jaccard\\tproduct_classes\\tcomplete_status\\tmibig_ids\\tknownclusterblast_hit\\tknownclusterblast_similarity\\tbigscape_gcf\\tbigslice_gcf\\tbigfam_gcf\\tbigfam_distance_d\\tis_putative_novel\\n' > consensus_loci.tsv
    printf 'sample_id\\tgtdb_classification\\tcompleteness\\tcontamination\\tassembly_length_bp\\tnum_contigs\\tn50\\tantismash_bgcs\\tgecco_bgcs\\tdeepbgc_bgcs\\tsanntis_bgcs\\tconsensus_loci\\thigh_support_loci\\tmedium_support_loci\\tsingle_tool_loci\\tputative_novel_bgcs\\tbgcs_per_mbp\\tcomplete_bgc_fraction\\tcontig_edge_rate\\tamp_count\\n' > genome_qc_bgc_summary.tsv
    printf 'sample_id\\tHigh_3_to_4_tools\\tMedium_2_tools\\tSingle_tool\\n' > consensus_support_mqc.tsv
    printf 'sample_id\\tgtdb_classification\\tcompleteness\\tcontamination\\tn50\\tantismash_bgcs\\tgecco_bgcs\\tdeepbgc_bgcs\\tsanntis_bgcs\\tconsensus_loci\\tputative_novel_bgcs\\tbgcs_per_mbp\\tcontig_edge_rate\\n' > genome_overview_mqc.tsv
    ${amp_stub}
    """
}
