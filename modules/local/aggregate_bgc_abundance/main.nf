/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AGGREGATE_BGC_ABUNDANCE
    Aggregates CoverM contig profiles into BGC/GCF TPM, RPKM, covered fraction,
    DNA vs RNA expression activity, and MultiQC summary tables.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process AGGREGATE_BGC_ABUNDANCE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    path coverm_tsvs, stageAs: 'coverm_inputs/*'
    path metadata_tsv
    path gcf_membership_tsv
    path assay_manifest_tsv

    output:
    path 'bgc_abundance_tpm.tsv', emit: bgc_tpm_tsv
    path 'bgc_abundance_rpkm.tsv', emit: bgc_rpkm_tsv
    path 'bgc_covered_fraction.tsv', emit: bgc_covered_fraction_tsv
    path 'gcf_abundance_tpm.tsv', emit: gcf_tpm_tsv
    path 'gcf_abundance_rpkm.tsv', emit: gcf_rpkm_tsv
    path 'bgc_activity_summary.tsv', emit: activity_summary_tsv
    path 'bgc_abundance_mqc.tsv', emit: abundance_mqc
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_aggregate_abundance

    when:
    task.ext.when == null || task.ext.when

    script:
    def meta_arg = metadata_tsv ? "--metadata ${metadata_tsv}" : ''
    def gcf_arg = gcf_membership_tsv ? "--gcf-membership ${gcf_membership_tsv}" : ''
    def manifest_arg = assay_manifest_tsv ? "--manifest ${assay_manifest_tsv}" : ''
    """
    aggregate_bgc_abundance.py \\
        --coverm coverm_inputs/*.tsv \\
        ${meta_arg} \\
        ${gcf_arg} \\
        ${manifest_arg} \\
        --outdir .
    """

    stub:
    """
    printf 'bgc_id\\tsample1\\n' > bgc_abundance_tpm.tsv
    printf 'bgc_id\\tsample1\\n' > bgc_abundance_rpkm.tsv
    printf 'bgc_id\\tsample1\\n' > bgc_covered_fraction.tsv
    printf 'gcf_id\\tsample1\\n' > gcf_abundance_tpm.tsv
    printf 'gcf_id\\tsample1\\n' > gcf_abundance_rpkm.tsv
    printf 'bgc_id\\tgcf_id\\tproduct\\tprevalence_samples\\tmean_dna_tpm\\tmean_rna_tpm\\trna_dna_tpm_ratio\\n' > bgc_activity_summary.tsv
    printf 'sample_id\\tassay_type\\tdetected_bgcs\\ttotal_tpm\\tmean_depth\\n' > bgc_abundance_mqc.tsv
    """
}
