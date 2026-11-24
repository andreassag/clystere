/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GROUP_GCF_GBKS
    Groups BGC region GBKs into multi-member GCF folders for clinker synteny plotting
    and selects representative/orphan GBKs for cblaster homology searching.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GROUP_GCF_GBKS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    path bgc_dirs, stageAs: 'bgc_inputs/?/*'
    path bigscape_db
    path bigslice_tsv_dir
    val clinker_max_bgcs
    val cblaster_max_queries

    output:
    path 'gcf_groups/*', type: 'dir', optional: true, emit: gcf_dirs
    path 'cblaster_queries/*.gbk', optional: true, emit: cblaster_queries
    path 'gcf_membership_summary.tsv', emit: membership_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_group_gcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def bs_arg = bigscape_db ? "--bigscape-db ${bigscape_db}" : ''
    def sl_arg = bigslice_tsv_dir ? "--bigslice-tsv ${bigslice_tsv_dir}" : ''
    """
    group_gcf_gbks.py \\
        --gbks bgc_inputs \\
        ${bs_arg} \\
        ${sl_arg} \\
        --clinker-max-bgcs ${clinker_max_bgcs} \\
        --cblaster-max-queries ${cblaster_max_queries} \\
        --outdir .
    """

    stub:
    """
    mkdir -p gcf_groups/bigscape_GCF_0001 cblaster_queries
    printf 'LOCUS       stub1                    100 bp    DNA     linear   BCT 01-JAN-2026\\n//\\n' > gcf_groups/bigscape_GCF_0001/sample1.region001.gbk
    printf 'LOCUS       stub2                    100 bp    DNA     linear   BCT 01-JAN-2026\\n//\\n' > gcf_groups/bigscape_GCF_0001/sample2.region001.gbk
    cp gcf_groups/bigscape_GCF_0001/sample1.region001.gbk cblaster_queries/sample1.region001.gbk
    printf 'bgc_file\\tsample_id\\ttool\\tgcf_id\\tgcf_size\\tis_representative\\tcutoff\\n' > gcf_membership_summary.tsv
    printf 'sample1.region001.gbk\\tsample1\\tBiG-SCAPE\\tbigscape_GCF_0001\\t2\\ttrue\\t0.3\\n' >> gcf_membership_summary.tsv
    """
}
