/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUMMARY
    Produces both per-region and per-genome summary TSVs across all antiSMASH results.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    path antismash_jsons, stageAs: 'antismash_inputs/?/*'

    output:
    path 'all_regions.tsv', emit: all_regions_tsv
    path 'region_counts.tsv', emit: region_counts_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    tabulate_regions.py antismash_inputs all_regions.tsv \\
        ${args} \\
        --threads ${task.cpus}

    count_regions.py antismash_inputs region_counts.tsv \\
        ${args2} \\
        --threads ${task.cpus}
    """

    stub:
    """
    printf 'file\\trecord_id\\tregion\\tstart\\tend\\tcontig_edge\\tproduct\\trecord_desc\\n' > all_regions.tsv
    printf 'record\\ttotal_count\\tdescription\\n' > region_counts.tsv
    """
}
