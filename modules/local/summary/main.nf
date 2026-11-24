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
    path 'ripp_candidates.tsv', optional: true, emit: ripp_candidates_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def extract_ripps = task.ext.args3?.contains('--extract-ripps') ?: false
    def ripp_cmd = extract_ripps
        ? "extract_ripp_candidates.py antismash_inputs ripp_candidates.tsv --threads ${task.cpus}"
        : ''
    """
    tabulate_regions.py antismash_inputs all_regions.tsv \\
        ${args} \\
        --threads ${task.cpus}

    count_regions.py antismash_inputs region_counts.tsv \\
        ${args2} \\
        --threads ${task.cpus}

    ${ripp_cmd}
    """

    stub:
    def extract_ripps = task.ext.args3?.contains('--extract-ripps') ?: false
    def ripp_stub = extract_ripps
        ? "printf 'sample_id\\trecord_id\\tregion\\tripp_subclass\\tlocus_tag\\trre_bitscore\\tleader_sequence\\tcore_sequence\\tpredicted_molecular_weight\\ttailoring_modifications\\n' > ripp_candidates.tsv"
        : ''
    """
    printf 'file\\trecord_id\\tregion\\tstart\\tend\\tcontig_edge\\tproduct\\trecord_desc\\n' > all_regions.tsv
    printf 'record\\ttotal_count\\tdescription\\n' > region_counts.tsv
    ${ripp_stub}
    """
}
