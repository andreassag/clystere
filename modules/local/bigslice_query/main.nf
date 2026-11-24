/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIGSLICE_QUERY
    Queries BGC region GBKs against a pre-computed BiG-SLiCE / BiG-FAM database using
    `bigslice --query` and extracts novelty d-scores into bigfam_novelty_scores.tsv.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BIGSLICE_QUERY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0'
        : 'quay.io/biocontainers/bigslice:2.0.2--pyh8ed023e_0'}"

    input:
    tuple val(meta), path(bgc, stageAs: 'query_in/s*/*')
    path hmmdb
    path bigfam_db
    val n_ranks
    val novelty_threshold

    output:
    tuple val(meta), path("${prefix}/bigfam_novelty_scores.tsv"), emit: novelty_tsv
    tuple val(meta), path("${prefix}/reports"), optional: true, emit: reports_dir
    tuple val("${task.process}"), val('bigslice'), val("2.0.2"), topic: versions, emit: versions_bigslice_query

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p query_gbks "${prefix}/result" "${prefix}/reports"
    find -L query_in -name '*.gbk' -size +0c -exec cp {} query_gbks/ \\;

    if [[ -f "${bigfam_db}/result/data.db" ]]; then
        cp "${bigfam_db}/result/data.db" "${prefix}/result/data.db"
        if [[ -d "${bigfam_db}/result/tmp" ]]; then
            ln -s "\$(realpath "${bigfam_db}/result/tmp")" "${prefix}/result/tmp"
        fi
        if [[ -d "${bigfam_db}/result/cache" ]]; then
            ln -s "\$(realpath "${bigfam_db}/result/cache")" "${prefix}/result/cache"
        else
            mkdir -p "${prefix}/result/cache"
        fi
    elif [[ -f "${bigfam_db}/data.db" ]]; then
        cp "${bigfam_db}/data.db" "${prefix}/result/data.db"
        mkdir -p "${prefix}/result/cache"
    else
        echo "ERROR: Could not find result/data.db or data.db inside --bigfam_db (${bigfam_db})" >&2
        exit 1
    fi

    if [[ -n "\$(ls -A query_gbks 2>/dev/null)" ]]; then
        bigslice \\
            --query query_gbks \\
            --query_name clystere_query \\
            --n_ranks ${n_ranks} \\
            --num_threads ${task.cpus} \\
            --program_db_folder "${hmmdb}" \\
            ${args} \\
            "${prefix}"
    fi

    extract_bigfam_scores.py \\
        "${prefix}/reports" \\
        "${prefix}/bigfam_novelty_scores.tsv" \\
        --novelty-threshold "${novelty_threshold}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/reports/1"
    printf 'bgc_id\\tbgc_filename\\tgcf_id\\tdistance_d\\trank\\tmembership_type\\tis_putative_novel\\n' > "${prefix}/bigfam_novelty_scores.tsv"
    printf 'sample1.region001\\tsample1.region001.gbk\\tGCF_1\\t1050.0\\t0\\torphan_novel\\ttrue\\n' >> "${prefix}/bigfam_novelty_scores.tsv"
    """
}
