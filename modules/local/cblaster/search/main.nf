/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CBLASTER_SEARCH
    Searches query BGC regions against a local cohort/custom DIAMOND+SQLite database
    or remotely against NCBI via `cblaster search`.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CBLASTER_SEARCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cblaster:1.3.20--pyhdfd78af_0'
        : 'quay.io/biocontainers/cblaster:1.3.20--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(query_gbk)
    path db_dmnd
    path db_sqlite
    val mode

    output:
    tuple val(meta), path("${prefix}_cblaster.tsv"), emit: summary_tsv
    tuple val(meta), path("${prefix}_binary.tsv"), optional: true, emit: binary_tsv
    tuple val("${task.process}"), val('cblaster'), eval("cblaster --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | tail -n 1"), topic: versions, emit: versions_cblaster_search

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def db_arg = mode == 'local' ? "-db ${db_dmnd}" : ''
    def cpu_arg = mode == 'local' ? "-cp ${task.cpus}" : ''
    """
    mkdir -p .config
    cblaster config --email clystere@local.invalid
    cblaster search \\
        -m ${mode} \\
        ${db_arg} \\
        -qf "${query_gbk}" \\
        ${cpu_arg} \\
        -o "${prefix}_cblaster.tsv" \\
        -ode \$'\\t' \\
        -b "${prefix}_binary.tsv" \\
        -bde \$'\\t' \\
        ${args}

    if [[ ! -f "${prefix}_cblaster.tsv" ]]; then
        printf 'Organism\\tScaffold\\tStart\\tEnd\\tScore\\n' > "${prefix}_cblaster.tsv"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Organism\\tScaffold\\tStart\\tEnd\\tScore\\n' > "${prefix}_cblaster.tsv"
    printf 'Organism\\tScaffold\\tStart\\tEnd\\tScore\\n' > "${prefix}_binary.tsv"
    """
}
