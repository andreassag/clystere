/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CBLASTER_MAKEDB
    Builds a local DIAMOND + SQLite cohort database from sample GenBank genomes for
    offline cblaster gene-neighborhood searches.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CBLASTER_MAKEDB {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cblaster:1.3.20--pyhdfd78af_0'
        : 'quay.io/biocontainers/cblaster:1.3.20--pyhdfd78af_0'}"

    input:
    path gbks, stageAs: 'cohort_gbks/?/*'

    output:
    path 'cohort.dmnd', emit: dmnd
    path 'cohort.sqlite3', emit: sqlite
    tuple val("${task.process}"), val('cblaster'), eval("cblaster --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | tail -n 1"), topic: versions, emit: versions_cblaster_makedb

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p .config
    cblaster config --email clystere@local.invalid
    find -L cohort_gbks -type f \\( -name '*.gbk' -o -name '*.gbff' -o -name '*.gb' \\) > gbk_list.txt
    cblaster makedb \\
        gbk_list.txt \\
        -n cohort \\
        -cp ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch cohort.dmnd
    touch cohort.sqlite3
    """
}
