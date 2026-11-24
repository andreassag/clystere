/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINKER
    Generates interactive HTML gene cluster alignment & synteny plots for multi-member GCFs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINKER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/clinker-py:0.0.31--pyhdfd78af_0'
        : 'quay.io/biocontainers/clinker-py:0.0.31--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(gcf_dir)

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val(meta), path("${prefix}_alignments.tsv"), emit: alignments
    tuple val("${task.process}"), val('clinker'), eval("clinker --version 2>&1 | sed 's/.*clinker v//; s/clinker //'"), topic: versions, emit: versions_clinker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    clinker \\
        "${gcf_dir}"/*.gbk \\
        -j ${task.cpus} \\
        -p "${prefix}.html" \\
        -o "${prefix}_alignments.tsv" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '<html><body>stub clinker plot for ${prefix}</body></html>\\n' > "${prefix}.html"
    printf 'Query\\tTarget\\tIdentity\\tSimilarity\\n' > "${prefix}_alignments.tsv"
    """
}
