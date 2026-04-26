/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GECCO process
    Runs GECCO v0.10.3 on a single genome file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GECCO_GECCO {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/gecco:0.10.3--pyhdfd78af_0'

    input:
    tuple val(meta), path(genome), path(annotation), path(gecco_compat)

    output:
    tuple val(meta), path("${meta.id}/"), emit: output_dir
    tuple val(meta), path("${meta.id}/*.clusters.tsv"), emit: clusters_tsv
    tuple val(meta), path("${meta.id}/gecco_bigslice/"), emit: bigslice_dir, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    def convert_cmd = params.bigscape_run || params.bigslice_run
        ? """
    mkdir -p ${prefix}/gecco_bigslice
    gecco convert gbk --input-dir ${prefix} --format bigslice --output-dir ${prefix}/gecco_bigslice
  """
        : ''
    """
    python ${gecco_compat} run \\
        --genome ${genome} \\
        --output-dir ${prefix} \\
        --jobs ${task.cpus} \\
        --force-tsv \\
        ${args}

    ${convert_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gecco: \$(python -c "import importlib.metadata as m; print(m.version('gecco'))")
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p ${prefix}
    printf 'sequence_id\tbgc_id\ttype\taverage_p\tstart\tend\tdomains\tproteins\n' > ${prefix}/${prefix}.clusters.tsv
    mkdir -p ${prefix}/gecco_bigslice
    touch ${prefix}/gecco_bigslice/${prefix}.region001.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gecco: 0.10.3
    END_VERSIONS
    """
}
