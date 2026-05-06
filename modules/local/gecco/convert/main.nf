process GECCO_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/gecco:0.10.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(clusters), path(gbk)

    output:
    tuple val(meta), path("${prefix}/*.region*.gbk"), emit: bigslice, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gecco \\
        convert \\
        $args \\
        --jobs $task.cpus \\
        gbk \\
        --input-dir ./ \\
        --format bigslice \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gecco: \$(echo \$(gecco --version) | cut -f 2 -d ' ' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.region001.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gecco: \$(echo \$(gecco --version) | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}