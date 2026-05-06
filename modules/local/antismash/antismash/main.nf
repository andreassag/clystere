process ANTISMASH_ANTISMASH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/antismash:8.0.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(sequence_input)
    path databases

    output:
    tuple val(meta), path("${prefix}/"), emit: output_dir
    tuple val(meta), path("${prefix}/*.json"), emit: json_results
    tuple val(meta), path("${prefix}/index.html"), emit: html, optional: true
    tuple val(meta), path("${prefix}/*.zip"), emit: zip, optional: true
    tuple val(meta), path("${prefix}/*region*.gbk"), emit: gbk_results, optional: true
    tuple val(meta), path("${prefix}/knownclusterblast/"), emit: knownclusterblast_dir, optional: true
    tuple val(meta), path("${prefix}/clusterblast/"), emit: clusterblast_dir, optional: true
    tuple val("${task.process}"), val('antismash'), eval("antismash --version | sed 's/antiSMASH //;s/-.*//g'"), emit: versions_antismash, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reuse_arg = params.antismash_reuse_results ? '--reuse-results' : ''

    """
    antismash \\
        ${args} \\
        -c ${task.cpus} \\
        --output-dir ${prefix} \\
        --output-basename ${prefix} \\
        --genefinding-tool none \\
        --logfile ${prefix}/${prefix}.log \\
        --databases ${databases} \\
        ${reuse_arg} \\
        ${sequence_input}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.json
    touch ${prefix}/index.html
    touch ${prefix}/${prefix}.zip
    touch ${prefix}/${prefix}.region001.gbk
    """
}
