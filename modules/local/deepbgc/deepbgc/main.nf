process DEEPBGC_PIPELINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'

    input:
    tuple val(meta), path(genome)
    path db

    output:
    tuple val(meta), path("${meta.id}/"), emit: output_dir
    tuple val(meta), path("${meta.id}/${meta.id}.bgc.tsv"), optional: true, emit: bgc_tsv
    tuple val(meta), path("${meta.id}/${meta.id}.bgc.gbk"), optional: true, emit: bgc_gbk
    tuple val("${task.process}"), val('deepbgc'), eval("deepbgc info 2>&1 | sed '6!d;s/.*= version //;s/ .*//'") , emit: versions_deepbgc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    export DEEPBGC_DOWNLOADS_DIR=${db}

    deepbgc \\
        pipeline \\
        ${args} \\
        ${genome}

    if [[ "${genome.baseName}/" != "${prefix}/" ]]; then
        mv "${genome.baseName}/" "${prefix}/"
    fi

    for i in \$(find -name '${genome.baseName}*' -type f); do
        mv \$i \${i/${genome.baseName}/${prefix}};
    done
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.bgc.tsv
    touch ${prefix}/${prefix}.bgc.gbk
    """
}
