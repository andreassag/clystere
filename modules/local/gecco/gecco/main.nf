process GECCO_GECCO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/gecco:0.10.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}/"), emit: output_dir
    tuple val(meta), path("${meta.id}/*.genes.tsv")    , optional: true, emit: genes
    tuple val(meta), path("${meta.id}/*.features.tsv")                 , emit: features
    tuple val(meta), path("${meta.id}/*.clusters.tsv") , optional: true, emit: clusters
    tuple val(meta), path("${meta.id}/*_cluster_*.gbk"), optional: true, emit: gbk
    tuple val(meta), path("${meta.id}/*.json")         , optional: true, emit: json
    tuple val("${task.process}"), val('gecco'), eval('gecco -V | cut -d" " -f2'), emit: versions_gecco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    gecco \\
        run \\
        $args \\
        -j $task.cpus \\
        -o ${prefix} \\
        -g ${input}

    for i in \$(find ${prefix} -name '${input.baseName}*' -type f); do
        mv \$i \${i/${input.baseName}/${prefix}};
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.genes.tsv
    touch ${prefix}/${prefix}.features.tsv
    touch ${prefix}/${prefix}.clusters.tsv
    touch ${prefix}/${prefix}_cluster_1.gbk
    """
}
