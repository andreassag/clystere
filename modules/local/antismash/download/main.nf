process ANTISMASH_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/antismash:8.0.1--pyhdfd78af_0'

    input:
    val(db_dest)

    output:
    val(db_dest), emit: databases
    tuple val("${task.process}"), val('antismash'), eval("antismash --version | sed 's/antiSMASH //;s/-.*//g'"), emit: versions_antismash, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download-antismash-databases \\
        --database-dir ${db_dest} \\
        ${args}

    chmod -R a+rX ${db_dest}
    """

    stub:
    """
    mkdir -p ${db_dest}

    """
}
