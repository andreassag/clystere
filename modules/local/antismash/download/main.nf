process ANTISMASH_DOWNLOAD {
    label 'process_single'

    storeDir "${db_dest}"

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/antismash:8.0.1--pyhdfd78af_0'

    input:
    val(db_dest)

    output:
    val(db_dest), emit: databases
    path(".antismash_db_complete"), emit: marker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download-antismash-databases \\
        --database-dir ${db_dest} \\
        ${args}

    chmod -R a+rX ${db_dest}
    touch .antismash_db_complete
    """

    stub:
    """
    mkdir -p ${db_dest}/pfam
    touch .antismash_db_complete
    """
}
