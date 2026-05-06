process DEEPBGC_DOWNLOAD {
    label 'process_single'

    storeDir "${data_dest}"

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'

    input:
    val(data_dest)

    output:
    val(data_dest), emit: data_dir
    path(".deepbgc_data_complete"), emit: marker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${data_dest}
    export DEEPBGC_DOWNLOADS_DIR="${data_dest}"

    deepbgc download \\
        ${args}

    chmod -R a+rX ${data_dest}
    touch .deepbgc_data_complete
    """

    stub:
    """
    mkdir -p ${data_dest}/detector
    touch ${data_dest}/detector/deepbgc.pkl
    touch .deepbgc_data_complete
    """
}
