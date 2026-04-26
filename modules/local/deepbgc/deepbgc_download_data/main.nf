/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEEPBGC_DOWNLOAD_DATA
    Downloads deepBGC models and Pfam data once per pipeline run.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEEPBGC_DOWNLOAD {
    label 'process_single'
    conda "${moduleDir}/../deepbgc/environment.yml"
    container 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'

    input:
    val(data_dest)

    output:
    val(data_dest), emit: data_dir
    path "versions.yml", emit: versions, topic: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(deepbgc info 2>&1 | sed -n '6p' | sed 's/.*= version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${data_dest}/detector
    touch ${data_dest}/detector/deepbgc.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: 0.1.31
    END_VERSIONS
    """
}
