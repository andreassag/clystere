/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GTDBTK_DOWNLOADDB
    Optional helper process to download and unpack the GTDB-Tk reference package
    when --gtdbtk_download_db is explicitly enabled.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GTDBTK_DOWNLOADDB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data'
        : 'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd'}"

    input:
    val db_url

    output:
    tuple val('gtdbtk_db'), path('gtdbtk_db'), emit: database
    tuple val("${task.process}"), val('aria2'), eval("aria2c --version 2>&1 | sed '1s/[^ ]* [^ ]* //; q'"), topic: versions, emit: versions_gtdbtk_downloaddb

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p gtdbtk_db
    aria2c ${args} -x 4 -s 4 -o gtdbtk_data.tar.gz "${db_url}"
    tar -xzf gtdbtk_data.tar.gz -C gtdbtk_db --strip-components=1
    rm -f gtdbtk_data.tar.gz
    """

    stub:
    """
    mkdir -p gtdbtk_db/metadata
    printf 'VERSION_DATA=r220\\n' > gtdbtk_db/metadata/metadata.txt
    """
}
