/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIGSCAPE
    Runs BiG-SCAPE v2 on all unified BGC region directories.
    Optionally runs `bigscape dereplicate` before clustering to remove near-identical
    representative sequences across tools.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BIGSCAPE {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigscape:2.0.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/bigscape:2.0.2--pyhdfd78af_0'}"

    input:
    path bgc_dirs, stageAs: 'bgc_inputs/*'
    path pfam
    val dereplicate
    val zip_output

    output:
    path "bigscape_dereplicate/", emit: dereplicate_dir, optional: true
    path "bigscape/bigscape.db", emit: db
    path "bigscape/", emit: results_dir
    path "bigscape.tar.gz", emit: archive, optional: true
    tuple val("${task.process}"), val('bigscape'), val("2.0.2"), topic: versions, emit: versions_bigscape

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def pfam_is_gz = pfam && pfam.name.endsWith('.gz')
    def pfam_resolved = pfam ? (pfam.isDirectory() ? "${pfam}/Pfam-A.hmm" : (pfam_is_gz ? 'Pfam-A.hmm' : "${pfam}")) : ''
    def pfam_arg = pfam_resolved ? "-p ${pfam_resolved}" : ''
    def unpack_pfam = pfam_is_gz
        ? "python3 -c \"import gzip, shutil; shutil.copyfileobj(gzip.open('${pfam}', 'rb'), open('Pfam-A.hmm', 'wb'))\""
        : ''
    def dereplicate_cmd = dereplicate
        ? """
    bigscape dereplicate \\
        -i bgc_inputs \\
        -o bigscape_dereplicate \\
        -c ${task.cpus} \\
        ${args2} \\
        --exclude-gbk final
    """
        : ''
    def cluster_input = dereplicate ? 'bigscape_dereplicate/representative_clusters' : 'bgc_inputs'
    def archive_cmd = zip_output
        ? "tar --exclude 'bigscape/cache' -zcf bigscape.tar.gz bigscape/"
        : ''
    """
    ${unpack_pfam}
    ${dereplicate_cmd}

    bigscape cluster \\
        -i ${cluster_input} \\
        -o bigscape \\
        -c ${task.cpus} \\
        ${pfam_arg} \\
        --force-gbk \\
        --exclude-gbk final \\
        ${args}

    ${archive_cmd}

    if [[ "${pfam_is_gz}" == "true" ]]; then
        rm -f Pfam-A.hmm
    fi
    """

    stub:
    def derep_stub = dereplicate ? 'mkdir -p bigscape_dereplicate/representative_clusters' : ''
    def zip_stub = zip_output ? 'touch bigscape.tar.gz' : ''
    """
    ${derep_stub}
    mkdir -p bigscape
    touch bigscape/bigscape.db
    ${zip_stub}
    """
}
