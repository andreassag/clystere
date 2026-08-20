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
  container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/bigscape:2.0.2--pyhdfd78af_0' :
      'quay.io/biocontainers/bigscape:2.0.2--pyhdfd78af_0' }"

  input:
  path bgc_dirs
  path pfam

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
  def pfam_arg = pfam
    ? (pfam.isDirectory() ? "-p ${pfam}/Pfam-A.hmm" : "-p ${pfam}")
    : ''
  def dereplicate_cmd = params.bigscape_dereplicate
    ? """
    bigscape dereplicate \\
        -i . \\
        -o bigscape_dereplicate \\
        -c ${task.cpus} \\
        --cutoff ${params.bigscape_dereplicate_cutoff} \\
        --exclude-gbk final,deepbgc
  """
    : ''
  def cluster_input = params.bigscape_dereplicate ? 'bigscape_dereplicate/representative_clusters' : '.'
  def archive_cmd = params.bigscape_zip_output
    ? "tar --exclude 'bigscape/cache' -zcf bigscape.tar.gz bigscape/"
    : ''
  """
    ${dereplicate_cmd}

    bigscape cluster \\
        -i ${cluster_input} \\
        -o bigscape \\
        -c ${task.cpus} \\
        ${pfam_arg} \\
        --force-gbk \\
        --exclude-gbk final,deepbgc \\
        ${args}

    ${archive_cmd}
    """

  stub:
  """
    mkdir -p bigscape_dereplicate/representative_clusters
    mkdir -p bigscape
    touch bigscape/bigscape.db
    """
}
