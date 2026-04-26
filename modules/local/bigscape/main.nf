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
  container 'quay.io/biocontainers/bigscape:2.0.2--pyhdfd78af_0'

  input:
  path bgc_dirs
  path pfam

  output:
  path "bigscape_dereplicate/", emit: dereplicate_dir, optional: true
  path "bigscape/bigscape.db", emit: db
  path "bigscape/", emit: results_dir
  path "bigscape.tar.gz", emit: archive, optional: true
  path "versions.yml", emit: versions

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
        --cutoff ${params.bigscape_dereplicate_cutoff}
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
        ${args}

    ${archive_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: \$(bigscape --version 2>&1 | head -1 | sed 's/bigscape //')
    END_VERSIONS
    """

  stub:
  """
    mkdir -p bigscape_dereplicate/representative_clusters
    mkdir -p bigscape
    touch bigscape/bigscape.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: 2.0.2
    END_VERSIONS
    """
}
