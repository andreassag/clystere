/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIGSCAPE
    Runs BiG-SCAPE v2 (bigscape cluster) across all antiSMASH output directories.
    All per-sample antiSMASH output directories are passed as a collected list and
    staged into the work directory, so BiG-SCAPE can scan them recursively.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BIGSCAPE {
  label 'process_high'
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bigscape:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/bigscape:2.0.2--pyhdfd78af_0' }"
  input:
  path(antismash_dirs)
  path(pfam)
  output:
  path "bigscape/bigscape.db", emit: db
  path "bigscape/", emit: results_dir
  path "bigscape.tar.gz", emit: archive, optional: true
  path "versions.yml", emit: versions
  when:
  task.ext.when == null || task.ext.when
  script:
  def args = task.ext.args ?: ''
  def pfam_arg = pfam
 ? (pfam.isDirectory() ? "-p ${pfam}/Pfam-A.hmm": "-p ${pfam}")
: ''
  def archive_cmd = params.bigscape_zip_output
 ? "tar --exclude 'bigscape/cache' -zcf bigscape.tar.gz bigscape/"
: ''
  """
    bigscape cluster \\
        -i . \\
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
    mkdir -p bigscape
    touch bigscape/bigscape.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: 2.0.2
    END_VERSIONS
    """
}
