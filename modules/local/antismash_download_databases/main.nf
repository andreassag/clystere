process ANTISMASH_DOWNLOAD_DATABASES {
  label 'process_single'
  conda "${moduleDir}/../antismash/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash:8.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash:8.0.1--pyhdfd78af_0' }"
  input:
  val(db_dest)
  output:
  val(db_dest), emit: databases
  path "versions.yml", emit: versions
  when:
  task.ext.when == null || task.ext.when
  script:
  def args = task.ext.args ?: ''
  """
    download-antismash-databases \\
        --database-dir ${db_dest} \\
        ${args}

    chmod -R a+rX ${db_dest}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version 2>&1 | head -1 | sed 's/antiSMASH //')
    END_VERSIONS
    """
  stub:
  """
    mkdir -p ${db_dest}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: 8.0.1
    END_VERSIONS
    """
}
