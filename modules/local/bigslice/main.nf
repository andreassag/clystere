/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIGSLICE
    Runs BiG-SLiCE clustering across all antiSMASH output directories.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BIGSLICE {
  label 'process_high'
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bigslice%3A2.0.2--pyh8ed023e_0' :
        'quay.io/biocontainers/bigslice:2.0.2--pyh8ed023e_0' }"
  input:
  path(antismash_dirs)
  output:
  path "bigslice/", emit: results_dir
  path "bigslice.tar.gz", emit: archive, optional: true
  path "versions.yml", emit: versions
  when:
  task.ext.when == null || task.ext.when
  script:
  def args = task.ext.args ?: ''
  def model_url = 'https://github.com/medema-group/bigslice/releases/download/v2.0.0rc/bigslice-models.2022-11-30.tar.gz'
  def model_md5 = 'aaabde911ec107d08e5c24f68aaf31d1'
  def archive_cmd = params.bigslice_zip_output
 ? "tar --exclude 'bigslice/cache' -zcf bigslice.tar.gz bigslice/"
: ''
  """
    model_dir="\$PWD/bigslice-models"
    model_tar="\$PWD/bigslice_models.tar.gz"

    # Ensure BiG-SLiCE HMM models are available in the task work directory.
    # The process then points --program_db_folder to this local model directory.
    if [[ ! -s "\$model_dir/biosynthetic_pfams/biopfam.md5sum" ]]; then
      rm -rf "\$model_dir"
      if [[ ! -s "\$model_tar" ]]; then
        python - <<'PY'
import urllib.request
urllib.request.urlretrieve('${model_url}', 'bigslice_models.tar.gz')
PY
      fi
      echo "${model_md5}  \$model_tar" | md5sum -c -
      mkdir -p "\$model_dir"
      tar -xzf "\$model_tar" -C "\$model_dir"
    fi

    printf 'antismash\t.\t\tantiSMASH outputs\n' > datasets.tsv

    bigslice -i . --program_db_folder "\$model_dir" bigslice ${args}

    ${archive_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigslice: \$(bigslice --version 2>&1 | grep -m1 -Eo '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """
  stub:
  """
    mkdir -p bigslice

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigslice: 2.0.2
    END_VERSIONS
    """
}
