/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TABULATE_REGIONS
    Produces a per-region TSV across all antiSMASH results.
    All per-sample antiSMASH output directories must be collected and passed together
    so the script can glob *\/*.json across all of them.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TABULATE_REGIONS {
  label 'process_low'
  conda 'conda-forge::python=3.11'
  container 'python:3.11'
  input:
  path(antismash_dirs)
  output:
  path "all_regions.tsv", emit: tsv
  path "versions.yml", emit: versions
  when:
  task.ext.when == null || task.ext.when
  script:
  def knownclusters_arg = params.antismash_cb_knownclusters ? '--knownclusters': ''
  """
    tabulate_regions.py . all_regions.tsv \\
        ${knownclusters_arg} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}