/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COUNT_REGIONS
    Produces a per-genome BGC count TSV across all antiSMASH results.
    All per-sample antiSMASH output directories must be collected and passed together
    so the script can glob *\/*.json across all of them.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COUNT_REGIONS {
  label 'process_low'
  conda 'conda-forge::python=3.11'
  container 'python:3.11'
  input:
  path(antismash_dirs)
  output:
  path "region_counts.tsv", emit: tsv
  path "versions.yml", emit: versions
  when:
  task.ext.when == null || task.ext.when
  script:
  def by_contig_arg = params.count_per_contig ? '--by_contig': ''
  def split_hybrid_arg = params.split_hybrids ? '--split_hybrids': ''
  """
    count_regions.py . region_counts.tsv \\
        ${by_contig_arg} \\
        ${split_hybrid_arg} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
  stub:
  """
    printf 'sample\tcount\n' > region_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
