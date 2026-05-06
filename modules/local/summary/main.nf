/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUMMARY
    Produces both per-region and per-genome summary TSVs across all antiSMASH results.
    All per-sample antiSMASH output directories must be collected and passed together
    so the summary scripts can glob *\/*.json across all of them.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SUMMARY {
    label 'process_low'
    conda 'conda-forge::python=3.11'
    container 'python:3.11'

    input:
    path antismash_dirs

    output:
    path 'all_regions.tsv', emit: all_regions_tsv
    path 'region_counts.tsv', emit: region_counts_tsv
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def knownclusters_arg = params.antismash_cb_knownclusters ? '--knownclusters' : ''
    def by_contig_arg = params.count_per_contig ? '--by_contig' : ''
    def split_hybrid_arg = params.split_hybrids ? '--split_hybrids' : ''
    """
    tabulate_regions.py . all_regions.tsv \
        ${knownclusters_arg} \
        --threads ${task.cpus}

    count_regions.py . region_counts.tsv \
        ${by_contig_arg} \
        ${split_hybrid_arg} \
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\tregion\tproduct\n' > all_regions.tsv
    printf 'sample\tcount\n' > region_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
