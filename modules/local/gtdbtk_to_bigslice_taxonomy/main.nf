/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GTDBTK_TO_BIGSLICE_TAXONOMY
    Converts GTDB-Tk summary TSV(s) into BiG-SLiCE's 9-column taxonomy.tsv format.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GTDBTK_TO_BIGSLICE_TAXONOMY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    path summaries, stageAs: 'gtdbtk_summaries/?/*'

    output:
    path 'taxonomy.tsv', emit: taxonomy_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_gtdbtk_to_bigslice

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gtdbtk_to_bigslice_taxonomy.py gtdbtk_summaries/*/* -o taxonomy.tsv
    """

    stub:
    """
    printf '# path_startswith\\tKingdom\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus\\tSpecies\\tOrganism\\n' > taxonomy.tsv
    """
}
