/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GBK_TO_FASTA
    Extracts nucleotide FASTA (and optional CDS FAA) from GenBank files or normalizes
    FASTA inputs to a uniform uncompressed <sample>.fna filename.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GBK_TO_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'
        : 'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'}"

    input:
    tuple val(meta), path(genome, stageAs: 'input_genome/*')
    val is_fasta

    output:
    tuple val(meta), path("${prefix}.fna"), emit: fasta
    tuple val(meta), path("${prefix}.faa"), optional: true, emit: faa
    tuple val("${task.process}"), val('biopython'), eval("python3 -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def mode_arg = is_fasta ? '--normalize-fasta' : "--faa ${prefix}.faa"
    """
    gbk_to_fasta.py "${genome}" "${prefix}.fna" ${mode_arg}
    if [[ -f "${prefix}.faa" && ! -s "${prefix}.faa" ]]; then
        rm -f "${prefix}.faa"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>stub_contig\\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG\\n' > "${prefix}.fna"
    """
}
