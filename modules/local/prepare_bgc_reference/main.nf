/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_BGC_REFERENCE
    Extracts non-redundant BGC nucleotide FASTA sequences and catalog metadata for
    metagenomic / metatranscriptomic read mapping.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PREPARE_BGC_REFERENCE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'
        : 'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'}"

    input:
    path bgc_dirs, stageAs: 'bgc_inputs/?/*'
    path bigscape_derep_dir
    path gcf_membership_tsv

    output:
    tuple val(meta), path('bgc_catalog.fna'), emit: fasta
    path 'bgc_catalog_metadata.tsv', emit: metadata_tsv
    tuple val("${task.process}"), val('biopython'), eval("python3 -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_prepare_bgc_ref

    when:
    task.ext.when == null || task.ext.when

    script:
    meta = [id: 'bgc_catalog']
    def derep_arg = bigscape_derep_dir ? "--derep-dir ${bigscape_derep_dir}" : ''
    def gcf_arg = gcf_membership_tsv ? "--gcf-membership ${gcf_membership_tsv}" : ''
    """
    prepare_bgc_reference.py \\
        --gbks bgc_inputs \\
        ${derep_arg} \\
        ${gcf_arg} \\
        --out-fasta bgc_catalog.fna \\
        --out-metadata bgc_catalog_metadata.tsv
    """

    stub:
    meta = [id: 'bgc_catalog']
    """
    printf '>sample1.region001 sample=sample1 gcf=NA product=NRPS\\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\\n' > bgc_catalog.fna
    printf 'bgc_id\\tsample_id\\tsource_file\\tlength_bp\\tproduct\\tgcf_id\\n' > bgc_catalog_metadata.tsv
    printf 'sample1.region001\\tsample1\\tsample1.region001.gbk\\t60\\tNRPS\\tNA\\n' >> bgc_catalog_metadata.tsv
    """
}
