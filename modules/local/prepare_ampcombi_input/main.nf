/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_AMPCOMBI_INPUT
    Decompresses Macrel prediction and protein ORF FASTA files to plain text with
    sample-prefixed filenames expected by AMPCOMBI2_PARSETABLES.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PREPARE_AMPCOMBI_INPUT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'
        : 'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'}"

    input:
    tuple val(meta), path(prediction_gz), path(orfs_faa_gz), path(genome_fasta, stageAs: 'genome_in/*'), path(antismash_gbk, stageAs: 'antismash_in/*')

    output:
    tuple val(meta), path("${prefix}.macrel.prediction"), emit: prediction
    tuple val(meta), path("${prefix}.faa"), emit: faa
    tuple val(meta), path("${prefix}.gbk"), emit: gbk
    tuple val("${task.process}"), val('biopython'), eval("python3 -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_prepare_ampcombi

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def gbk_arg = antismash_gbk ? "--antismash-gbk \"${antismash_gbk}\"" : ''
    """
    prepare_ampcombi_inputs.py \\
        --prediction "${prediction_gz}" \\
        --orfs-faa "${orfs_faa_gz}" \\
        --fasta "${genome_fasta}" \\
        ${gbk_arg} \\
        --out-prediction "${prefix}.macrel.prediction" \\
        --out-faa "${prefix}.faa" \\
        --out-gbk "${prefix}.gbk"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '# Prediction from macrel v1.4.0\\n' > "${prefix}.macrel.prediction"
    printf 'Access\\tSequence\\tAMP_family\\tAMP_probability\\tHemolytic\\tHemolytic_probability\\n' >> "${prefix}.macrel.prediction"
    printf '>stub_contig_1 # 1 # 30 # 1 # ID=1_1\\nMKKLLVVAA\\n' > "${prefix}.faa"
    printf 'LOCUS       stub_contig              100 bp    DNA     linear   BCT 01-JAN-2026\\n//\\n' > "${prefix}.gbk"
    """
}
