/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SANNTIS
    Runs SanntiS neural-network BGC predictor on CDS-annotated GenBank genomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SANNTIS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sanntis:0.9.3.5--pyhdfd78af_0'
        : 'quay.io/biocontainers/sanntis:0.9.3.5--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(gbk), path(ip_file)

    output:
    tuple val(meta), path("${prefix}/${prefix}.sanntis.full.gff"), emit: gff
    tuple val(meta), path("${prefix}"), emit: outdir
    tuple val("${task.process}"), val('sanntis'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"sanntis\"))'"), topic: versions, emit: versions_sanntis

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ip_arg = ip_file ? "--ip-file \"${ip_file}\"" : ''
    """
    if [[ -z "${ip_file ?: ''}" ]] && ! command -v interproscan.sh >/dev/null 2>&1; then
        echo "ERROR: SanntiS requires InterProScan (interproscan.sh in PATH) or pre-computed InterProScan GFF3/TSV files via --sanntis_interpro_dir <dir>." >&2
        echo "The standard SanntiS BioContainer does not bundle InterProScan. Provide per-sample InterProScan annotations matching <sample_id>* in --sanntis_interpro_dir." >&2
        exit 1
    fi

    mkdir -p "${prefix}"
    sanntis \\
        --cpu ${task.cpus} \\
        --outdir "${prefix}" \\
        --outfile "${prefix}/${prefix}.sanntis.full.gff" \\
        ${ip_arg} \\
        ${args} \\
        "${gbk}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}"
    printf '##gff-version 3\\n' > "${prefix}/${prefix}.sanntis.full.gff"
    printf 'ctg1\\tSanntiSv0.9.3.5\\tCLUSTER\\t1\\t100\\t.\\t.\\t.\\tID=${prefix}_sanntis_1;nearest_MiBIG=BGC0000349;nearest_MiBIG_class=NRP;nearest_MiBIG_diceDistance=0.150;score=0.920;partial=00\\n' >> "${prefix}/${prefix}.sanntis.full.gff"
    """
}
