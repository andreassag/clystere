/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEEPBGC process
    Runs deepBGC v0.1.31 on a single genome file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEEPBGC_PIPELINE {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'

    input:
    tuple val(meta), path(genome), path(annotation), path(downloads_dir)

    output:
    tuple val(meta), path("${meta.id}/"), emit: output_dir
    tuple val(meta), path("${meta.id}/${meta.id}.bgc.tsv"), emit: bgc_tsv
    tuple val(meta), path("${meta.id}/deepbgc_bigslice/"), emit: bigslice_dir, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def convert_cmd = params.bigscape_run || params.bigslice_run
        ? """
    mkdir -p ${prefix}/deepbgc_bigslice
    python deepbgc_to_antismash_region_gbk.py \\
      --full-gbk ${prefix}/${prefix}.full.gbk \\
      --bgc-tsv ${prefix}/${prefix}.bgc.tsv \\
      --output-dir ${prefix}/deepbgc_bigslice \\
      --prefix ${prefix}
  """
        : ''
    """
    export DEEPBGC_DOWNLOADS_DIR="${downloads_dir}"

    deepbgc pipeline \\
        --output ${prefix} \\
        ${args} \\
        ${genome}

    if [[ ! -f "${prefix}/${prefix}.bgc.tsv" ]]; then
      printf 'sequence_id\tdetector\tdetector_version\tdetector_label\tbgc_candidate_id\tnucl_start\tnucl_end\tnucl_length\tnum_proteins\tnum_domains\tnum_bio_domains\tdeepbgc_score\tproduct_activity\tantibacterial\tcytotoxic\tinhibitor\tantifungal\tproduct_class\tAlkaloid\tNRP\tOther\tPolyketide\tRiPP\tSaccharide\tTerpene\tprotein_ids\tbio_pfam_ids\tpfam_ids\n' > ${prefix}/${prefix}.bgc.tsv
    fi

    ${convert_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      deepbgc: \$(deepbgc info 2>&1 | sed -n '6p' | sed 's/.*= version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p ${prefix}
    printf 'sequence_id\tdetector\tdetector_version\tdetector_label\tbgc_candidate_id\tnucl_start\tnucl_end\tnucl_length\tnum_proteins\tnum_domains\tnum_bio_domains\tdeepbgc_score\tproduct_activity\tantibacterial\tcytotoxic\tinhibitor\tantifungal\tproduct_class\tAlkaloid\tNRP\tOther\tPolyketide\tRiPP\tSaccharide\tTerpene\tprotein_ids\tbio_pfam_ids\tpfam_ids\n' > ${prefix}/${prefix}.bgc.tsv
    mkdir -p ${prefix}/deepbgc_bigslice
    touch ${prefix}/deepbgc_bigslice/${prefix}.region001.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: 0.1.31
    END_VERSIONS
    """
}
