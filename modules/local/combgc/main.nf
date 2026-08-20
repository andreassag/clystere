/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMBGC_FILTER
    Runs comBGC on antiSMASH + GECCO + deepBGC outputs and selects
    representative antiSMASH-like region GBKs for downstream clustering.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COMBGC_FILTER {
  tag "${meta.id}"
  label 'process_medium'
  conda "conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2"
  container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
      'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

  input:
  tuple val(meta),
    path(antismash_dir, stageAs: 'antismash_dir'),
    path(gecco_dir, stageAs: 'gecco_dir'),
    path(deepbgc_bgc_gbk, stageAs: 'deepbgc_bgc.gbk'),
    path(gecco_clusters_tsv, stageAs: 'gecco.clusters.tsv'),
    path(deepbgc_bgc_tsv, stageAs: 'deepbgc.bgc.tsv')

  output:
  tuple val(meta), path("combgc_out/combgc_summary.tsv"), emit: summary_tsv
  tuple val(meta), path("combgc_out/${meta.id}_regions/"), emit: combined_regions_dir
  tuple val("${task.process}"), val('combgc'), val("0.6.9"), topic: versions, emit: versions_combgc

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = meta.id
  """
    input_args=()

    if [[ -f "${antismash_dir}/${prefix}.gbk" ]]; then
      input_args+=("${antismash_dir}/${prefix}.gbk")
    fi

    if [[ -f "${gecco_clusters_tsv}" ]]; then
      input_args+=("${gecco_clusters_tsv}")
    fi

    if [[ -f "${deepbgc_bgc_tsv}" ]]; then
      input_args+=("${deepbgc_bgc_tsv}")
    fi

    if [[ -d "${antismash_dir}/knownclusterblast" ]]; then
      input_args+=("${antismash_dir}/knownclusterblast/")
    fi

    combgc \\
        ${args} \\
        --outdir combgc_out \\
        --input "\${input_args[@]}"

    mkdir -p combgc_out/${prefix}_regions

    # Build a unified region directory using only native module outputs.
    shopt -s nullglob
    for gbk in "${antismash_dir}"/*region*.gbk; do
      cp "\${gbk}" "combgc_out/${prefix}_regions/"
    done

    for gbk in "${gecco_dir}"/*_cluster_*.gbk "${gecco_dir}"/*.region*.gbk; do
      if [[ -f "\${gbk}" ]]; then
        cp "\${gbk}" "combgc_out/${prefix}_regions/"
      fi
    done

    if [[ -f "${deepbgc_bgc_gbk}" ]]; then
      cp "${deepbgc_bgc_gbk}" "combgc_out/${prefix}_regions/${prefix}.deepbgc.region001.gbk"
    fi

    if [[ -z "\$(ls -A combgc_out/${prefix}_regions)" ]]; then
      touch "combgc_out/${prefix}_regions/${prefix}.region001.gbk"
    fi
    """

  stub:
  def prefix = meta.id
  """
    mkdir -p combgc_out/${prefix}_regions

    printf 'sample_id\tcontig_id\tBGC_start\tBGC_end\tTool_representative\\n' > combgc_out/combgc_summary.tsv
    printf '${prefix}\t${prefix}\t1\t100\tantiSMASH\\n' >> combgc_out/combgc_summary.tsv
    touch combgc_out/${prefix}_regions/${prefix}.region001.gbk
    """
}
