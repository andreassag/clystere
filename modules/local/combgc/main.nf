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
  conda "${moduleDir}/environment.yml"
  container 'quay.io/biocontainers/python:3.11.12--h9e4cc4f_0'

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
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = meta.id
  """
    if ! command -v combgc >/dev/null 2>&1; then
      python -m pip install --no-cache-dir "combgc==0.6.9"
    fi

    input_args=(
      "${antismash_dir}/${prefix}.gbk"
      "${gecco_clusters_tsv}"
      "${deepbgc_bgc_tsv}"
    )

    if [[ -d "${antismash_dir}/knownclusterblast" ]]; then
      input_args+=("${antismash_dir}/knownclusterblast/")
    fi

    combgc \\
        --cores ${task.cpus} \\
      ${args} \\
          --outdir combgc_out \\
        --input "\${input_args[@]}"

    mkdir -p combgc_out/${prefix}_regions

    # Build a unified region directory using only native module outputs.
    shopt -s nullglob
    for gbk in "${antismash_dir}"/*region*.gbk; do
      cp "\${gbk}" "combgc_out/${prefix}_regions/"
    done

    for gbk in "${gecco_dir}"/*_cluster_*.gbk; do
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combgc: \$(python -c "import importlib.metadata as m; print(m.version('combgc'))")
    END_VERSIONS
    """

  stub:
  def prefix = meta.id
  """
    mkdir -p combgc_out/${prefix}_regions

    printf 'sample_id\tcontig_id\tBGC_start\tBGC_end\tTool_representative\n' > combgc_out/combgc_summary.tsv
    printf '${prefix}\t${prefix}\t1\t100\tantiSMASH\n' >> combgc_out/combgc_summary.tsv
      touch combgc_out/${prefix}_regions/${prefix}.region001.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combgc: 0.6.9
    END_VERSIONS
    """
}
