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
  tuple val(meta), path(antismash_dir), path(gecco_bigslice_dir), path(deepbgc_bigslice_dir), path(gecco_clusters_tsv), path(deepbgc_bgc_tsv)

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

    python select_combined_regions.py \\
        --summary combgc_out/combgc_summary.tsv \\
        --sample-id ${prefix} \\
        --antismash-dir ${antismash_dir} \\
        --gecco-dir ${gecco_bigslice_dir} \\
        --deepbgc-dir ${deepbgc_bigslice_dir} \\
          --out-dir combgc_out/${prefix}_regions

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
