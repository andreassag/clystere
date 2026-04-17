/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANTISMASH process
    Runs antiSMASH v8 on a single genome file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANTISMASH {
  tag "${meta.id}"
  label 'process_medium'
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash:8.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash:8.0.1--pyhdfd78af_0' }"
  input:
  tuple val(meta), path(genome), path(annotation), path(databases)
  output:
  tuple val(meta), path("${meta.id}/"), emit: output_dir
  tuple val(meta), path("${meta.id}/${meta.id}.json"), emit: json_results
  tuple val(meta), path("${meta.id}/index.html"), emit: html, optional: true
  tuple val(meta), path("${meta.id}/${meta.id}.zip"), emit: zip, optional: true
  tuple val(meta), path("${meta.id}/region*.gbk"), emit: gbk_results, optional: true
  tuple val(meta), path("${meta.id}/knownclusterblast/"), emit: knownclusterblast_dir, optional: true
  tuple val(meta), path("${meta.id}/clusterblast/"), emit: clusterblast_dir, optional: true
  path "versions.yml", emit: versions
  when:
  task.ext.when == null || task.ext.when
  script:
  def args = task.ext.args ?: ''
  def prefix = meta.id
  def annotation_arg = annotation ? "--genefinding-gff3 ${annotation}": ''
  def reuse_arg = params.antismash_reuse_results ? '--reuse-results': ''
  def genefinding_arg = annotation ? '': "--genefinding-tool ${params.antismash_genefinding_tool}"
  def accept_fail = params.antismash_accept_failure
  // Detect pre-annotated input formats (GenBank / EMBL) when no GFF3 is supplied
  if (!annotation) {
    def ext = genome.name.replaceAll(/\.gz$/, '').tokenize('.')[-1]
    if (ext in ['gbk', 'gb', 'gbff', 'gbf', 'embl', 'emb']) {
      genefinding_arg = '--genefinding-tool none'
    }
  }
  """
    antismash \\
        --output-dir ${prefix} \\
        --output-basename ${prefix} \\
        --databases ${databases} \\
        -c ${task.cpus} \\
        --logfile ${prefix}/${prefix}.log \\
        ${genefinding_arg} \\
        ${annotation_arg} \\
        ${reuse_arg} \\
        ${args} \\
        ${genome} \\
        || { rc=\$?
             if [ "${accept_fail}" = "true" ]; then
                 mkdir -p ${prefix}
                 echo '{}' > ${prefix}/${prefix}.json
                 echo "WARNING: antiSMASH failed for ${prefix} (exit \$rc). Continuing due to --antismash_accept_failure." >&2
             else
                 exit \$rc
             fi
           }

    # make outputs readable by host user
    chmod -R a+rX ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version 2>&1 | head -1 | sed 's/antiSMASH //')
    END_VERSIONS
    """
}