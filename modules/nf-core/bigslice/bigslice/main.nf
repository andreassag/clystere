process BIGSLICE_BIGSLICE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0'
        : 'quay.io/biocontainers/bigslice:2.0.2--pyh8ed023e_0'}"

    input:
    tuple val(meta), path(bgc, stageAs: 'bgc_files/s*/*')
    path hmmdb
    val export_tsv
    path taxonomy_tsv

    output:
    tuple val(meta), path("${prefix}/result"), emit: output
    tuple val(meta), path("${prefix}/result/tsv_export"), emit: tsv, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('bigslice'), val("2.0.2"), topic: versions, emit: versions_bigslice

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample = meta.id
    def export_tsv_cmd = export_tsv ? "bigslice --export-tsv ${prefix}/result/tsv_export --program_db_folder ${hmmdb} ${args2} ${prefix}" : ''
    def taxonomy_cmd = taxonomy_tsv ? "cp ${taxonomy_tsv} input/taxonomy/taxonomy.tsv" : "touch input/taxonomy/taxonomy.tsv"
    """
    mkdir -p input/dataset input/taxonomy
    shopt -s nullglob
    staged_dirs=0
    for d in bgc_files/s*/*_regions; do
        if [[ -d "\${d}" ]]; then
            genome_id="\$(basename "\${d}" _regions)"
            mkdir -p "input/dataset/\${genome_id}"
            for g in "\${d}"/*.gbk; do
                cp "\${g}" "input/dataset/\${genome_id}/"
            done
            staged_dirs=1
        fi
    done
    if [[ \${staged_dirs} -eq 0 ]]; then
        mkdir -p "input/dataset/${sample}"
        find -L bgc_files -name '*.gbk' -exec cp {} "input/dataset/${sample}/" \\;
    fi

    printf "# dataset_name\\tdataset_path\\ttaxonomy_path\\tdescription\\n" > input/datasets.tsv
    printf "dataset\\tdataset\\ttaxonomy/taxonomy.tsv\\tBGC dataset\\n" >> input/datasets.tsv

    ${taxonomy_cmd}

    bigslice \\
        ${args} \\
        --num_threads ${task.cpus} \\
        -i input \\
        --program_db_folder ${hmmdb} \\
        ${prefix}

    ${export_tsv_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}/result/tmp/2e555308dfc411186cf012334262f127
    touch ${prefix}/result/data.db
    touch ${prefix}/result/tmp/2e555308dfc411186cf012334262f127/test.fa
    if ${export_tsv}; then
        mkdir -p ${prefix}/result/tsv_export
        touch ${prefix}/result/tsv_export/bgcs.tsv
    fi
    """
}
