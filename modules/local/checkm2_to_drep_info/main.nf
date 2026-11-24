/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKM2_TO_DREP_INFO
    Converts CheckM2 quality_report.tsv into dRep --genomeInfo CSV format.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHECKM2_TO_DREP_INFO {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    tuple val(meta), path(checkm2_tsv)
    path genomes, stageAs: 'fastas/*'

    output:
    path 'drep_genome_info.csv', emit: genome_info_csv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_checkm2_to_drep_info

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    checkm2_to_drep_info.py \\
        --checkm2 "${checkm2_tsv}" \\
        --genomes fastas/* \\
        --out drep_genome_info.csv
    """

    stub:
    """
    python3 - fastas/* <<'PY'
    import csv
    import sys
    from pathlib import Path

    with Path("drep_genome_info.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\\n")
        writer.writerow(["genome", "completeness", "contamination"])
        for raw in sorted(sys.argv[1:]):
            writer.writerow([Path(raw).name, "99.00", "0.50"])
    PY
    """
}
