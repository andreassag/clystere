/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GUNZIP_GFF
    Decompresses gzipped GFF3 annotations (e.g., from Pyrodigal or user input) into
    plain-text .gff3 files required by antiSMASH 8 --genefinding-gff3.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GUNZIP_GFF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    tuple val(meta), path(gff_in, stageAs: 'input_gff/*')

    output:
    tuple val(meta), path("${prefix}.gff3"), emit: gff
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_gunzip_gff

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${gff_in}" "${prefix}.gff3" <<'PY'
import gzip
import shutil
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
if src.suffix == ".gz":
    with gzip.open(src, "rb") as fin, dst.open("wb") as fout:
        shutil.copyfileobj(fin, fout)
else:
    shutil.copyfile(src, dst)
PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '##gff-version 3\\n' > "${prefix}.gff3"
    """
}
