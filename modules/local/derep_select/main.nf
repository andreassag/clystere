/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEREP_SELECT
    Selects non-redundant representative genomes from skani triangle ANI output.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEREP_SELECT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'python:3.11'}"

    input:
    tuple val(meta), path(skani_tsv)
    path genomes, stageAs: 'genomes/*'
    val ani_threshold
    path checkm2_tsv

    output:
    path 'dereplicated_genome_ids.txt', emit: ids
    path 'dereplication_clusters.tsv', emit: clusters_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_derep_select

    when:
    task.ext.when == null || task.ext.when

    script:
    def checkm2_arg = checkm2_tsv ? "--checkm2 \"${checkm2_tsv}\"" : ''
    """
    select_derep_genomes.py \\
        --skani "${skani_tsv}" \\
        --genomes genomes/* \\
        --ani "${ani_threshold}" \\
        ${checkm2_arg} \\
        --out-ids dereplicated_genome_ids.txt \\
        --out-clusters dereplication_clusters.tsv
    """

    stub:
    """
    python3 - genomes/* <<'PY'
import sys
from pathlib import Path

ids = []
for raw in sorted(sys.argv[1:]):
    name = Path(raw).name
    for suf in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.lower().endswith(suf):
            name = name[:-len(suf)]
            break
    ids.append(name)

Path("dereplicated_genome_ids.txt").write_text("\\n".join(ids) + ("\\n" if ids else ""))
with Path("dereplication_clusters.tsv").open("w") as fh:
    fh.write("cluster_id\\tsample_id\\tis_representative\\tscore\\tcompleteness\\tcontamination\\tn50\\ttotal_length\\tani_to_rep\\n")
    for idx, sid in enumerate(ids, start=1):
        fh.write(f"DEREP_{idx:04d}\\t{sid}\\ttrue\\t1.0\\tNA\\tNA\\t1000\\t1000\\t100.0\\n")
PY
    """
}
