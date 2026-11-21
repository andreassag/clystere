/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMBGC_FILTER
    Runs comBGC.py on antiSMASH + GECCO + deepBGC outputs and stages
    sample-prefixed region GBKs for downstream clustering.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COMBGC_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'
        : 'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0'}"

    input:
    tuple val(meta), path(antismash_gbk, stageAs: 'antismash_input/*'), path(antismash_region_gbks, stageAs: 'antismash_regions/*'), path(knownclusterblast_dir, stageAs: 'knownclusterblast'), path(gecco_clusters_tsv, stageAs: 'gecco_tsv/*'), path(gecco_gbks, stageAs: 'gecco_gbks/*'), path(deepbgc_bgc_gbk, stageAs: 'deepbgc_gbk/*'), path(deepbgc_bgc_tsv, stageAs: 'deepbgc_tsv/*')

    output:
    tuple val(meta), path("${prefix}/combgc_summary.tsv"), emit: summary_tsv
    tuple val(meta), path("${prefix}/${prefix}_regions/"), emit: combined_regions_dir
    tuple val(meta), path("${prefix}/${prefix}_regions/*.gbk"), optional: true, emit: combined_region_gbks
    tuple val("${task.process}"), val('combgc'), eval("comBGC.py --version 2>&1 | sed 's/comBGC //'"), topic: versions, emit: versions_combgc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    shopt -s nullglob
    input_args=()

    if [[ -f "antismash_input/${prefix}.gbk" ]]; then
        input_args+=("antismash_input/${prefix}.gbk")
    else
        for f in antismash_input/*.gbk; do
            if [[ ! "\${f}" =~ region[0-9]{3}\\.gbk\$ ]]; then
                input_args+=("\${f}")
            fi
        done
    fi

    if [[ -d "knownclusterblast" ]]; then
        input_args+=("knownclusterblast/")
    fi

    for f in gecco_tsv/*.tsv; do
        input_args+=("\${f}")
    done

    for f in gecco_gbks/*.gbk; do
        input_args+=("\${f}")
    done

    for f in deepbgc_tsv/*.tsv; do
        input_args+=("\${f}")
    done

    mkdir -p "${prefix}/${prefix}_regions"

    if [[ \${#input_args[@]} -gt 0 ]]; then
        comBGC.py \\
            ${args} \\
            --outdir "${prefix}" \\
            --input "\${input_args[@]}"
    else
        printf 'sample_id\\tcontig_id\\tPrediction_tool\\tProduct_class\\tBGC_probability\\tBGC_complete\\tBGC_start\\tBGC_end\\tBGC_length\\tBGC_region_contig_ids\\tCDS_count\\tPFAM_domains\\tMIBiG_ID\\tInterPro_ID\\n' > "${prefix}/combgc_summary.tsv"
    fi

    for gbk in antismash_regions/*.gbk; do
        base=\$(basename "\${gbk}")
        if [[ "\${base}" == "${prefix}"* ]]; then
            cp "\${gbk}" "${prefix}/${prefix}_regions/\${base}"
        else
            cp "\${gbk}" "${prefix}/${prefix}_regions/${prefix}_\${base}"
        fi
    done

    python3 - "${prefix}/${prefix}_regions" "${prefix}" <<'PY'
import re
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

outdir, sample_prefix = Path(sys.argv[1]), sys.argv[2]

for src in sorted(Path("gecco_gbks").glob("*.gbk")):
    base = src.name
    if base.startswith(sample_prefix):
        dest = outdir / base
    else:
        dest = outdir / f"{sample_prefix}_{base}"
    if re.search(r"[.]region[0-9]{3}[.]gbk\$", base) or dest.exists():
        dest = outdir / f"{sample_prefix}.gecco.{base}"
    if src.stat().st_size == 0:
        continue
    records = list(SeqIO.parse(src, "genbank"))
    for rec in records:
        sc = rec.annotations.get("structured_comment", {})
        gecco_meta = sc.get("GECCO-Data", {})
        product = gecco_meta.get("cluster_type") or "other"
        if "antiSMASH-Data" in sc and not any(f.type == "region" for f in rec.features):
            sc.pop("antiSMASH-Data", None)
        clusters = [f for f in rec.features if f.type == "cluster"]
        if not clusters:
            cluster_feat = SeqFeature(
                FeatureLocation(0, len(rec.seq)),
                type="cluster",
                qualifiers={
                    "note": ["Cluster number: 1"],
                    "product": [product],
                    "contig_edge": ["False"],
                },
            )
            rec.features.insert(0, cluster_feat)
        else:
            for c_feat in clusters:
                c_feat.qualifiers["note"] = ["Cluster number: 1"]
                c_feat.qualifiers.setdefault("product", [product])
                c_feat.qualifiers.setdefault("contig_edge", ["False"])
    SeqIO.write(records, dest, "genbank")

for src in sorted(Path("deepbgc_gbk").glob("*.gbk")):
    if src.stat().st_size == 0:
        continue
    for idx, rec in enumerate(SeqIO.parse(src, "genbank"), start=1):
        clusters = [f for f in rec.features if f.type == "cluster"]
        if not clusters:
            cluster_feat = SeqFeature(
                FeatureLocation(0, len(rec.seq)),
                type="cluster",
                qualifiers={
                    "note": [f"Cluster number: {idx}"],
                    "product": ["other"],
                    "contig_edge": ["False"],
                },
            )
            rec.features.insert(0, cluster_feat)
        else:
            for c_feat in clusters:
                c_feat.qualifiers["note"] = [f"Cluster number: {idx}"]
                if "product_class" in c_feat.qualifiers:
                    c_feat.qualifiers["product"] = c_feat.qualifiers["product_class"]
                c_feat.qualifiers.setdefault("product", ["other"])
                c_feat.qualifiers.setdefault("contig_edge", ["False"])
        out = outdir / f"{sample_prefix}.deepbgc.region{idx:03d}.gbk"
        SeqIO.write(rec, out, "genbank")
PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/${prefix}_regions"

    printf 'sample_id\\tcontig_id\\tPrediction_tool\\tProduct_class\\tBGC_probability\\tBGC_complete\\tBGC_start\\tBGC_end\\tBGC_length\\tBGC_region_contig_ids\\tCDS_count\\tPFAM_domains\\tMIBiG_ID\\tInterPro_ID\\n' > "${prefix}/combgc_summary.tsv"
    printf '${prefix}\\t${prefix}\\tantiSMASH\\tNRPS\\tNA\\tYes\\t1\\t100\\t100\\tctg1_1\\t1\\tAMP-binding\\tNA\\tNA\\n' >> "${prefix}/combgc_summary.tsv"
    printf 'LOCUS       stub                    100 bp    DNA     linear   BCT 01-JAN-2026\\n//\\n' > "${prefix}/${prefix}_regions/${prefix}.region001.gbk"
    """
}
