/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    clystere Main Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome preprocessing & QC/taxonomy -> antiSMASH + GECCO + deepBGC + SanntiS + AMPs
    -> comBGC + BiG-SCAPE / BiG-SLiCE / BiG-FAM -> clinker / cblaster
    -> Metagenomic/Metatranscriptomic BGC abundance -> Consensus Summary & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPROCESSING   } from '../subworkflows/local/genome_preprocessing'
include { GENOME_QC_TAXONOMY     } from '../subworkflows/local/genome_qc_taxonomy'
include { BGC_ANNOTATION         } from '../subworkflows/local/bgc_annotation'
include { BGC_CLUSTERING         } from '../subworkflows/local/bgc_clustering'
include { BGC_COMPARATIVE        } from '../subworkflows/local/bgc_comparative'
include { BGC_ABUNDANCE          } from '../subworkflows/local/bgc_abundance'
include { SUMMARY                } from '../modules/local/summary/main'
include { CONSENSUS_SUMMARY      } from '../modules/local/consensus_summary/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

workflow CLYSTERE {
    take:
    ch_samples // channel: [ val(meta), path(genome), path(annotation) ]
    ch_reads   // channel: [ val(meta_reads), path(reads) ]

    main:
    ch_multiqc_files = channel.empty()

    //
    // 1. Genome Preprocessing (GenBank->FASTA bridge, optional dereplication, Pyrodigal/Bakta pre-annotation)
    //
    GENOME_PREPROCESSING(ch_samples)
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_PREPROCESSING.out.multiqc_files)

    //
    // 2. Assembly QC & Taxonomic Context (QUAST, CheckM2, GTDB-Tk)
    //
    GENOME_QC_TAXONOMY(GENOME_PREPROCESSING.out.genome_fastas)
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_QC_TAXONOMY.out.multiqc_files)

    //
    // 3. BGC & Peptide Annotation (antiSMASH 8, GECCO, deepBGC, SanntiS, Macrel + AMPcombi2)
    //
    BGC_ANNOTATION(
        GENOME_PREPROCESSING.out.samples_for_bgc,
        GENOME_PREPROCESSING.out.genome_fastas
    )

    //
    // 4. comBGC Unification & GCF Clustering / Novelty Query (BiG-SCAPE, BiG-SLiCE, BiG-FAM)
    //
    BGC_CLUSTERING(
        BGC_ANNOTATION.out.antismash_gbk_input,
        BGC_ANNOTATION.out.antismash_gbk_results,
        BGC_ANNOTATION.out.antismash_knownclusterblast,
        BGC_ANNOTATION.out.gecco_clusters,
        BGC_ANNOTATION.out.gecco_gbks,
        BGC_ANNOTATION.out.deepbgc_bgc_gbk,
        BGC_ANNOTATION.out.deepbgc_tsv,
        BGC_ANNOTATION.out.antismash_db,
        BGC_ANNOTATION.out.sanntis_gff,
        GENOME_QC_TAXONOMY.out.gtdbtk_summary
    )

    //
    // 5. Comparative Genomics & Synteny Visualization (clinker & cblaster)
    //
    BGC_COMPARATIVE(
        BGC_CLUSTERING.out.combined_regions_dir,
        BGC_ANNOTATION.out.antismash_gbk_input,
        BGC_CLUSTERING.out.bigscape_db,
        BGC_CLUSTERING.out.bigslice_tsv
    )

    //
    // 6. Metagenomic & Metatranscriptomic BGC/GCF Abundance & Expression
    //
    BGC_ABUNDANCE(
        ch_reads,
        BGC_CLUSTERING.out.combined_regions_dir,
        BGC_CLUSTERING.out.bigscape_dereplicate_dir,
        BGC_COMPARATIVE.out.gcf_membership_tsv
    )
    ch_multiqc_files = ch_multiqc_files.mix(BGC_ABUNDANCE.out.abundance_mqc)

    //
    // 7. Tabulation, RiPP Candidate Extraction & Cross-Tool Consensus Summary
    //
    if (params.run_tabulation) {
        ch_all_jsons = BGC_ANNOTATION.out.antismash_json
            .map { _meta, json -> json }
            .collect()
        SUMMARY(ch_all_jsons)
        ch_multiqc_files = ch_multiqc_files.mix(SUMMARY.out.all_regions_tsv.collect())
        ch_multiqc_files = ch_multiqc_files.mix(SUMMARY.out.region_counts_tsv.collect())

        ch_combgc_tsvs = BGC_CLUSTERING.out.summary_tsv
            .map { _meta, tsv -> tsv }
            .collect()
        ch_gcf_opt = BGC_COMPARATIVE.out.gcf_membership_tsv.ifEmpty([])
        ch_bigfam_opt = BGC_CLUSTERING.out.bigfam_novelty_tsv.map { _meta, tsv -> tsv }.ifEmpty([])
        ch_quast_opt = GENOME_QC_TAXONOMY.out.quast_tsv.map { _meta, tsv -> tsv }.collect().ifEmpty([])
        ch_checkm2_opt = GENOME_PREPROCESSING.out.checkm2_tsv
            .mix(GENOME_QC_TAXONOMY.out.checkm2_tsv)
            .map { _meta, tsv -> tsv }
            .collect()
            .ifEmpty([])
        ch_gtdbtk_opt = GENOME_QC_TAXONOMY.out.gtdbtk_summary.map { _meta, files -> files }.flatten().collect().ifEmpty([])
        ch_amp_opt = BGC_ANNOTATION.out.amp_summary_tsv.ifEmpty([])

        CONSENSUS_SUMMARY(
            ch_combgc_tsvs,
            SUMMARY.out.all_regions_tsv,
            ch_gcf_opt,
            ch_bigfam_opt,
            ch_quast_opt,
            ch_checkm2_opt,
            ch_gtdbtk_opt,
            ch_amp_opt,
            params.consensus_min_overlap
        )
        ch_multiqc_files = ch_multiqc_files.mix(CONSENSUS_SUMMARY.out.consensus_support_mqc)
        ch_multiqc_files = ch_multiqc_files.mix(CONSENSUS_SUMMARY.out.genome_overview_mqc)
        ch_multiqc_files = ch_multiqc_files.mix(CONSENSUS_SUMMARY.out.amp_summary_mqc)
    }

    //
    // 8. Software Version Collection & MultiQC Reporting
    //
    ch_collated_versions = softwareVersionsToYAML(channel.topic('versions'))
    ch_versions_file = ch_collated_versions
        .collectFile(name: 'versions.yml', newLine: true)
    ch_multiqc_files = ch_multiqc_files.mix(ch_versions_file)

    ch_multiqc_config = params.multiqc_config ? file(params.multiqc_config) : []
    ch_multiqc_logo = params.multiqc_logo ? file(params.multiqc_logo) : []
    ch_multiqc_input = ch_multiqc_files
        .collect()
        .map { files -> [[id: 'multiqc_report'], files, ch_multiqc_config, ch_multiqc_logo, [], []] }

    MULTIQC(ch_multiqc_input)

    emit:
    multiqc_report = MULTIQC.out.report
    versions       = ch_versions_file
}
