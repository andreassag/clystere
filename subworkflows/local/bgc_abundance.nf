/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_ABUNDANCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Tier 4 metagenomic & metatranscriptomic BGC/GCF read mapping and abundance/expression
    quantification using official nf-core modules (bowtie2/minimap2, samtools, coverm).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_BGC_REFERENCE   } from '../../modules/local/prepare_bgc_reference/main'
include { BOWTIE2_BUILD           } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN           } from '../../modules/nf-core/bowtie2/align/main'
include { MINIMAP2_INDEX          } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN          } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort/main'
include { COVERM_CONTIG           } from '../../modules/nf-core/coverm/contig/main'
include { AGGREGATE_BGC_ABUNDANCE } from '../../modules/local/aggregate_bgc_abundance/main'

workflow BGC_ABUNDANCE {
    take:
    ch_reads                    // channel: [ val(meta_reads), path(reads) ]
    ch_combined_regions_dir     // channel: [ val(meta), path(regions_dir) ]
    ch_bigscape_dereplicate_dir // channel: path(derep_dir) (optional/empty)
    ch_gcf_membership_tsv       // channel: path(gcf_membership_summary.tsv) (optional/empty)

    main:
    if (!(params.abundance_mapper in ['bowtie2', 'minimap2'])) {
        error "Invalid --abundance_mapper '${params.abundance_mapper}'. Must be one of: bowtie2, minimap2."
    }

    // Gate entire subworkflow so zero processes launch when no reads are supplied
    ch_reads_gate = ch_reads.map { true }.first()

    ch_bgc_dirs_gated = ch_combined_regions_dir
        .map { _meta, dir -> dir }
        .collect()
        .map { dirs -> [dirs] }
        .combine(ch_reads_gate)
        .map { dirs, _gate -> dirs }

    ch_derep_opt = ch_bigscape_dereplicate_dir.ifEmpty([])
    ch_gcf_opt = ch_gcf_membership_tsv.ifEmpty([])

    PREPARE_BGC_REFERENCE(ch_bgc_dirs_gated, ch_derep_opt, ch_gcf_opt)
    ch_nonempty_ref = PREPARE_BGC_REFERENCE.out.fasta
        .filter { _meta, fasta -> workflow.stubRun || fasta.size() > 0 }

    ch_sorted_bams = channel.empty()
    if (params.abundance_mapper == 'bowtie2') {
        BOWTIE2_BUILD(ch_nonempty_ref)
        BOWTIE2_ALIGN(
            ch_reads,
            BOWTIE2_BUILD.out.index.first(),
            ch_nonempty_ref.first(),
            false,
            false
        )
        SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam, [[], [], []], 'bai')
        ch_sorted_bams = SAMTOOLS_SORT.out.bam
    } else {
        MINIMAP2_INDEX(ch_nonempty_ref)
        MINIMAP2_ALIGN(
            ch_reads,
            MINIMAP2_INDEX.out.index.first(),
            true,
            'bai',
            false,
            false
        )
        ch_sorted_bams = MINIMAP2_ALIGN.out.bam
    }

    COVERM_CONTIG(ch_sorted_bams, [[], []], true, false, false)

    ch_assay_manifest = ch_reads
        .map { meta, _reads -> "${meta.id}\t${meta.assay_type ?: 'metagenome'}" }
        .collectFile(name: 'assay_manifest.tsv', seed: "sample_id\tassay_type\n", newLine: true)

    ch_coverm_tsvs = COVERM_CONTIG.out.coverage
        .map { _meta, tsv -> tsv }
        .collect()

    AGGREGATE_BGC_ABUNDANCE(
        ch_coverm_tsvs,
        PREPARE_BGC_REFERENCE.out.metadata_tsv,
        ch_gcf_opt,
        ch_assay_manifest
    )

    emit:
    bgc_catalog_fasta        = PREPARE_BGC_REFERENCE.out.fasta
    bgc_tpm_tsv              = AGGREGATE_BGC_ABUNDANCE.out.bgc_tpm_tsv
    bgc_rpkm_tsv             = AGGREGATE_BGC_ABUNDANCE.out.bgc_rpkm_tsv
    bgc_covered_fraction_tsv = AGGREGATE_BGC_ABUNDANCE.out.bgc_covered_fraction_tsv
    gcf_tpm_tsv              = AGGREGATE_BGC_ABUNDANCE.out.gcf_tpm_tsv
    gcf_rpkm_tsv             = AGGREGATE_BGC_ABUNDANCE.out.gcf_rpkm_tsv
    activity_summary_tsv     = AGGREGATE_BGC_ABUNDANCE.out.activity_summary_tsv
    abundance_mqc            = AGGREGATE_BGC_ABUNDANCE.out.abundance_mqc
}
