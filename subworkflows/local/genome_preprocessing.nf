/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: GENOME_PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Handles:
    1. GenBank -> nucleotide FASTA extraction (and FASTA normalization) when any
       FASTA-requiring module is active.
    2. Optional genome dereplication via skani (default) or dRep before heavy annotation.
    3. Automatic gene-calling / pre-annotation of unannotated FASTA inputs using
       Pyrodigal (default) or Bakta, plus GFF3 decompression for antiSMASH 8.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GBK_TO_FASTA                    } from '../../modules/local/gbk_to_fasta/main'
include { GUNZIP_GFF as GUNZIP_USER_GFF   } from '../../modules/local/gunzip_gff/main'
include { GUNZIP_GFF as GUNZIP_PYRO_GFF   } from '../../modules/local/gunzip_gff/main'
include { DEREP_SELECT                    } from '../../modules/local/derep_select/main'
include { CHECKM2_TO_DREP_INFO            } from '../../modules/local/checkm2_to_drep_info/main'
include { CHECKM2_QC                      } from './checkm2_qc'
include { PYRODIGAL                       } from '../../modules/nf-core/pyrodigal/main'
include { BAKTA_BAKTA                     } from '../../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD           } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { SKANI_TRIANGLE                  } from '../../modules/nf-core/skani/triangle/main'
include { DREP_DEREPLICATE                } from '../../modules/nf-core/drep/dereplicate/main'

def isFastaFile(path) {
    return path.name.toLowerCase() ==~ /.*\.(fa|fna|fasta)(\.gz)?$/
}

def isGzippedFile(path) {
    return path.name.toLowerCase().endsWith('.gz')
}

def stemSampleId(String filename) {
    def lower = filename.toLowerCase()
    def matched = ['.fasta.gz', '.fna.gz', '.fa.gz', '.fasta', '.fna', '.fa'].find { suffix -> lower.endsWith(suffix) }
    if (matched) {
        return filename.substring(0, filename.length() - matched.length())
    }
    return filename.replaceFirst(/\.[^.]+$/, '')
}

def resolveBaktaDbDir(basePath) {
    def base = file(basePath)
    return [base, file("${base}/db-light"), file("${base}/db")].find { candidate ->
        file("${candidate}/bakta.db").exists()
    }
}

workflow GENOME_PREPROCESSING {
    take:
    ch_samples // channel: [ val(meta), path(genome), path(annotation) ]

    main:
    if (!(params.annotation_tool in ['pyrodigal', 'bakta', 'none'])) {
        error "Invalid --annotation_tool '${params.annotation_tool}'. Must be one of: pyrodigal, bakta, none."
    }
    if (params.dereplicate_genomes && !(params.dereplicate_tool in ['skani', 'drep'])) {
        error "Invalid --dereplicate_tool '${params.dereplicate_tool}'. Must be one of: skani, drep."
    }

    def need_fastas = params.quast_run || params.checkm2_run || params.gtdbtk_run || params.amp_run || params.dereplicate_genomes

    //
    // 1. Build unified uncompressed <sample_id>.fna channel when FASTA-requiring tools are enabled
    //
    ch_genome_fastas = channel.empty()
    if (need_fastas) {
        ch_fasta_split = ch_samples.branch { meta, genome, _ann ->
            ready_fna: isFastaFile(genome) && !isGzippedFile(genome) && genome.name == "${meta.id}.fna"
            needs_convert: true
        }
        ch_to_convert = ch_fasta_split.needs_convert.multiMap { meta, genome, _ann ->
            genome_in: [meta, genome]
            is_fasta: isFastaFile(genome)
        }
        GBK_TO_FASTA(ch_to_convert.genome_in, ch_to_convert.is_fasta)
        ch_genome_fastas = ch_fasta_split.ready_fna
            .map { meta, genome, _ann -> [meta, genome] }
            .mix(GBK_TO_FASTA.out.fasta)
    }

    //
    // 2. Optional Tier 4 Genome Dereplication (skani or dRep, weighted by CheckM2 when enabled)
    //
    ch_active_samples = ch_samples
    ch_active_fastas = ch_genome_fastas
    ch_derep_clusters = channel.empty()
    ch_checkm2_tsv = channel.empty()
    ch_checkm2_output = channel.empty()
    ch_multiqc_files = channel.empty()

    if (params.dereplicate_genomes) {
        ch_all_fastas_list = ch_genome_fastas
            .map { _meta, fa -> fa }
            .collect()
        ch_all_fastas_tuple = ch_all_fastas_list.map { files -> [[id: 'cohort'], files] }

        ch_checkm2_for_derep = channel.value([])
        ch_drep_genome_info = channel.value([])
        if (params.checkm2_run) {
            CHECKM2_QC(ch_genome_fastas)
            ch_checkm2_tsv = CHECKM2_QC.out.checkm2_tsv
            ch_checkm2_output = CHECKM2_QC.out.checkm2_output
            ch_multiqc_files = ch_multiqc_files.mix(CHECKM2_QC.out.multiqc_files)
            ch_checkm2_for_derep = CHECKM2_QC.out.checkm2_tsv.map { _meta, tsv -> tsv }.first()
            if (params.dereplicate_tool == 'drep') {
                CHECKM2_TO_DREP_INFO(CHECKM2_QC.out.checkm2_tsv, ch_all_fastas_list)
                ch_drep_genome_info = CHECKM2_TO_DREP_INFO.out.genome_info_csv.first()
            }
        }

        ch_winning_ids = channel.empty()
        if (params.dereplicate_tool == 'skani') {
            SKANI_TRIANGLE(ch_all_fastas_tuple)
            DEREP_SELECT(SKANI_TRIANGLE.out.triangle, ch_all_fastas_list, params.dereplicate_ani, ch_checkm2_for_derep)
            ch_derep_clusters = DEREP_SELECT.out.clusters_tsv
            ch_winning_ids = DEREP_SELECT.out.ids
                .map { ids_file ->
                    ids_file.readLines().collect { line -> line.trim() }.findAll { line -> line } as Set
                }
        } else {
            DREP_DEREPLICATE(ch_all_fastas_tuple, [[], []], ch_drep_genome_info)
            ch_derep_clusters = DREP_DEREPLICATE.out.summary_tables
                .map { _meta, tables -> (tables instanceof List) ? tables : [tables] }
                .flatten()
                .filter { tbl -> tbl.name == 'Cdb.csv' || tbl.name == 'Wdb.csv' }
                .first()
            ch_winning_ids = DREP_DEREPLICATE.out.fastas
                .map { _meta, files ->
                    def list = (files instanceof List) ? files : [files]
                    list.collect { f -> stemSampleId(f.name) } as Set
                }
        }

        ch_active_samples = ch_samples
            .combine(ch_winning_ids)
            .filter { meta, _genome, _ann, winners ->
                (meta.id in winners) || (workflow.stubRun && winners.every { w -> w.startsWith('test') })
            }
            .map { meta, genome, ann, _winners -> [meta, genome, ann] }

        ch_active_fastas = ch_genome_fastas
            .combine(ch_winning_ids)
            .filter { meta, _fa, winners ->
                (meta.id in winners) || (workflow.stubRun && winners.every { w -> w.startsWith('test') })
            }
            .map { meta, fa, _winners -> [meta, fa] }
    }

    //
    // 3. Pre-annotate unannotated FASTA genomes & decompress gzipped GFF3 annotations
    //
    ch_split = ch_active_samples.branch { _meta, genome, annotation ->
        unannotated_fasta: isFastaFile(genome) && !annotation
        annotated_gz_gff: annotation && isGzippedFile(annotation)
        annotated_ready: true
    }

    GUNZIP_USER_GFF(ch_split.annotated_gz_gff.map { meta, _genome, ann -> [meta, ann] })
    ch_annotated_unpacked = ch_split.annotated_gz_gff
        .map { meta, genome, _ann -> [meta, genome] }
        .join(GUNZIP_USER_GFF.out.gff)
        .map { meta, genome, gff3 -> [meta, genome, gff3] }

    ch_preannotated = channel.empty()
    if (params.annotation_tool == 'pyrodigal') {
        PYRODIGAL(ch_split.unannotated_fasta.map { meta, genome, _ann -> [meta, genome] }, 'gff')
        GUNZIP_PYRO_GFF(PYRODIGAL.out.annotations)
        ch_preannotated = ch_split.unannotated_fasta
            .map { meta, genome, _ann -> [meta, genome] }
            .join(GUNZIP_PYRO_GFF.out.gff)
            .map { meta, genome, gff3 -> [meta, genome, gff3] }
    } else if (params.annotation_tool == 'bakta') {
        ch_bakta_trigger = ch_split.unannotated_fasta.map { true }.first()
        def bakta_base = params.bakta_db ?: "${params.clystere_db_dir}/bakta_db"
        def resolved_bakta = resolveBaktaDbDir(bakta_base)
        if (resolved_bakta != null) {
            log.info "Using existing Bakta database at: ${resolved_bakta}"
            ch_bakta_db = channel.value(resolved_bakta)
        } else {
            log.info "No existing Bakta database at ${bakta_base} — downloading (${params.bakta_db_type}) if unannotated FASTAs are present."
            BAKTA_BAKTADBDOWNLOAD(ch_bakta_trigger)
            ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db.first()
        }
        BAKTA_BAKTA(
            ch_split.unannotated_fasta.map { meta, genome, _ann -> [meta, genome] },
            ch_bakta_db,
            [],
            [],
            [],
            []
        )
        ch_preannotated = BAKTA_BAKTA.out.gbff.map { meta, gbff -> [meta, gbff, []] }
    } else {
        ch_preannotated = ch_split.unannotated_fasta
    }

    ch_samples_for_bgc = ch_split.annotated_ready
        .mix(ch_annotated_unpacked)
        .mix(ch_preannotated)

    emit:
    samples_for_bgc        = ch_samples_for_bgc
    genome_fastas          = ch_active_fastas
    dereplication_clusters = ch_derep_clusters
    checkm2_tsv            = ch_checkm2_tsv
    checkm2_output         = ch_checkm2_output
    multiqc_files          = ch_multiqc_files
}
