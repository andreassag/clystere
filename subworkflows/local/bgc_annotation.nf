/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs antiSMASH 8, GECCO, deepBGC, SanntiS, and optional Macrel + AMPcombi2 on input
    genome samples.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANTISMASH_ANTISMASH                  } from '../../modules/nf-core/antismash/antismash/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES } from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { GECCO_RUN                            } from '../../modules/nf-core/gecco/run/main'
include { GECCO_CONVERT                        } from '../../modules/nf-core/gecco/convert/main'
include { DEEPBGC_PIPELINE                     } from '../../modules/nf-core/deepbgc/pipeline/main'
include { DEEPBGC_DOWNLOAD                     } from '../../modules/nf-core/deepbgc/download/main'
include { SANNTIS                              } from '../../modules/local/sanntis/main'
include { AMP_DISCOVERY                        } from './amp_discovery'

workflow BGC_ANNOTATION {
    take:
    ch_samples       // channel: [ val(meta), path(genome), path(annotation) ]
    ch_genome_fastas // channel: [ val(meta), path(fasta) ] (populated when --amp_run is active)

    main:
    // ── Resolve antiSMASH database (hardened sentinel check) ──
    ch_antismash_db = channel.empty()
    def db_path = params.antismash_db ? file(params.antismash_db) : file("${params.clystere_db_dir}/antismash_db")
    def db_exists = file("${db_path}/nrps_pks/transATor").isDirectory() &&
                    (file("${db_path}/pfam").listFiles()?.any { entry -> entry.isDirectory() || entry.name.endsWith('.hmm') } ?: false)
    if (db_exists) {
        log.info "Using existing antiSMASH databases at: ${db_path}"
        ch_antismash_db = channel.value(db_path)
    } else {
        if (params.antismash_db) {
            log.info "antiSMASH database not found at ${db_path} — downloading there."
        } else {
            log.info "No --antismash_db provided — downloading to ${db_path}"
        }
        ANTISMASH_ANTISMASHDOWNLOADDATABASES()
        ch_antismash_db = ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.database.first()
    }

    // ── Run antiSMASH 8 ──
    ch_antismash_in = ch_samples.multiMap { meta, genome, annotation ->
        sequence: [meta, genome]
        gff: annotation ?: []
    }
    ANTISMASH_ANTISMASH(ch_antismash_in.sequence, ch_antismash_db, ch_antismash_in.gff)

    // ── Run GECCO (optional) ──
    ch_gecco_clusters = channel.empty()
    ch_gecco_gbks = channel.empty()
    if (params.gecco_run) {
        ch_gecco_input = ch_samples.map { meta, genome, _annotation -> [meta, genome, []] }
        GECCO_RUN(ch_gecco_input, [])
        ch_gecco_clusters = GECCO_RUN.out.clusters

        if (params.bigscape_run || params.bigslice_run || params.bigfam_db) {
            ch_gecco_convert_input = GECCO_RUN.out.clusters
                .join(GECCO_RUN.out.gbk)
                .map { meta, clusters, gbk -> [meta, clusters, gbk] }
            GECCO_CONVERT(ch_gecco_convert_input, 'gbk', 'bigslice')
            ch_gecco_gbks = GECCO_CONVERT.out.bigslice
        } else {
            ch_gecco_gbks = GECCO_RUN.out.gbk
        }
    }

    // ── Run deepBGC (optional) ──
    ch_deepbgc_bgc_gbk = channel.empty()
    ch_deepbgc_tsv = channel.empty()
    if (params.deepbgc_run) {
        ch_deepbgc_data = channel.empty()
        def data_path = params.deepbgc_data_dir ? file(params.deepbgc_data_dir) : file("${params.clystere_db_dir}/deepbgc_db")
        def data_exists = file("${data_path}/0.1.0/detector/deepbgc.pkl").exists() ||
                          file("${data_path}/detector/deepbgc.pkl").exists()

        if (data_exists) {
            log.info "Using existing deepBGC data at: ${data_path}"
            ch_deepbgc_data = channel.value(data_path)
        } else {
            if (params.deepbgc_data_dir) {
                log.info "deepBGC data not found at ${data_path} — downloading there."
            } else {
                log.info "No --deepbgc_data_dir provided — downloading deepBGC data to ${data_path}"
            }
            DEEPBGC_DOWNLOAD()
            ch_deepbgc_data = DEEPBGC_DOWNLOAD.out.db.first()
        }

        DEEPBGC_PIPELINE(ch_samples.map { meta, genome, _annotation -> [meta, genome] }, ch_deepbgc_data)
        ch_deepbgc_bgc_gbk = DEEPBGC_PIPELINE.out.bgc_gbk
        ch_deepbgc_tsv = DEEPBGC_PIPELINE.out.bgc_tsv
    }

    // ── Run SanntiS as 4th BGC predictor (optional) ──
    ch_sanntis_gff = channel.empty()
    if (params.sanntis_run) {
        def interpro_dir = params.sanntis_interpro_dir ? file(params.sanntis_interpro_dir, checkIfExists: true) : null
        ch_sanntis_in = ANTISMASH_ANTISMASH.out.gbk_input.map { meta, gbk ->
            def ip_file = []
            if (interpro_dir != null && interpro_dir.isDirectory()) {
                def match = interpro_dir.listFiles()?.find { entry ->
                    entry.isFile() &&
                    entry.name.startsWith(meta.id) &&
                    (entry.name.endsWith('.gff3') || entry.name.endsWith('.gff') || entry.name.endsWith('.tsv'))
                }
                if (match != null) {
                    ip_file = match
                }
            }
            [meta, gbk, ip_file]
        }
        SANNTIS(ch_sanntis_in)
        ch_sanntis_gff = SANNTIS.out.gff
    }

    // ── Run Macrel + AMPcombi2 antimicrobial peptide discovery (optional) ──
    ch_amp_summary_tsv = channel.empty()
    ch_amp_cluster_tsv = channel.empty()
    if (params.amp_run) {
        AMP_DISCOVERY(ch_genome_fastas, ANTISMASH_ANTISMASH.out.gbk_input)
        ch_amp_summary_tsv = AMP_DISCOVERY.out.amp_summary_tsv
        ch_amp_cluster_tsv = AMP_DISCOVERY.out.amp_cluster_tsv
    }

    emit:
    antismash_gbk_input         = ANTISMASH_ANTISMASH.out.gbk_input
    antismash_gbk_results       = ANTISMASH_ANTISMASH.out.gbk_results
    antismash_knownclusterblast = ANTISMASH_ANTISMASH.out.knownclusterblast_dir
    antismash_json              = ANTISMASH_ANTISMASH.out.json_results
    antismash_db                = ch_antismash_db
    gecco_clusters              = ch_gecco_clusters
    gecco_gbks                  = ch_gecco_gbks
    deepbgc_bgc_gbk             = ch_deepbgc_bgc_gbk
    deepbgc_tsv                 = ch_deepbgc_tsv
    sanntis_gff                 = ch_sanntis_gff
    amp_summary_tsv             = ch_amp_summary_tsv
    amp_cluster_tsv             = ch_amp_cluster_tsv
}
