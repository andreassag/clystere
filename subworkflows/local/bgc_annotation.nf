/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs antiSMASH, GECCO, and deepBGC annotation tools on input genome samples
    using official nf-core modules.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANTISMASH_ANTISMASH                  } from '../../modules/nf-core/antismash/antismash/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES } from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { GECCO_RUN                            } from '../../modules/nf-core/gecco/run/main'
include { GECCO_CONVERT                        } from '../../modules/nf-core/gecco/convert/main'
include { DEEPBGC_PIPELINE                     } from '../../modules/nf-core/deepbgc/pipeline/main'
include { DEEPBGC_DOWNLOAD                     } from '../../modules/nf-core/deepbgc/download/main'

workflow BGC_ANNOTATION {
    take:
    ch_samples // channel: [ val(meta), path(genome), path(annotation) ]

    main:
    // ── Resolve antiSMASH database ──
    ch_antismash_db = channel.empty()
    def db_path = params.antismash_db ? file(params.antismash_db) : file("${params.clystere_db_dir}/antismash_db")
    def db_exists = file("${db_path}/pfam").isDirectory()
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

    // ── Run antiSMASH ──
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

        if (params.bigscape_run || params.bigslice_run) {
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

    emit:
    antismash_gbk_input           = ANTISMASH_ANTISMASH.out.gbk_input
    antismash_gbk_results         = ANTISMASH_ANTISMASH.out.gbk_results
    antismash_knownclusterblast   = ANTISMASH_ANTISMASH.out.knownclusterblast_dir
    antismash_json                = ANTISMASH_ANTISMASH.out.json_results
    antismash_db                  = ch_antismash_db
    gecco_clusters                = ch_gecco_clusters
    gecco_gbks                    = ch_gecco_gbks
    deepbgc_bgc_gbk               = ch_deepbgc_bgc_gbk
    deepbgc_tsv                   = ch_deepbgc_tsv
}
