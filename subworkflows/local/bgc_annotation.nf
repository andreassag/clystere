/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs antiSMASH, GECCO, and DeepBGC annotation tools on input genome samples
    using official nf-core modules.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {
  ANTISMASH_ANTISMASH as ANTISMASH
} from '../../modules/nf-core/antismash/antismash/main'
include {
  ANTISMASH_ANTISMASHDOWNLOADDATABASES as ANTISMASH_DOWNLOAD_DATABASES
} from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include {
  GECCO_RUN as GECCO
} from '../../modules/nf-core/gecco/run/main'
include {
  GECCO_CONVERT
} from '../../modules/nf-core/gecco/convert/main'
include {
  DEEPBGC_PIPELINE as DEEPBGC
} from '../../modules/nf-core/deepbgc/pipeline/main'
include {
  DEEPBGC_DOWNLOAD as DEEPBGC_DOWNLOAD_DATA
} from '../../modules/nf-core/deepbgc/download/main'

workflow BGC_ANNOTATION {
  take:
  ch_samples // channel: [ val(meta), path(genome) ]

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
    ANTISMASH_DOWNLOAD_DATABASES()
    ch_antismash_db = ANTISMASH_DOWNLOAD_DATABASES.out.database
  }

  // ── Run antiSMASH ──
  ANTISMASH(ch_samples, ch_antismash_db, [])
  ch_antismash_dirs = ANTISMASH.out.json_results.map { meta, json -> [meta, json.parent] }

  // ── Run GECCO ──
  ch_gecco_input = ch_samples.map { meta, genome -> [meta, genome, []] }
  GECCO(ch_gecco_input, [])
  ch_gecco_output_dir = GECCO.out.features.map { meta, feat -> [meta, feat.parent] }

  if (params.bigscape_run || params.bigslice_run) {
    ch_gecco_convert_input = GECCO.out.clusters
      .join(GECCO.out.gbk)
      .map { meta, clusters, gbk -> [meta, clusters, gbk] }
    GECCO_CONVERT(ch_gecco_convert_input, 'gbk', 'bigslice')
  }

  // ── Resolve DeepBGC database ──
  ch_deepbgc_data = channel.empty()
  def data_path = params.deepbgc_data_dir ? file(params.deepbgc_data_dir) : file("${params.clystere_db_dir}/deepbgc_data")
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
    DEEPBGC_DOWNLOAD_DATA()
    ch_deepbgc_data = DEEPBGC_DOWNLOAD_DATA.out.db
  }

  // ── Run DeepBGC ──
  DEEPBGC(ch_samples, ch_deepbgc_data)

  emit:
  antismash_dirs      = ch_antismash_dirs
  antismash_db        = ch_antismash_db
  gecco_output_dir    = ch_gecco_output_dir
  gecco_clusters      = GECCO.out.clusters
  deepbgc_bgc_gbk     = DEEPBGC.out.bgc_gbk
  deepbgc_tsv         = DEEPBGC.out.bgc_tsv
}
