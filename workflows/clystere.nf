/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    antiSMASH + BiG-SCAPE / BiG-SLiCE workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {
  ANTISMASH_ANTISMASH as ANTISMASH
} from '../modules/local/antismash/antismash/main'
include {
  ANTISMASH_ANTISMASHDOWNLOADDATABASES as ANTISMASH_DOWNLOAD_DATABASES
} from '../modules/local/antismash/antismash_download_databases/main'
include {
  GECCO_GECCO as GECCO
} from '../modules/local/gecco/main'
include {
  DEEPBGC_PIPELINE as DEEPBGC
} from '../modules/local/deepbgc/deepbgc/main'
include {
  DEEPBGC_DOWNLOAD as DEEPBGC_DOWNLOAD_DATA
} from '../modules/local/deepbgc/deepbgc_download_data/main'
include {
  COMBGC_FILTER
} from '../modules/local/combgc/main'
include {
  BIGSCAPE
} from '../modules/local/bigscape/main'
include {
  BIGSLICE
} from '../modules/local/bigslice/main'
include {
  TABULATE_REGIONS
} from '../modules/local/tabulate_regions/main'
include {
  COUNT_REGIONS
} from '../modules/local/count_regions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate samplesheet and build input channel
    Expected CSV columns: sample, genome[, annotation]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def validateAndParseSamplesheet(samplesheet) {
  def ssDir = file(samplesheet).parent
  return channel
    .fromPath(samplesheet, checkIfExists: true)
    .splitCsv(header: true, strip: true)
    .map {
    row -> // Required fields
    if (!row.sample) {
      error "Samplesheet row missing 'sample' column: ${row}"
    }
    if (!row.genome) {
      error "Samplesheet row missing 'genome' column: ${row}"
    }
    def meta = [id: row.sample]
    def genomePath = row.genome.startsWith('/') ? row.genome: "${ssDir}/${row.genome}"
    def genome = file(genomePath, checkIfExists: true)
    def annotation = row.annotation ? file(row.annotation.startsWith('/') ? row.annotation: "${ssDir}/${row.annotation}", checkIfExists: true): []
    return [meta, genome, annotation]
  }
}

// Locate Pfam-A.hmm within an antiSMASH database root.
// antiSMASH stores it at <db>/pfam/Pfam-A.hmm or in a versioned subdirectory
// <db>/pfam/<version>/Pfam-A.hmm (e.g. pfam/35.0/Pfam-A.hmm).
def findPfamHmm(db) {
  def pfam_dir = file("${db}/pfam")
  if (!pfam_dir.isDirectory()) {
    error "No pfam/ directory found in ${db}. Provide it via --bigscape_pfam_path."
  }
  // Flat: <db>/pfam/Pfam-A.hmm
  def flat = file("${pfam_dir}/Pfam-A.hmm")
  if (flat.exists()) {
    return flat
  }
  // Versioned: <db>/pfam/<version>/Pfam-A.hmm
  def versioned = pfam_dir.listFiles()
    ?.findAll { entry ->
    entry.isDirectory()
  }
    ?.collect { entry ->
    file("${entry}/Pfam-A.hmm")
  }
    ?.find { candidate ->
    candidate.exists()
  }
  if (versioned) {
    return versioned
  }
  error "Could not find Pfam-A.hmm under ${pfam_dir}. Provide it via --bigscape_pfam_path."
}

workflow CLYSTERE {
  //
  // Parse samplesheet
  //
  if (!params.input) {
    error "Please provide a samplesheet with --input."
  }
  if (params.bigscape_run && params.bigslice_run) {
    error "Options --bigscape_run and --bigslice_run are mutually exclusive. Enable only one."
  }
  if ((params.bigscape_run || params.bigslice_run) && (params.disable_gecco || params.disable_deepbgc)) {
    error "BiG-SCAPE/BiG-SLiCE runs require both GECCO and deepBGC to be enabled (not disabled) for comBGC-based unification."
  }
  ch_samples = validateAndParseSamplesheet(params.input)
  //
  // Resolve antiSMASH databases
  // ─────────────────────────────────────────────────────────────────────────
  // 1. --antismash_db provided AND directory is non-empty  → use it as-is
  // 2. --antismash_db provided but missing / empty         → download there
  // 3. --antismash_db not provided                         → download to $TMPDIR/antismash_db
  // ─────────────────────────────────────────────────────────────────────────
  ch_antismash_db = channel.empty()
  def db_path = params.antismash_db ? file(params.antismash_db): null
  def pfam_dir = db_path ? file("${db_path}/pfam"): null
  def db_exists = pfam_dir && pfam_dir.isDirectory()
  if (db_exists) {
    log.info "Using existing antiSMASH databases at: ${db_path}"
    ch_antismash_db = channel.value(db_path)
  } else {
    def download_dest
    if (db_path) {
      log.info "antiSMASH database not found at ${db_path} — downloading there."
      download_dest = db_path.toString()
    } else {
      def tmp = System.getenv('TMPDIR') ?: System.getenv('TMP') ?: '/tmp'
      download_dest = "${tmp}/antismash_db"
      log.info "No --antismash_db provided — downloading to ${download_dest}"
    }
    ANTISMASH_DOWNLOAD_DATABASES(channel.value(download_dest))
    ch_antismash_db = ANTISMASH_DOWNLOAD_DATABASES.out.databases
      .map { dbDir ->
      file(dbDir)
    }
  }
  //
  // Run antiSMASH (one task per genome)
  //
  ch_antismash_input = ch_samples.combine(ch_antismash_db)
    .map {
    meta, genome, annotation, db -> [meta, genome, annotation, db]
  }
  ANTISMASH(ch_antismash_input)
  ch_antismash_dirs = ANTISMASH.out.output_dir
  //
  // Run GECCO (one task per genome)
  //
  ch_gecco_clusters = channel.empty()
  ch_gecco_bigslice_dirs = channel.empty()
  if (!params.disable_gecco) {
    ch_gecco_compat = channel.value(file("${projectDir}/bin/gecco_run_compat.py", checkIfExists: true))
    ch_gecco_input = ch_samples.combine(ch_gecco_compat)
      .map {
      meta, genome, annotation, gecco_compat -> [meta, genome, annotation, gecco_compat]
    }
    GECCO(ch_gecco_input)
    ch_gecco_clusters = GECCO.out.clusters_tsv
    ch_gecco_bigslice_dirs = GECCO.out.bigslice_dir
  }
  //
  // Resolve deepBGC models and run deepBGC (one task per genome)
  //
  ch_deepbgc_bgc_tsv = channel.empty()
  ch_deepbgc_bigslice_dirs = channel.empty()
  if (!params.disable_deepbgc) {
    ch_deepbgc_data = channel.empty()
    def data_path = params.deepbgc_data_dir ? file(params.deepbgc_data_dir): null
    def detector = data_path ? file("${data_path}/detector/deepbgc.pkl"): null
    def data_exists = detector && detector.exists()
    if (data_exists) {
      log.info "Using existing deepBGC data at: ${data_path}"
      ch_deepbgc_data = channel.value(data_path)
    } else {
      def download_dest
      if (data_path) {
        log.info "deepBGC data not found at ${data_path} — downloading there."
        download_dest = data_path.toString()
      } else {
        def tmp = System.getenv('TMPDIR') ?: System.getenv('TMP') ?: '/tmp'
        download_dest = "${tmp}/deepbgc_data"
        log.info "No --deepbgc_data_dir provided — downloading deepBGC data to ${download_dest}"
      }
      DEEPBGC_DOWNLOAD_DATA(channel.value(download_dest))
      ch_deepbgc_data = DEEPBGC_DOWNLOAD_DATA.out.data_dir
        .map { deepbgcDir ->
        file(deepbgcDir)
      }
    }
    ch_deepbgc_input = ch_samples.combine(ch_deepbgc_data)
      .map {
      meta, genome, annotation, data -> [meta, genome, annotation, data]
    }
    DEEPBGC(ch_deepbgc_input)
    ch_deepbgc_bgc_tsv = DEEPBGC.out.bgc_tsv
    ch_deepbgc_bigslice_dirs = DEEPBGC.out.bigslice_dir
  }
  //
  // comBGC unification (runs if GECCO AND deepBGC are enabled, or if clustering is enabled)
  //
  ch_bgc_dirs_for_clustering = ch_antismash_dirs
  if ((!params.disable_gecco && !params.disable_deepbgc) || params.bigscape_run || params.bigslice_run) {
    ch_combgc_input = ch_antismash_dirs
      .join(ch_gecco_bigslice_dirs)
      .join(ch_deepbgc_bigslice_dirs)
      .join(ch_gecco_clusters)
      .join(ch_deepbgc_bgc_tsv)
      .map {
      meta, antismash_dir, gecco_bigslice_dir, deepbgc_bigslice_dir, gecco_clusters_tsv, deepbgc_tsv ->
      [meta, antismash_dir, gecco_bigslice_dir, deepbgc_bigslice_dir, gecco_clusters_tsv, deepbgc_tsv]
    }
    COMBGC_FILTER(ch_combgc_input)
    ch_bgc_dirs_for_clustering = COMBGC_FILTER.out.combined_regions_dir
  }
  //
  // Summary tables (all samples collected into one task each)
  //
  if (params.run_tabulation) {
    ch_all_dirs = ch_antismash_dirs.map {
      entry -> entry[1]
    }.collect()
    TABULATE_REGIONS(ch_all_dirs)
    COUNT_REGIONS(ch_all_dirs)
  }
  if (params.bigscape_run) {
    // Collect all per-genome unified region directories
    ch_bgc_input = ch_bgc_dirs_for_clustering.map {
      entry -> entry[1]
    }.collect()
    ch_pfam = params.bigscape_pfam_path
 ? channel.value(file(params.bigscape_pfam_path, checkIfExists: true))
: ch_antismash_db.map {
      db -> findPfamHmm(db)
    }
    BIGSCAPE(ch_bgc_input, ch_pfam)
  }
  //
  // BiG-SLiCE
  //
  if (params.bigslice_run) {
    ch_bgc_input = ch_bgc_dirs_for_clustering.map {
      entry -> entry[1]
    }.collect()
    BIGSLICE(ch_bgc_input)
  }
}
