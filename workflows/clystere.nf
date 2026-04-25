/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    antiSMASH + BiG-SCAPE / BiG-SLiCE workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {
  ANTISMASH
} from '../modules/local/antismash/main'
include {
  ANTISMASH_DOWNLOAD_DATABASES
} from '../modules/local/antismash_download_databases/main'
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
  return Channel
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
    ?.findAll {
    it.isDirectory()
  }
    ?.collect {
    file("${it}/Pfam-A.hmm")
  }
    ?.find {
    it.exists()
  }
  if (versioned) {
    return versioned
  }
  error "Could not find Pfam-A.hmm under ${pfam_dir}. Provide it via --bigscape_pfam_path."
}

workflow ANTISMASH_BIGSCAPE {
  //
  // Parse samplesheet
  //
  if (!params.input) {
    error "Please provide a samplesheet with --input."
  }
  if (params.bigscape_run && params.bigslice_run) {
    error "Options --bigscape_run and --bigslice_run are mutually exclusive. Enable only one."
  }
  ch_samples = validateAndParseSamplesheet(params.input)
  //
  // Resolve antiSMASH databases
  // ─────────────────────────────────────────────────────────────────────────
  // 1. --antismash_db provided AND directory is non-empty  → use it as-is
  // 2. --antismash_db provided but missing / empty         → download there
  // 3. --antismash_db not provided                         → download to $TMPDIR/antismash_db
  // ─────────────────────────────────────────────────────────────────────────
  ch_antismash_db = Channel.empty()
  def db_path = params.antismash_db ? file(params.antismash_db): null
  def pfam_dir = db_path ? file("${db_path}/pfam"): null
  def db_exists = pfam_dir && pfam_dir.isDirectory()
  if (db_exists) {
    log.info "Using existing antiSMASH databases at: ${db_path}"
    ch_antismash_db = Channel.value(db_path)
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
    ANTISMASH_DOWNLOAD_DATABASES(Channel.value(download_dest))
    ch_antismash_db = ANTISMASH_DOWNLOAD_DATABASES.out.databases
      .map {
      file(it)
    }
      .first()
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
  // Summary tables (all samples collected into one task each)
  //
  if (params.run_tabulation) {
    ch_all_dirs = ch_antismash_dirs.map {
      it[1]
    }.collect()
    TABULATE_REGIONS(ch_all_dirs)
    COUNT_REGIONS(ch_all_dirs)
  }
  //
  // BiG-SCAPE
  //
  if (params.bigscape_run) {
    // Collect all per-genome antiSMASH output directories
    ch_bgc_input = ch_antismash_dirs.map {
      it[1]
    }.collect()
    ch_pfam = params.bigscape_pfam_path
 ? Channel.value(file(params.bigscape_pfam_path, checkIfExists: true))
: ch_antismash_db.map {
      db -> findPfamHmm(db)
    }
    BIGSCAPE(ch_bgc_input, ch_pfam)
  }
  //
  // BiG-SLiCE
  //
  if (params.bigslice_run) {
    ch_bgc_input = ch_antismash_dirs.map {
      it[1]
    }.collect()
    BIGSLICE(ch_bgc_input)
  }
}
