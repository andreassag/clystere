/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_CLUSTERING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unifies BGC predictions via comBGC and clusters GCFs using BiG-SCAPE or BiG-SLiCE.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {
  COMBGC_FILTER
} from '../../modules/local/combgc/main'
include {
  BIGSCAPE
} from '../../modules/local/bigscape/main'
include {
  BIGSLICE_BIGSLICE as BIGSLICE
} from '../../modules/nf-core/bigslice/bigslice/main'

def findPfamHmm(db) {
  def pfam_dir = file("${db}/pfam")
  if (!pfam_dir.isDirectory()) {
    error "No pfam/ directory found in ${db}. Provide it via --bigscape_pfam_path."
  }
  def flat = file("${pfam_dir}/Pfam-A.hmm")
  if (flat.exists()) {
    return flat
  }
  def versioned = pfam_dir.listFiles()
    ?.findAll { entry -> entry.isDirectory() }
    ?.collect { entry -> file("${entry}/Pfam-A.hmm") }
    ?.find { candidate -> candidate.exists() }
  if (versioned) {
    return versioned
  }
  error "Could not find Pfam-A.hmm under ${pfam_dir}. Provide it via --bigscape_pfam_path."
}

workflow BGC_CLUSTERING {
  take:
  ch_antismash_dirs
  ch_gecco_output_dir
  ch_deepbgc_bgc_gbk
  ch_gecco_clusters
  ch_deepbgc_tsv
  ch_antismash_db

  main:
  ch_combgc_input = ch_antismash_dirs
    .join(ch_gecco_output_dir, remainder: true)
    .join(ch_deepbgc_bgc_gbk, remainder: true)
    .join(ch_gecco_clusters, remainder: true)
    .join(ch_deepbgc_tsv, remainder: true)
    .map { item ->
      def meta = item[0]
      def antismash_dir = item.size() > 1 && item[1] ? item[1] : []
      def gecco_dir = item.size() > 2 && item[2] ? item[2] : []
      def deepbgc_bgc_gbk = item.size() > 3 && item[3] ? item[3] : []
      def gecco_clusters_tsv = item.size() > 4 && item[4] ? item[4] : []
      def deepbgc_tsv = item.size() > 5 && item[5] ? item[5] : []
      [meta, antismash_dir, gecco_dir, deepbgc_bgc_gbk, gecco_clusters_tsv, deepbgc_tsv]
    }

  COMBGC_FILTER(ch_combgc_input)
  ch_bgc_dirs_for_clustering = COMBGC_FILTER.out.combined_regions_dir

  if (params.bigscape_run) {
    ch_bgc_input = ch_bgc_dirs_for_clustering.map { entry -> entry[1] }.collect()
    ch_pfam = params.bigscape_pfam_path
      ? channel.value(file(params.bigscape_pfam_path, checkIfExists: true))
      : ch_antismash_db.map { db -> findPfamHmm(db) }
    BIGSCAPE(ch_bgc_input, ch_pfam)
  }

  if (params.bigslice_run) {
    ch_bigslice_input = ch_bgc_dirs_for_clustering
      .collect()
      .map { list -> [[id: 'all_samples'], list.collect { item -> item[1] }] }
    BIGSLICE(ch_bigslice_input, [], true)
  }

  emit:
  combined_regions_dir = ch_bgc_dirs_for_clustering
}
