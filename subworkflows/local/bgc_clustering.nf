/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_CLUSTERING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unifies BGC predictions via comBGC and clusters GCFs using BiG-SCAPE or BiG-SLiCE.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COMBGC_FILTER       } from '../../modules/local/combgc/main'
include { BIGSCAPE            } from '../../modules/local/bigscape/main'
include { BIGSLICE_BIGSLICE   } from '../../modules/nf-core/bigslice/bigslice/main'
include { BIGSLICE_DOWNLOADDB } from '../../modules/nf-core/bigslice/downloaddb/main'

def findPfamHmm(db) {
    def pfam_dir = file("${db}/pfam")
    if (!pfam_dir.isDirectory()) {
        error "No pfam/ directory found in ${db}. Provide Pfam-A.hmm via --bigscape_pfam_path."
    }
    def flat = ['Pfam-A.hmm', 'Pfam-A.hmm.gz']
        .collect { name -> file("${pfam_dir}/${name}") }
        .find { candidate -> candidate.exists() }
    if (flat) {
        return flat
    }
    def versioned = pfam_dir.listFiles()
        ?.findAll { entry -> entry.isDirectory() }
        ?.collectMany { entry -> [file("${entry}/Pfam-A.hmm"), file("${entry}/Pfam-A.hmm.gz")] }
        ?.find { candidate -> candidate.exists() }
    if (versioned) {
        return versioned
    }
    error "Could not find Pfam-A.hmm under ${pfam_dir}. Provide it via --bigscape_pfam_path."
}

workflow BGC_CLUSTERING {
    take:
    ch_antismash_gbk_input
    ch_antismash_gbk_results
    ch_antismash_knownclusterblast
    ch_gecco_clusters
    ch_gecco_gbks
    ch_deepbgc_bgc_gbk
    ch_deepbgc_tsv
    ch_antismash_db

    main:
    ch_combgc_input = ch_antismash_gbk_input
        .join(ch_antismash_gbk_results, remainder: true)
        .join(ch_antismash_knownclusterblast, remainder: true)
        .join(ch_gecco_clusters, remainder: true)
        .join(ch_gecco_gbks, remainder: true)
        .join(ch_deepbgc_bgc_gbk, remainder: true)
        .join(ch_deepbgc_tsv, remainder: true)
        .map { item ->
            def meta = item[0]
            def as_gbk = item[1] ?: []
            def as_regions = item[2] ?: []
            def as_kcb = item[3] ?: []
            def gecco_tsv = item[4] ?: []
            def gecco_gbk = item[5] ?: []
            def deepbgc_gbk = item[6] ?: []
            def deepbgc_tsv = item[7] ?: []
            [meta, as_gbk, as_regions, as_kcb, gecco_tsv, gecco_gbk, deepbgc_gbk, deepbgc_tsv]
        }

    COMBGC_FILTER(ch_combgc_input)
    ch_bgc_dirs_for_clustering = COMBGC_FILTER.out.combined_regions_dir

    if (params.bigscape_run) {
        ch_bgc_input = ch_bgc_dirs_for_clustering
            .filter { _meta, dir ->
                dir.listFiles()?.any { f -> f.name.endsWith('.gbk') && (workflow.stubRun || f.size() > 0) }
            }
            .map { _meta, dir -> dir }
            .collect()
        ch_pfam = params.bigscape_pfam_path
            ? channel.value(file(params.bigscape_pfam_path, checkIfExists: true))
            : ch_antismash_db.map { db -> findPfamHmm(db) }.first()
        BIGSCAPE(ch_bgc_input, ch_pfam, params.bigscape_dereplicate, params.bigscape_zip_output)
    }

    if (params.bigslice_run) {
        ch_bigslice_db = channel.empty()
        def bs_db_path = params.bigslice_db ? file(params.bigslice_db) : file("${params.clystere_db_dir}/bigslice-models")
        def bs_db_exists = file("${bs_db_path}/biosyn_pfam/Biosyn_pfams.hmm").exists() ||
                           file("${bs_db_path}/biosynthetic_pfams/biopfam.md5sum").exists() ||
                           file("${bs_db_path}/sub_pfams").isDirectory()
        if (bs_db_exists) {
            log.info "Using existing BiG-SLiCE HMM database at: ${bs_db_path}"
            ch_bigslice_db = channel.value(bs_db_path)
        } else {
            if (params.bigslice_db) {
                log.info "BiG-SLiCE database not found at ${bs_db_path} — downloading there."
            } else {
                log.info "No --bigslice_db provided — downloading BiG-SLiCE HMM database to ${bs_db_path}"
            }
            BIGSLICE_DOWNLOADDB([id: 'bigslice_models'])
            ch_bigslice_db = BIGSLICE_DOWNLOADDB.out.db.map { _meta, db -> db }.first()
        }

        ch_bigslice_input = COMBGC_FILTER.out.combined_region_gbks
            .map { _meta, gbks -> (gbks instanceof List) ? gbks : [gbks] }
            .flatten()
            .filter { gbk -> workflow.stubRun || gbk.size() > 0 }
            .collect()
            .map { files -> [[id: 'all_samples'], files] }
        BIGSLICE_BIGSLICE(ch_bigslice_input, ch_bigslice_db, true)
    }

    emit:
    combined_regions_dir = ch_bgc_dirs_for_clustering
    summary_tsv          = COMBGC_FILTER.out.summary_tsv
}
