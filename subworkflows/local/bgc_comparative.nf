/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_COMPARATIVE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Tier 3 comparative genomics & visualization:
    - GROUP_GCF_GBKS: Organizes multi-member GCFs from BiG-SCAPE and/or BiG-SLiCE and
      extracts representative/orphan queries plus gcf_membership_summary.tsv.
    - CLINKER (--clinker_run): Generates interactive HTML synteny plots per GCF.
    - CBLASTER (--cblaster_run): Offline cohort-wide (default) or external/remote
      gene-neighborhood homology searches.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GROUP_GCF_GBKS  } from '../../modules/local/group_gcf_gbks/main'
include { CLINKER         } from '../../modules/local/clinker/main'
include { CBLASTER_MAKEDB } from '../../modules/local/cblaster/makedb/main'
include { CBLASTER_SEARCH } from '../../modules/local/cblaster/search/main'

workflow BGC_COMPARATIVE {
    take:
    ch_combined_regions_dir // channel: [ val(meta), path(regions_dir) ]
    ch_antismash_gbk_input  // channel: [ val(meta), path(gbk) ]
    ch_bigscape_db          // channel: path(bigscape.db) (optional/empty)
    ch_bigslice_tsv         // channel: [ val(meta), path(tsv_export) ] (optional/empty)

    main:
    if (params.cblaster_run && !(params.cblaster_mode in ['local', 'remote'])) {
        error "Invalid --cblaster_mode '${params.cblaster_mode}'. Must be one of: local, remote."
    }

    ch_gcf_membership_tsv = channel.empty()
    ch_clinker_html = channel.empty()
    ch_clinker_alignments = channel.empty()
    ch_cblaster_summary = channel.empty()

    def need_grouping = params.clinker_run || params.cblaster_run || params.bigscape_run || params.bigslice_run
    if (need_grouping) {
        ch_all_bgc_dirs = ch_combined_regions_dir
            .map { _meta, dir -> dir }
            .collect()
        ch_bs_db_opt = ch_bigscape_db.ifEmpty([])
        ch_sl_tsv_opt = ch_bigslice_tsv.map { _meta, dir -> dir }.ifEmpty([])

        GROUP_GCF_GBKS(
            ch_all_bgc_dirs,
            ch_bs_db_opt,
            ch_sl_tsv_opt,
            params.clinker_max_bgcs,
            params.cblaster_max_queries
        )
        ch_gcf_membership_tsv = GROUP_GCF_GBKS.out.membership_tsv

        if (params.clinker_run) {
            ch_gcf_dirs = GROUP_GCF_GBKS.out.gcf_dirs
                .flatten()
                .map { dir -> [[id: dir.name], dir] }
            CLINKER(ch_gcf_dirs)
            ch_clinker_html = CLINKER.out.html
            ch_clinker_alignments = CLINKER.out.alignments
        }

        if (params.cblaster_run) {
            ch_db_dmnd = channel.value([])
            ch_db_sqlite = channel.value([])
            if (params.cblaster_mode == 'local') {
                if (params.cblaster_db) {
                    def raw_db = params.cblaster_db.toString()
                    def base_db = raw_db.endsWith('.dmnd') ? raw_db.replaceFirst(/\.dmnd$/, '') : raw_db
                    ch_db_dmnd = channel.value(file("${base_db}.dmnd", checkIfExists: true))
                    ch_db_sqlite = channel.value(file("${base_db}.sqlite3", checkIfExists: true))
                } else {
                    ch_cohort_gbks = ch_antismash_gbk_input
                        .map { _meta, gbk -> gbk }
                        .collect()
                    CBLASTER_MAKEDB(ch_cohort_gbks)
                    ch_db_dmnd = CBLASTER_MAKEDB.out.dmnd.first()
                    ch_db_sqlite = CBLASTER_MAKEDB.out.sqlite.first()
                }
            }

            ch_queries = GROUP_GCF_GBKS.out.cblaster_queries
                .flatten()
                .map { gbk -> [[id: gbk.baseName], gbk] }
            CBLASTER_SEARCH(ch_queries, ch_db_dmnd, ch_db_sqlite, params.cblaster_mode)
            ch_cblaster_summary = CBLASTER_SEARCH.out.summary_tsv
        }
    }

    emit:
    gcf_membership_tsv = ch_gcf_membership_tsv
    clinker_html       = ch_clinker_html
    clinker_alignments = ch_clinker_alignments
    cblaster_summary   = ch_cblaster_summary
}
