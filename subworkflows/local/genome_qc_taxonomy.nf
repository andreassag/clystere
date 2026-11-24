/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: GENOME_QC_TAXONOMY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Optional assembly QC and taxonomic classification branches:
    - QUAST (--quast_run)
    - CheckM2 (--checkm2_run) with auto database download to clystere_db_dir
    - GTDB-Tk (--gtdbtk_run) requiring --gtdbtk_db (or optional --gtdbtk_download_db)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUAST             } from '../../modules/nf-core/quast/main'
include { CHECKM2_QC        } from './checkm2_qc'
include { GTDBTK_CLASSIFYWF } from '../../modules/nf-core/gtdbtk/classifywf/main'
include { GTDBTK_DOWNLOADDB } from '../../modules/local/gtdbtk_downloaddb/main'

def resolveGtdbtkDb(pathStr) {
    def p = file(pathStr)
    if (!p.isDirectory()) {
        return null
    }
    if (file("${p}/metadata/metadata.txt").exists()) {
        return p
    }
    def nested = p.listFiles()?.find { entry -> entry.isDirectory() && file("${entry}/metadata/metadata.txt").exists() }
    return nested ?: null
}

workflow GENOME_QC_TAXONOMY {
    take:
    ch_genome_fastas // channel: [ val(meta), path(fasta) ]

    main:
    ch_quast_tsv = channel.empty()
    ch_quast_results = channel.empty()
    ch_checkm2_tsv = channel.empty()
    ch_checkm2_output = channel.empty()
    ch_gtdbtk_summary = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_collected_fastas = ch_genome_fastas
        .map { _meta, fasta -> fasta }
        .collect()

    //
    // 1. QUAST assembly contiguity QC (optional)
    //
    if (params.quast_run) {
        ch_quast_in = ch_collected_fastas.map { files -> [[id: 'quast'], files] }
        QUAST(ch_quast_in, [[], []], [[], []])
        ch_quast_tsv = QUAST.out.tsv
        ch_quast_results = QUAST.out.results
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.map { _meta, dir -> dir })
    }

    //
    // 2. CheckM2 genome completeness & contamination QC (skipped here if already run pre-dereplication)
    //
    if (params.checkm2_run && !params.dereplicate_genomes) {
        CHECKM2_QC(ch_genome_fastas)
        ch_checkm2_tsv = CHECKM2_QC.out.checkm2_tsv
        ch_checkm2_output = CHECKM2_QC.out.checkm2_output
        ch_multiqc_files = ch_multiqc_files.mix(CHECKM2_QC.out.multiqc_files)
    }

    //
    // 3. GTDB-Tk taxonomic classification (optional)
    //
    if (params.gtdbtk_run) {
        def gtdb_base = params.gtdbtk_db ?: "${params.clystere_db_dir}/gtdbtk_db"
        def resolved_gtdb = resolveGtdbtkDb(gtdb_base)
        if (resolved_gtdb != null) {
            log.info "Using existing GTDB-Tk database at: ${resolved_gtdb}"
            ch_gtdbtk_db = channel.value(['gtdbtk_db', resolved_gtdb])
        } else if (params.gtdbtk_download_db) {
            log.info "GTDB-Tk database not found at ${gtdb_base} — downloading from ${params.gtdbtk_db_url}"
            GTDBTK_DOWNLOADDB(params.gtdbtk_db_url)
            ch_gtdbtk_db = GTDBTK_DOWNLOADDB.out.database.first()
        } else {
            error "Option --gtdbtk_run requires a valid GTDB-Tk database via --gtdbtk_db <path> (containing metadata/metadata.txt), or set --gtdbtk_download_db true."
        }

        ch_gtdb_in = ch_collected_fastas.map { files -> [[id: 'gtdbtk'], files] }
        GTDBTK_CLASSIFYWF(ch_gtdb_in, ch_gtdbtk_db, false)
        ch_gtdbtk_summary = GTDBTK_CLASSIFYWF.out.summary
    }

    emit:
    quast_tsv      = ch_quast_tsv
    quast_results  = ch_quast_results
    checkm2_tsv    = ch_checkm2_tsv
    checkm2_output = ch_checkm2_output
    gtdbtk_summary = ch_gtdbtk_summary
    multiqc_files  = ch_multiqc_files
}
