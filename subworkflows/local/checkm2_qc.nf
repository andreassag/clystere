/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: CHECKM2_QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared CheckM2 database resolution and completeness/contamination prediction.
    Invoked once per run: either prior to dereplication (in GENOME_PREPROCESSING when
    both --dereplicate_genomes and --checkm2_run are enabled) or in GENOME_QC_TAXONOMY.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECKM2_PREDICT          } from '../../modules/nf-core/checkm2/predict/main'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'

def resolveCheckm2Dmnd(pathStr) {
    def p = file(pathStr)
    if (p.isFile() && p.name.endsWith('.dmnd')) {
        return p
    }
    if (p.isDirectory()) {
        def found = p.listFiles()?.find { entry -> entry.isFile() && entry.name.endsWith('.dmnd') }
        if (found) {
            return found
        }
    }
    return null
}

workflow CHECKM2_QC {
    take:
    ch_genome_fastas // channel: [ val(meta), path(fasta) ]

    main:
    ch_collected_fastas = ch_genome_fastas
        .map { _meta, fasta -> fasta }
        .collect()

    ch_trigger = ch_collected_fastas.map { true }.first()

    def checkm2_base = params.checkm2_db ?: "${params.clystere_db_dir}/checkm2_db"
    def dmnd_file = resolveCheckm2Dmnd(checkm2_base)
    if (dmnd_file != null) {
        log.info "Using existing CheckM2 database at: ${dmnd_file}"
        ch_checkm2_db = channel.value([[id: 'checkm2_db'], dmnd_file])
    } else {
        log.info "CheckM2 database not found at ${checkm2_base} — downloading automatically."
        ch_zenodo_id = ch_trigger.map { params.checkm2_db_zenodo_id }
        CHECKM2_DATABASEDOWNLOAD(ch_zenodo_id)
        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database.first()
    }

    ch_checkm2_in = ch_collected_fastas.map { files -> [[id: 'checkm2'], files] }
    CHECKM2_PREDICT(ch_checkm2_in, ch_checkm2_db)

    emit:
    checkm2_tsv    = CHECKM2_PREDICT.out.checkm2_tsv
    checkm2_output = CHECKM2_PREDICT.out.checkm2_output
    multiqc_files  = CHECKM2_PREDICT.out.checkm2_tsv.map { _meta, tsv -> tsv }
}
