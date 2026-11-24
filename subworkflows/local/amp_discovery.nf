/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: AMP_DISCOVERY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Antimicrobial peptide (AMP) discovery using Macrel and AMPcombi2, with single-pass
    AMP reference database download/indexing via PREPARE_AMP_DB to prevent parallel
    write races and duplicate downloads across samples.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MACREL_CONTIGS        } from '../../modules/nf-core/macrel/contigs/main'
include { PREPARE_AMP_DB        } from '../../modules/local/prepare_amp_db/main'
include { PREPARE_AMPCOMBI_INPUT } from '../../modules/local/prepare_ampcombi_input/main'
include { AMPCOMBI2_PARSETABLES } from '../../modules/nf-core/ampcombi2/parsetables/main'
include { AMPCOMBI2_COMPLETE    } from '../../modules/nf-core/ampcombi2/complete/main'
include { AMPCOMBI2_CLUSTER     } from '../../modules/nf-core/ampcombi2/cluster/main'

workflow AMP_DISCOVERY {
    take:
    ch_genome_fastas      // channel: [ val(meta), path(fasta) ]
    ch_antismash_gbk_in   // channel: [ val(meta), path(gbk) ] (optional/empty)

    main:
    MACREL_CONTIGS(ch_genome_fastas)

    ch_macrel_joined = MACREL_CONTIGS.out.amp_prediction
        .join(MACREL_CONTIGS.out.all_orfs)
        .join(ch_genome_fastas)
        .join(ch_antismash_gbk_in, remainder: true)
        .filter { _meta, pred, orfs, fasta, _as_gbk -> pred != null && orfs != null && fasta != null }
        .map { meta, pred, orfs, fasta, as_gbk -> [meta, pred, orfs, fasta, as_gbk ?: []] }

    PREPARE_AMPCOMBI_INPUT(ch_macrel_joined)

    ch_ampcombi_in = PREPARE_AMPCOMBI_INPUT.out.prediction
        .join(PREPARE_AMPCOMBI_INPUT.out.faa)
        .join(PREPARE_AMPCOMBI_INPUT.out.gbk)
        .multiMap { meta, pred, faa, gbk ->
            amp_input: [meta, [pred]]
            faa: faa
            gbk: gbk
        }

    ch_amp_db_dir = channel.empty()
    ch_amp_trigger = ch_ampcombi_in.amp_input.map { true }.first()

    if (params.amp_database_dir) {
        def user_amp_db = file(params.amp_database_dir, checkIfExists: true)
        if (user_amp_db.name.contains(params.amp_database) && file("${user_amp_db}/mmseqs2/ref_DB.dbtype").isFile()) {
            log.info "Using pre-indexed AMP database at: ${user_amp_db}"
            ch_amp_db_dir = channel.value(user_amp_db)
        } else {
            log.info "Preparing indexed AMP database from ${user_amp_db} via PREPARE_AMP_DB"
            PREPARE_AMP_DB(params.amp_database, user_amp_db, ch_amp_trigger)
            ch_amp_db_dir = PREPARE_AMP_DB.out.db.first()
        }
    } else {
        def cached_amp_db = file("${params.clystere_db_dir}/amp_${params.amp_database}_database")
        if (file("${cached_amp_db}/mmseqs2/ref_DB.dbtype").isFile()) {
            log.info "Using cached AMP database at: ${cached_amp_db}"
            ch_amp_db_dir = channel.value(cached_amp_db)
        } else {
            log.info "Preparing AMP (${params.amp_database}) database once via PREPARE_AMP_DB"
            PREPARE_AMP_DB(params.amp_database, [], ch_amp_trigger)
            ch_amp_db_dir = PREPARE_AMP_DB.out.db.first()
        }
    }

    AMPCOMBI2_PARSETABLES(
        ch_ampcombi_in.amp_input,
        ch_ampcombi_in.faa,
        ch_ampcombi_in.gbk,
        params.amp_database,
        ch_amp_db_dir,
        []
    )

    ch_amp_tables = AMPCOMBI2_PARSETABLES.out.tsv
        .map { _meta, tsv -> tsv }
        .collect()
    AMPCOMBI2_COMPLETE(ch_amp_tables)
    AMPCOMBI2_CLUSTER(AMPCOMBI2_COMPLETE.out.tsv)

    emit:
    amp_summary_tsv = AMPCOMBI2_COMPLETE.out.tsv
    amp_cluster_tsv = AMPCOMBI2_CLUSTER.out.cluster_tsv
}
