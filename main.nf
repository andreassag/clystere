#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANTISMASH_BIGSCAPE } from './workflows/clystere'

//
// MAIN WORKFLOW
//
workflow {
    ANTISMASH_BIGSCAPE()
}

workflow.onComplete {
    log.info ( workflow.success ? "\nPipeline completed successfully!\n" : "Pipeline failed." )
}
