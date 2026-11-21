#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLYSTERE                } from './workflows/clystere'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_clystere_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_clystere_pipeline'

//
// MAIN WORKFLOW ENTRYPOINT
//
workflow {
    PIPELINE_INITIALISATION(
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    CLYSTERE(PIPELINE_INITIALISATION.out.samplesheet)

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.monochrome_logs,
        params.hook_url,
        CLYSTERE.out.multiqc_report
    )
}
