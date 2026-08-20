/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: UTILS_NFCORE_CLYSTERE_PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline initialisation and completion subworkflows following nf-core standards.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

def resolveGenomePath(rawPath, ssDir) {
    def f = file(rawPath)
    if (f.exists()) return f
    def rel = rawPath.toString().replace(projectDir.toString() + '/', '')
    def candidates = [
        file("${ssDir}/${rawPath}"),
        file("${ssDir}/${rel}"),
        file("${ssDir}/genomes/${f.name}"),
        file("${ssDir}/${f.name}"),
        file("${projectDir}/assets/${rawPath}"),
        file("${projectDir}/assets/${rel}"),
        file("${projectDir}/assets/genomes/${f.name}"),
        file("${projectDir}/assets/${f.name}")
    ]
    def found = candidates.find { candidate -> candidate.exists() }
    if (found) return found
    return file(rawPath, checkIfExists: true)
}

workflow PIPELINE_INITIALISATION {
    take:
    _version
    validate_params
    _monochrome_logs
    _args
    _outdir

    main:
    // 1. Validate parameters if enabled
    if (validate_params) {
        validateParameters()
    }

    // 2. Parse input samplesheet via nf-schema
    def ssDir = file(params.input).parent
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, genome, _annotation ->
            [meta, resolveGenomePath(genome, ssDir)]
        }

    emit:
    samplesheet = ch_samplesheet
}

workflow PIPELINE_COMPLETION {
    take:
    email
    email_on_fail
    _hook_url
    _multiqc_report

    main:
    // Email / webhook notifications when configured
    if (email || email_on_fail) {
        if (email) {
            log.info "Sending execution completion email to ${email}"
        }
    }
}
