/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: UTILS_NFCORE_CLYSTERE_PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline initialisation and completion subworkflows following nf-core standards.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFCORE_PIPELINE; completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE                  } from '../../nf-core/utils_nextflow_pipeline'
include { paramsHelp; paramsSummaryLog; validateParameters; samplesheetToList } from 'plugin/nf-schema'

def resolveInputPath(rawPath, ssDir) {
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
    version
    help
    validate_params
    _monochrome_logs
    nextflow_cli_args
    outdir
    input

    main:
    // 1. Display help text and exit when --help is passed
    if (help) {
        def help_param = ['true', 'false'].contains(help.toString()) ? null : help.toString()
        log.info paramsHelp(
            parameter: help_param,
            command: "nextflow run ${workflow.manifest.name} --input samplesheet.csv --outdir results -profile <docker/singularity/apptainer/conda>"
        )
        System.exit(0)
    }

    // 2. Run standard Nextflow and nf-core utility checks (--version, params dump, profile check)
    def check_conda = workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    UTILS_NEXTFLOW_PIPELINE(version, true, outdir, check_conda)
    UTILS_NFCORE_PIPELINE(nextflow_cli_args)

    // 3. Validate parameters if enabled and log summary
    if (validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog()

    // 4. Parse input samplesheet via nf-schema
    def ssDir = file(input).parent
    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, genome, annotation ->
            def resolved_genome = resolveInputPath(genome, ssDir)
            def resolved_annotation = annotation ? resolveInputPath(annotation, ssDir) : []
            [meta, resolved_genome, resolved_annotation]
        }

    emit:
    samplesheet = ch_samplesheet
}

workflow PIPELINE_COMPLETION {
    take:
    email
    email_on_fail
    monochrome_logs
    _hook_url
    _multiqc_report

    main:
    workflow.onComplete {
        completionSummary(monochrome_logs)
        if ((email || email_on_fail) && email) {
            log.info "Completion notification target configured for ${email}"
        }
    }
}
