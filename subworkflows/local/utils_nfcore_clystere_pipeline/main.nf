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
        file("${ssDir}/reads/${f.name}"),
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
    reads_input

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

    // 4. Parse input samplesheet via nf-schema (supports optional annotation, fastq_1, fastq_2, assay_type)
    def ssDir = file(input).parent
    ch_parsed_input = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))

    ch_samplesheet = ch_parsed_input.map { row ->
        def meta = row[0]
        def resolved_genome = resolveInputPath(row[1], ssDir)
        def resolved_annotation = (row.size() > 2 && row[2]) ? resolveInputPath(row[2], ssDir) : []
        [meta, resolved_genome, resolved_annotation]
    }

    ch_inline_reads = ch_parsed_input
        .filter { row -> row.size() > 3 && row[3] }
        .map { row ->
            def meta = row[0]
            def fq1 = resolveInputPath(row[3], ssDir)
            def fq2 = (row.size() > 4 && row[4]) ? resolveInputPath(row[4], ssDir) : null
            def assay = (row.size() > 5 && row[5]) ? row[5].toString() : 'metagenome'
            def read_meta = [id: meta.id, single_end: (fq2 == null), assay_type: assay]
            def reads_list = fq2 ? [fq1, fq2] : [fq1]
            [read_meta, reads_list]
        }

    // 5. Parse optional standalone --reads_input samplesheet
    ch_extra_reads = channel.empty()
    if (reads_input) {
        def readsDir = file(reads_input).parent
        ch_extra_reads = channel.fromList(samplesheetToList(reads_input, "${projectDir}/assets/schema_reads.json"))
            .map { row ->
                def meta = row[0]
                def fq1 = resolveInputPath(row[1], readsDir)
                def fq2 = (row.size() > 2 && row[2]) ? resolveInputPath(row[2], readsDir) : null
                def assay = (row.size() > 3 && row[3]) ? row[3].toString() : 'metagenome'
                def read_meta = [id: meta.id, single_end: (fq2 == null), assay_type: assay]
                def reads_list = fq2 ? [fq1, fq2] : [fq1]
                [read_meta, reads_list]
            }
    }

    ch_reads = ch_inline_reads.mix(ch_extra_reads)

    emit:
    samplesheet = ch_samplesheet
    reads       = ch_reads
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
