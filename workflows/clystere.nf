/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    clystere Main Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    antiSMASH + GECCO + DeepBGC + comBGC + BiG-SCAPE / BiG-SLiCE pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BGC_ANNOTATION                                    } from '../subworkflows/local/bgc_annotation'
include { BGC_CLUSTERING                                    } from '../subworkflows/local/bgc_clustering'
include { SUMMARY                                           } from '../modules/local/summary/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'

workflow CLYSTERE {
  take:
  ch_samples // channel: [ val(meta), path(genome) ]

  main:
  if (params.bigscape_run && params.bigslice_run) {
    error "Options --bigscape_run and --bigslice_run are mutually exclusive. Enable only one."
  }

  ch_multiqc_files = channel.empty()

  //
  // 1. Run BGC Annotations (antiSMASH, GECCO, DeepBGC)
  //
  BGC_ANNOTATION(ch_samples)

  //
  // 2. Run comBGC Unification & GCF Clustering (BiG-SCAPE / BiG-SLiCE)
  //
  BGC_CLUSTERING(
    BGC_ANNOTATION.out.antismash_dirs,
    BGC_ANNOTATION.out.gecco_output_dir,
    BGC_ANNOTATION.out.deepbgc_bgc_gbk,
    BGC_ANNOTATION.out.gecco_clusters,
    BGC_ANNOTATION.out.deepbgc_tsv,
    BGC_ANNOTATION.out.antismash_db
  )

  //
  // 3. Tabulation / Summary tables across all samples
  //
  if (params.run_tabulation) {
    ch_all_dirs = BGC_ANNOTATION.out.antismash_dirs
      .map { entry -> entry[1] }
      .collect()
    SUMMARY(ch_all_dirs)
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARY.out.all_regions_tsv.collect())
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARY.out.region_counts_tsv.collect())
  }

  //
  // 4. Software Version Collection
  //
  ch_collated_versions = softwareVersionsToYAML(channel.topic('versions'))
  ch_versions_file = ch_collated_versions
    .collectFile(name: 'versions.yml', newLine: true)
  ch_multiqc_files = ch_multiqc_files.mix(ch_versions_file)

  //
  // 5. MultiQC Reporting
  //
  ch_multiqc_config = params.multiqc_config ? file(params.multiqc_config) : []
  ch_multiqc_logo = params.multiqc_logo ? file(params.multiqc_logo) : []
  ch_multiqc_input = ch_multiqc_files
    .collect()
    .map { files -> [[id: 'multiqc_report'], files, ch_multiqc_config, ch_multiqc_logo, [], []] }

  MULTIQC(ch_multiqc_input)

  emit:
  multiqc_report = MULTIQC.out.report
}
