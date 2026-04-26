#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLYSTERE } from './workflows/clystere'

//
// MAIN WORKFLOW
//
workflow {
    CLYSTERE()
}
