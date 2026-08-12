#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QTL_LIFTOVER_FLOW } from './workflows/qtl-liftover'

workflow {
    QTL_LIFTOVER_FLOW()
}
