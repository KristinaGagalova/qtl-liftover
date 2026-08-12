/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: LIFT_WGA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Route 1 - transfer QTL intervals by whole-genome DNA alignment:
      1. minimap2 source vs target        (source is the QUERY, so the PAF
                                           encodes source -> target)
      2. lift each interval through the CIGAR by tiled anchor points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'
include { PAF_LIFT_BED   } from '../../modules/local/paf_lift_bed'

workflow LIFT_WGA {

    take:
    ch_source // channel: [ val(meta), path(source genome fasta) ]
    ch_target // channel: [ val(meta), path(target genome fasta) ]
    ch_qtl    // channel: path(QTL bed, source coordinates)

    main:
    def ch_versions = Channel.empty()

    //
    // STEP 1: whole-genome alignment, source as query
    //
    MINIMAP2_ALIGN(ch_source, ch_target)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //
    // STEP 2: lift the intervals through the alignment
    //
    PAF_LIFT_BED(MINIMAP2_ALIGN.out.paf, ch_qtl)
    ch_versions = ch_versions.mix(PAF_LIFT_BED.out.versions.first())

    emit:
    paf       = MINIMAP2_ALIGN.out.paf     // channel: [ val(meta), path(paf) ]
    intervals = PAF_LIFT_BED.out.intervals // channel: [ val(meta), path(tsv) ]
    segments  = PAF_LIFT_BED.out.segments  // channel: [ val(meta), path(bed) ]
    versions  = ch_versions                // channel: path(versions.yml)
}
