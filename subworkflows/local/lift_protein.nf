/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: LIFT_PROTEIN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Route 2 - transfer QTL intervals by protein anchors:
      1. source GFF -> representative transcript BED + protein lengths
      2. index the target genome for miniprot
      3. spliced-align the source proteins onto the target
      4. genes inside each QTL become anchors; confident, near-unique,
         collinear anchors define the transferred interval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GFF2PROTEIN_BED } from '../../modules/local/gff2protein_bed'
include { MINIPROT_INDEX  } from '../../modules/local/miniprot_index'
include { MINIPROT_ALIGN  } from '../../modules/local/miniprot_align'
include { MINIPROT_LIFT   } from '../../modules/local/miniprot_lift'

workflow LIFT_PROTEIN {

    take:
    ch_gff    // channel: [ val(meta), path(source annotation gff3) ]
    ch_pep    // channel: [ val(meta), path(source proteins fasta) ]
    ch_target // channel: [ val(meta), path(target genome fasta) ]
    ch_qtl    // channel: path(QTL bed, source coordinates)

    main:
    def ch_versions = Channel.empty()

    //
    // STEP 1: representative transcript coordinates in the source assembly
    //
    GFF2PROTEIN_BED(ch_gff)
    ch_versions = ch_versions.mix(GFF2PROTEIN_BED.out.versions.first())

    //
    // STEP 2 + 3: index the target and place the source proteins on it
    //
    MINIPROT_INDEX(ch_target)
    MINIPROT_ALIGN(ch_pep, MINIPROT_INDEX.out.index)
    ch_versions = ch_versions.mix(MINIPROT_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(MINIPROT_ALIGN.out.versions.first())

    //
    // STEP 4: anchors -> transferred intervals
    //   join on meta.id so the gff, gene bed and length table stay paired
    //
    def ch_lift_in = MINIPROT_ALIGN.out.gff
        .join(GFF2PROTEIN_BED.out.bed)
        .join(GFF2PROTEIN_BED.out.protein_len)

    MINIPROT_LIFT(ch_lift_in, ch_qtl)
    ch_versions = ch_versions.mix(MINIPROT_LIFT.out.versions.first())

    emit:
    gff       = MINIPROT_ALIGN.out.gff     // channel: [ val(meta), path(gff) ]
    intervals = MINIPROT_LIFT.out.intervals // channel: [ val(meta), path(tsv) ]
    anchors   = MINIPROT_LIFT.out.anchors  // channel: [ val(meta), path(tsv) ]
    versions  = ch_versions                // channel: path(versions.yml)
}
