/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BENCHMARK_PROBES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Probes are mapped INDEPENDENTLY to both assemblies. Probes placing
    uniquely and confidently in both give (source pos, target pos) pairs
    that owe nothing to either liftover route - a method-free answer key.

      1. map probes to source and to target                (2 parallel jobs)
      2. keep only confident, single-locus placements
      3. lift the source probe positions through the WGA PAF
      4. score base-level error and per-QTL recall / precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN as MINIMAP2_PROBES } from '../../modules/local/minimap2_align'
include { PROBE_POSITIONS                   } from '../../modules/local/probe_positions'
include { PAF_LIFT_POINTS                   } from '../../modules/local/paf_lift_points'
include { BENCHMARK_LIFTOVER                } from '../../modules/local/benchmark_liftover'

workflow BENCHMARK_PROBES {

    take:
    ch_probes    // channel: [ val(meta), path(probe fasta) ]
    ch_source    // channel: [ val(meta), path(source genome fasta) ]
    ch_target    // channel: [ val(meta), path(target genome fasta) ]
    ch_wga_paf   // channel: [ val(meta), path(source -> target paf) ]
    ch_qtl       // channel: path(QTL bed)
    ch_wga_tsv   // channel: path(route 1 intervals, or assets/NO_FILE)
    ch_mp_tsv    // channel: path(route 2 intervals, or assets/NO_FILE)

    main:
    def ch_versions = Channel.empty()

    //
    // STEP 1: probes vs each assembly, tagged so the two runs stay distinct
    //
    def ch_probe_jobs = ch_probes
        .combine(ch_source.map { meta, fa -> ['src', fa] }
                 .mix(ch_target.map { meta, fa -> ['tgt', fa] }))
        .map { pmeta, probes, tag, genome ->
            tuple([id: "${pmeta.id}_${tag}", tag: tag], probes, genome)
        }

    MINIMAP2_PROBES(
        ch_probe_jobs.map { meta, probes, genome -> [meta, probes] },
        ch_probe_jobs.map { meta, probes, genome -> [meta, genome] }
    )
    ch_versions = ch_versions.mix(MINIMAP2_PROBES.out.versions.first())

    //
    // STEP 2: confident, uniquely placed probes only
    //
    PROBE_POSITIONS(MINIMAP2_PROBES.out.paf)
    ch_versions = ch_versions.mix(PROBE_POSITIONS.out.versions.first())

    def ch_pos = PROBE_POSITIONS.out.positions.branch { meta, tsv ->
        src: meta.tag == 'src'
        tgt: meta.tag == 'tgt'
    }
    def ch_src_pos = ch_pos.src.map { meta, tsv -> tsv }
    def ch_tgt_pos = ch_pos.tgt.map { meta, tsv -> tsv }

    //
    // STEP 3: lift the source probe positions through the WGA alignment
    //
    PAF_LIFT_POINTS(ch_wga_paf, ch_src_pos)
    ch_versions = ch_versions.mix(PAF_LIFT_POINTS.out.versions.first())

    //
    // STEP 4: score both routes against the truth set
    //
    def ch_bench_in = PAF_LIFT_POINTS.out.lifted
        .map { meta, lifted -> [meta, lifted] }
        .combine(ch_src_pos)
        .combine(ch_tgt_pos)
        .map { meta, lifted, src, tgt -> tuple(meta, src, tgt, lifted) }

    BENCHMARK_LIFTOVER(ch_bench_in, ch_qtl, ch_wga_tsv, ch_mp_tsv)
    ch_versions = ch_versions.mix(BENCHMARK_LIFTOVER.out.versions.first())

    emit:
    positions     = PROBE_POSITIONS.out.positions        // channel: [ val(meta), path(tsv) ]
    lifted        = PAF_LIFT_POINTS.out.lifted           // channel: [ val(meta), path(tsv) ]
    summary       = BENCHMARK_LIFTOVER.out.summary       // channel: [ val(meta), path(txt) ]
    interval_eval = BENCHMARK_LIFTOVER.out.interval_eval // channel: [ val(meta), path(tsv) ]
    versions      = ch_versions                          // channel: path(versions.yml)
}
