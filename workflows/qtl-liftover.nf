nextflow.enable.dsl = 2

include { LIFT_WGA         } from '../subworkflows/local/lift_wga'
include { LIFT_PROTEIN     } from '../subworkflows/local/lift_protein'
include { BENCHMARK_PROBES } from '../subworkflows/local/benchmark_probes'

workflow QTL_LIFTOVER_FLOW {

    main:
    if (!params.source_genome) { error "Missing required parameter: --source_genome" }
    if (!params.target_genome) { error "Missing required parameter: --target_genome" }
    if (!params.qtl_bed)       { error "Missing required parameter: --qtl_bed" }

    def ch_versions = Channel.empty()
    def no_file     = file("${projectDir}/assets/NO_FILE")

    def meta_src = [id: file(params.source_genome).baseName]
    def meta_tgt = [id: file(params.target_genome).baseName]

    def ch_source = Channel.of([meta_src, file(params.source_genome, checkIfExists: true)])
    def ch_target = Channel.of([meta_tgt, file(params.target_genome, checkIfExists: true)])
    def ch_qtl    = Channel.value(file(params.qtl_bed, checkIfExists: true))

    //
    // ROUTE 1: whole-genome alignment
    //
    LIFT_WGA(ch_source, ch_target, ch_qtl)
    ch_versions = ch_versions.mix(LIFT_WGA.out.versions)

    //
    // ROUTE 2: protein anchors
    //
    def ch_mp_intervals = Channel.value(no_file)
    if (!params.skip_miniprot) {
        if (!params.source_gff || !params.source_proteins) {
            error "--source_gff and --source_proteins are required unless --skip_miniprot"
        }
        def ch_gff = Channel.of([meta_src, file(params.source_gff,      checkIfExists: true)])
        def ch_pep = Channel.of([meta_src, file(params.source_proteins, checkIfExists: true)])

        LIFT_PROTEIN(ch_gff, ch_pep, ch_target, ch_qtl)
        ch_versions     = ch_versions.mix(LIFT_PROTEIN.out.versions)
        ch_mp_intervals = LIFT_PROTEIN.out.intervals.map { meta, tsv -> tsv }
    }

    //
    // BENCHMARK: probe-derived truth set
    //
    if (!params.skip_benchmark) {
        if (!params.probes) { error "--probes is required unless --skip_benchmark" }
        def ch_probes = Channel.of([[id: 'probes'], file(params.probes, checkIfExists: true)])

        BENCHMARK_PROBES(
            ch_probes,
            ch_source,
            ch_target,
            LIFT_WGA.out.paf,
            ch_qtl,
            LIFT_WGA.out.intervals.map { meta, tsv -> tsv },
            ch_mp_intervals
        )
        ch_versions = ch_versions.mix(BENCHMARK_PROBES.out.versions)
    }

    emit:
    wga_intervals      = LIFT_WGA.out.intervals // channel: [ val(meta), path(tsv) ]
    wga_segments       = LIFT_WGA.out.segments  // channel: [ val(meta), path(bed) ]
    miniprot_intervals = ch_mp_intervals        // channel: path(tsv)
    versions           = ch_versions            // channel: path(versions.yml)
}
