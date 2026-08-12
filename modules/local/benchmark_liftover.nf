process BENCHMARK_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(probes_src), path(probes_tgt), path(lifted)
    path bed
    path wga_intervals       // assets/NO_FILE when the WGA route was skipped
    path miniprot_intervals  // assets/NO_FILE when the protein route was skipped

    output:
    tuple val(meta), path("summary.txt")            , emit: summary
    tuple val(meta), path("interval_eval_*.tsv")    , emit: interval_eval, optional: true
    tuple val(meta), path("probe_level_errors.tsv") , emit: probe_errors , optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def wga_arg = wga_intervals.name      == 'NO_FILE' ? '' : "--wga-intervals ${wga_intervals}"
    def mp_arg  = miniprot_intervals.name == 'NO_FILE' ? '' : "--miniprot-intervals ${miniprot_intervals}"
    """
    benchmark.py \\
        --probes-src ${probes_src} \\
        --probes-tgt ${probes_tgt} \\
        --lifted-probes ${lifted} \\
        --qtl ${bed} \\
        ${wga_arg} \\
        ${mp_arg} \\
        --outdir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch summary.txt versions.yml
    """
}
