process PROBE_POSITIONS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(paf)

    output:
    tuple val(meta), path("*.probe_pos.tsv")   , emit: positions
    tuple val(meta), path("*.probe_filter.txt"), emit: report
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    probe_positions.py \\
        --paf ${paf} \\
        --out ${prefix}.probe_pos.tsv \\
        --report ${prefix}.probe_filter.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.probe_pos.tsv ${prefix}.probe_filter.txt versions.yml
    """
}
