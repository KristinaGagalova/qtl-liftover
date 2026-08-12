process PAF_LIFT_POINTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(paf)
    path positions   // TSV from PROBE_POSITIONS: probe chrom pos strand ident cov n_hits

    output:
    tuple val(meta), path("*.lifted.tsv"), emit: lifted
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tail -n +2 ${positions} \\
        | awk 'BEGIN{OFS="\\t"}{print \$2,\$3,\$1}' \\
        > points.tsv

    paf_lift.py points \\
        --paf ${paf} \\
        --tsv points.tsv \\
        --out ${prefix}.lifted.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lifted.tsv versions.yml
    """
}
