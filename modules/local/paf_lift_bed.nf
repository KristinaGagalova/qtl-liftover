process PAF_LIFT_BED {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(paf)
    path bed

    output:
    tuple val(meta), path("*.qtl_target.tsv") , emit: intervals
    tuple val(meta), path("*.qtl_segments.bed"), emit: segments
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paf_lift.py bed \\
        --paf ${paf} \\
        --bed ${bed} \\
        --out ${prefix}.qtl_target.tsv \\
        --segments-out ${prefix}.qtl_segments.bed \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qtl_target.tsv ${prefix}.qtl_segments.bed versions.yml
    """
}
