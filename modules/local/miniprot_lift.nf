process MINIPROT_LIFT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(gff), path(genes_bed), path(protein_len)
    path bed

    output:
    tuple val(meta), path("*.qtl_target.tsv"), emit: intervals
    tuple val(meta), path("*.anchors.tsv")   , emit: anchors
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    miniprot_lift.py \\
        --miniprot ${gff} \\
        --gene-bed ${genes_bed} \\
        --prot-len ${protein_len} \\
        --qtl ${bed} \\
        --out ${prefix}.qtl_target.tsv \\
        --anchors-out ${prefix}.anchors.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qtl_target.tsv ${prefix}.anchors.tsv versions.yml
    """
}
