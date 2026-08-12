process GFF2PROTEIN_BED {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.genes.bed")      , emit: bed
    tuple val(meta), path("*.protein_len.tsv"), emit: protein_len
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gff_to_protein_bed.py \\
        --gff ${gff} \\
        --out-bed ${prefix}.genes.bed \\
        --out-len ${prefix}.protein_len.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genes.bed ${prefix}.protein_len.tsv versions.yml
    """
}
