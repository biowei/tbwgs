process SNP_CLUSTER {
    tag "${meta.id}"
    label 'process_single'

    // Plain Python 3; the TB-Profiler image already carries it
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/tb-profiler:6.7.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/tb-profiler:6.7.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(dists)

    output:
    tuple val(meta), path("*.clusters.tsv"), emit: tsv
    tuple val("${task.process}"), val('snp_cluster.py'), eval("snp_cluster.py --version"), emit: versions_snp_cluster, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snp_cluster.py \\
        ${dists} \\
        ${prefix}.clusters.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clusters.tsv
    """
}
