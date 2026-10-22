process MTB_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    // Plain Python 3; the TB-Profiler image already carries it
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/tb-profiler:6.7.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/tb-profiler:6.7.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(samples)
    path jsons, stageAs: 'tbprofiler/*'
    path clusters
    path tbdb

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.tsv"), emit: tsv
    tuple val(meta), path("*_mqc.tsv")                       , emit: mqc
    tuple val("${task.process}"), val('mtb_summary.py'), eval("mtb_summary.py --version"), emit: versions_mtb_summary, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tbprofiler = jsons ? "--tbprofiler ${jsons.join(' ')}" : ''
    def cluster_arg = clusters ? "--clusters ${clusters}" : ''
    """
    mtb_summary.py \\
        --samples ${samples} \\
        --genes ${tbdb}/genes.bed \\
        ${tbprofiler} \\
        ${cluster_arg} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv ${prefix}_mqc.tsv
    """
}
