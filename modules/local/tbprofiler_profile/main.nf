// Adapted from nf-core/modules tbprofiler/profile (TB-Profiler 3.0.8) for TB-Profiler 6.7.0 and a frozen database
process TBPROFILER_PROFILE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/tb-profiler:6.7.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/tb-profiler:6.7.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    path db, stageAs: 'tbdb_dir/*'

    output:
    tuple val(meta), path("results/*.results.json"), emit: json
    tuple val(meta), path("results/*.results.txt") , emit: txt, optional: true
    tuple val(meta), path("results/*.results.csv") , emit: csv, optional: true
    tuple val(meta), path("bam/*.bam")             , emit: bam
    tuple val(meta), path("vcf/*.vcf.gz")          , emit: vcf
    tuple val("${task.process}"), val('tb-profiler'), eval("tb-profiler --version | sed 's/.* //'"), emit: versions_tb_profiler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The database is a directory named after itself inside --db_dir
    """
    tb-profiler \\
        profile \\
        ${args} \\
        --read1 ${reads[0]} \\
        --read2 ${reads[1]} \\
        --db_dir tbdb_dir \\
        --db ${db.name} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        --temp .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bam results vcf
    touch bam/${prefix}.bam
    echo '{}' > results/${prefix}.results.json
    touch results/${prefix}.results.txt
    touch results/${prefix}.results.csv
    echo "" | gzip > vcf/${prefix}.targets.vcf.gz
    """
}
