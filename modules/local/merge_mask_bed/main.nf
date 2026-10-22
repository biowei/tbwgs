process MERGE_MASK_BED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(quantized_bed)
    path mask_bed

    output:
    tuple val(meta), path("*.mask.bed"), emit: bed
    tuple val("${task.process}"), val('grep'), eval("grep --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_grep

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // LOW bins of the mosdepth quantized BED, plus the fixed mask; grep exits 1 when no bin is LOW
    """
    {
        zcat ${quantized_bed} | { grep -P '\\tLOW\$' || true; } | cut -f 1-3
        grep -v '^#' ${mask_bed} | cut -f 1-3
    } | sort -k1,1 -k2,2n > ${prefix}.mask.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mask.bed
    """
}
