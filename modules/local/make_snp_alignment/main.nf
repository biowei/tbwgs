process MAKE_SNP_ALIGNMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), val(names), path(fastas, stageAs: 'seq*/*')

    output:
    tuple val(meta), path("*.aln.fa"), emit: fasta
    tuple val("${task.process}"), val('coreutils'), eval("tr --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Every sequence is on H37Rv coordinates: one record per input, renamed, sequence joined onto one line
    def records = [ names, fastas instanceof List ? fastas : [ fastas ] ].transpose().collect { name, fasta ->
        "    printf '>%s\\n' '${name}'; grep -v '^>' '${fasta}' | tr -d '\\n\\r'; printf '\\n'"
    }.join('\n')
    """
    {
${records}
    } > ${prefix}.aln.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln.fa
    """
}
