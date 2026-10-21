//
// Map trimmed reads to H37Rv, then sort, mark duplicates and index
//

include { BWA_MEM           } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORMADUP } from '../../../modules/nf-core/samtools/sormadup/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'

workflow MAP_READS_TO_H37RV {

    take:
    ch_reads // channel: [ val(meta), [ trimmed_1, trimmed_2 ] ]
    ch_index // channel: [ val(meta), path(bwa) ]   value channel
    ch_fasta // channel: [ val(meta), path(fasta) ] value channel
    ch_fai   // channel: [ val(meta), path(fai) ]   value channel

    main:

    //
    // Unsorted BAM straight from bwa mem; samtools collate in SORMADUP groups it by name
    //
    BWA_MEM(ch_reads, ch_index, ch_fasta, false)

    def ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { meta, fasta, _meta_fai, fai -> [ meta, fasta, fai ] }
        .first()

    SAMTOOLS_SORMADUP(BWA_MEM.out.bam, ch_fasta_fai)

    SAMTOOLS_INDEX(SAMTOOLS_SORMADUP.out.bam)

    def ch_bam = SAMTOOLS_SORMADUP.out.bam
        .join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)

    emit:
    bam     = ch_bam                         // channel: [ val(meta), path(bam), path(bai) ]
    metrics = SAMTOOLS_SORMADUP.out.metrics  // channel: [ val(meta), path(metrics) ]
}
