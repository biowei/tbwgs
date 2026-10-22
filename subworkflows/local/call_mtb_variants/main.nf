//
// Call haploid variants against H37Rv and build a masked consensus sequence per sample
//

include { BCFTOOLS_MPILEUP   } from '../../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM      } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER    } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_STATS     } from '../../../modules/nf-core/bcftools/stats/main'
include { MERGE_MASK_BED     } from '../../../modules/local/merge_mask_bed/main'
include { BCFTOOLS_CONSENSUS } from '../../../modules/nf-core/bcftools/consensus/main'

workflow CALL_MTB_VARIANTS {

    take:
    ch_bam           // channel: [ val(meta), path(bam), path(bai) ]
    ch_quantized_bed // channel: [ val(meta), path(quantized.bed.gz) ] mosdepth LOW / PASS bins
    ch_fasta         // channel: [ val(meta), path(fasta) ] value channel
    ch_fai           // channel: [ val(meta), path(fai) ]   value channel
    ch_mask          // channel: path(mask.bed)              value channel

    main:

    def ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { meta, fasta, _meta_fai, fai -> [ meta, fasta, fai ] }
        .first()

    //
    // Every sample is called, PASS or FAIL; only M5 looks at mtb_qc
    //
    BCFTOOLS_MPILEUP(
        ch_bam.map { meta, bam, _bai -> [ meta, bam, [], [] ] },
        ch_fasta_fai,
        false
    )

    BCFTOOLS_NORM(
        BCFTOOLS_MPILEUP.out.vcf.join(BCFTOOLS_MPILEUP.out.index, failOnMismatch: true, failOnDuplicate: true),
        ch_fasta
    )

    BCFTOOLS_FILTER(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.index, failOnMismatch: true, failOnDuplicate: true)
    )

    def ch_vcf = BCFTOOLS_FILTER.out.vcf
        .join(BCFTOOLS_FILTER.out.index, failOnMismatch: true, failOnDuplicate: true)

    BCFTOOLS_STATS(ch_vcf, [[], []], [[], []], [[], []], [[], []], [[], []])

    //
    // Mask = fixed mask.bed + sites below mtb_min_site_depth
    //
    MERGE_MASK_BED(ch_quantized_bed, ch_mask)

    def ch_consensus_input = ch_vcf
        .map { meta, vcf, tbi -> [ meta.id, meta, vcf, tbi ] }
        .join(MERGE_MASK_BED.out.bed.map { meta, bed -> [ meta.id, bed ] }, failOnMismatch: true, failOnDuplicate: true)
        .combine(ch_fasta.map { _meta, fasta -> fasta })
        .map { _id, meta, vcf, tbi, bed, fasta -> [ meta, vcf, tbi, fasta, bed ] }

    BCFTOOLS_CONSENSUS(ch_consensus_input)

    def ch_multiqc_files = BCFTOOLS_STATS.out.stats.map { _meta, stats -> stats }

    emit:
    vcf           = ch_vcf                        // channel: [ val(meta), path(vcf.gz), path(tbi) ]
    consensus     = BCFTOOLS_CONSENSUS.out.fasta  // channel: [ val(meta), path(fa) ]
    mask          = MERGE_MASK_BED.out.bed        // channel: [ val(meta), path(mask.bed) ]
    multiqc_files = ch_multiqc_files              // channel: [ path(stats) ]
}
