//
// Measure coverage and flag each sample PASS or FAIL; no sample is dropped
//

include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth/main'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'

workflow CHECK_MAPPING_COVERAGE {

    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta        // channel: [ val(meta), path(fasta) ] value channel
    ch_fai          // channel: [ val(meta), path(fai) ]   value channel
    min_mean_depth  // number: minimum mean depth (x)
    min_breadth     // number: minimum percent of the genome at min_site_depth or more
    min_site_depth  // integer: depth below which a site counts as low coverage

    main:

    //
    // Quantized BED splits the genome into LOW (< min_site_depth) and PASS bins; M3 masks the LOW bins
    //
    MOSDEPTH(
        ch_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] },
        ch_fasta,
        [ 'LOW', 'PASS' ]
    )

    def ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { meta, fasta, _meta_fai, fai -> [ meta, fasta, fai ] }
        .first()

    SAMTOOLS_STATS(ch_bam, ch_fasta_fai)

    //
    // Add mtb_mean_depth, mtb_breadth and mtb_qc to meta
    //
    def ch_qc = MOSDEPTH.out.summary_txt
        .join(MOSDEPTH.out.global_txt, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, summary, global_dist ->
            def mean_depth = mosdepthMeanDepth(summary)
            def breadth    = mosdepthBreadth(global_dist, min_site_depth)
            def qc         = mean_depth >= min_mean_depth && breadth >= min_breadth ? 'PASS' : 'FAIL'
            [ meta.id, meta + [ mtb_mean_depth: mean_depth, mtb_breadth: breadth, mtb_qc: qc ] ]
        }

    def ch_bam_qc = ch_bam
        .map { meta, bam, bai -> [ meta.id, bam, bai ] }
        .join(ch_qc, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, bam, bai, meta -> [ meta, bam, bai ] }

    def ch_lowcov_bed = MOSDEPTH.out.quantized_bed
        .map { meta, bed -> [ meta.id, bed ] }
        .join(ch_qc, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, bed, meta -> [ meta, bed ] }

    def ch_multiqc_files = channel.empty()
        .mix(MOSDEPTH.out.global_txt.map { _meta, txt -> txt })
        .mix(MOSDEPTH.out.summary_txt.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS.out.stats.map { _meta, stats -> stats })

    emit:
    bam           = ch_bam_qc         // channel: [ val(meta + qc), path(bam), path(bai) ]
    quantized_bed = ch_lowcov_bed     // channel: [ val(meta + qc), path(quantized.bed.gz) ] LOW / PASS bins
    multiqc_files = ch_multiqc_files  // channel: [ path(file) ]
}

//
// Mean depth from the "total" row of mosdepth summary.txt (columns: chrom, length, bases, mean, min, max)
//
def mosdepthMeanDepth(summary) {
    def total = summary.readLines().find { line -> line.startsWith('total\t') }
    return total ? (total.split('\t')[3] as BigDecimal).round(2) : 0.0
}

//
// Percent of the genome at min_depth or more, from the "total" rows of mosdepth global.dist.txt
// (columns: chrom, depth, cumulative fraction of bases at that depth or more)
//
def mosdepthBreadth(global_dist, min_depth) {
    def row = global_dist.readLines()
        .collect { line -> line.split('\t') }
        .find { cols -> cols[0] == 'total' && (cols[1] as int) == (min_depth as int) }
    return row ? ((row[2] as BigDecimal) * 100).round(2) : 0.0
}
