//
// Merge lanes, trim reads and estimate the fraction of MTBC reads
//

include { CAT_FASTQ                    } from '../../../modules/nf-core/cat/fastq/main'
include { FASTQC                       } from '../../../modules/nf-core/fastqc/main'
include { FASTP                        } from '../../../modules/nf-core/fastp/main'
include { UNTAR as UNTAR_KRAKEN2DB     } from '../../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2              } from '../../../modules/nf-core/kraken2/kraken2/main'

workflow PREPROCESS_READS {

    take:
    ch_reads  // channel: [ val(meta), [ fastq_1, fastq_2, ... ] ] one item per sample, lanes in R1/R2 order
    kraken2db // string: path to a Kraken2 database directory or .tar.gz

    main:

    //
    // Merge lanes only for samples sequenced on more than one lane
    //
    def ch_reads_by_lanes = ch_reads
        .branch { _meta, fastqs ->
            single: fastqs.size() == 2
            multiple: true
        }

    CAT_FASTQ(ch_reads_by_lanes.multiple)

    def ch_merged_reads = CAT_FASTQ.out.reads.mix(ch_reads_by_lanes.single)

    //
    // Read QC and trimming
    //
    FASTQC(ch_merged_reads)

    FASTP(
        ch_merged_reads.map { meta, fastqs -> [ meta, fastqs, [] ] },
        false,
        false,
        false
    )

    //
    // Kraken2 on the trimmed reads
    //
    def ch_kraken2db
    if (kraken2db.endsWith('.tar.gz') || kraken2db.endsWith('.tgz')) {
        UNTAR_KRAKEN2DB(channel.of([ [ id: 'kraken2db' ], file(kraken2db, checkIfExists: true) ]))
        ch_kraken2db = UNTAR_KRAKEN2DB.out.untar.map { _meta, db -> db }.first()
    } else {
        ch_kraken2db = channel.value(file(kraken2db, checkIfExists: true))
    }

    KRAKEN2_KRAKEN2(FASTP.out.reads, ch_kraken2db, false, false)

    def ch_mtbc_fraction = KRAKEN2_KRAKEN2.out.report
        .map { meta, report -> [ meta, mtbcFraction(report) ] }

    def ch_multiqc_files = channel.empty()
        .mix(FASTQC.out.zip.map { _meta, zip -> zip })
        .mix(FASTP.out.json.map { _meta, json -> json })
        .mix(KRAKEN2_KRAKEN2.out.report.map { _meta, report -> report })

    emit:
    reads         = FASTP.out.reads   // channel: [ val(meta), [ trimmed_1, trimmed_2 ] ]
    mtbc_fraction = ch_mtbc_fraction  // channel: [ val(meta), val(percent of reads classified as MTBC) ]
    multiqc_files = ch_multiqc_files  // channel: [ path(file) ]
}

//
// Percent of reads assigned to the Mycobacterium tuberculosis complex (taxid 77643) in a Kraken2 report.
// Report columns: percent, clade reads, direct reads, rank code, taxid, name.
//
def mtbcFraction(report) {
    def clade_reads = [:]
    report.eachLine { line ->
        def cols = line.split('\t')
        if (cols.size() >= 6) {
            clade_reads[cols[4].trim()] = cols[1].trim() as long
        }
    }
    def total = (clade_reads['0'] ?: 0L) + (clade_reads['1'] ?: 0L)
    if (total == 0) {
        return 0.0
    }
    return ((clade_reads['77643'] ?: 0L) * 100.0 / total).round(2)
}
