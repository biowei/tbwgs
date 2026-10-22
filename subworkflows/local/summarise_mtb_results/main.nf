//
// One row per sample: coverage QC, MTBC fraction, lineage, per-drug calls and SNP clusters
//

include { MTB_SUMMARY } from '../../../modules/local/mtb_summary/main'

workflow SUMMARISE_MTB_RESULTS {

    take:
    ch_mtbc_fraction // channel: [ val(meta), val(percent MTBC reads) ] every sample
    ch_qc            // channel: [ val(meta + qc), ... ] meta carries mtb_mean_depth, mtb_breadth, mtb_qc
    ch_tbprofiler    // channel: [ val(meta), path(results.json) ]
    ch_clusters      // channel: [ val(meta), path(clusters.tsv) ]; empty when clustering was skipped
    ch_tbdb          // channel: path(tbdb_dir) value channel

    main:

    //
    // Every sample from preprocessing keeps its row; a sample without a coverage result gets empty QC fields
    //
    def ch_samples = ch_mtbc_fraction
        .map { meta, fraction -> [ meta.id, fraction ] }
        .join(ch_qc.map { entry -> [ entry[0].id, entry[0] ] }, remainder: true)
        .map { id, fraction, meta -> sampleRow(id, fraction, meta) }
        .collectFile(name: 'samples.tsv', seed: sampleHeader(), newLine: true, sort: true)
        .map { tsv -> [ [ id: 'mtb_summary' ], tsv ] }

    MTB_SUMMARY(
        ch_samples,
        ch_tbprofiler.map { _meta, json -> json }.collect().ifEmpty([]),
        ch_clusters.map { _meta, tsv -> tsv }.collect().ifEmpty([]),
        ch_tbdb
    )

    emit:
    tsv           = MTB_SUMMARY.out.tsv                                // channel: [ val(meta), path(mtb_summary.tsv) ]
    multiqc_files = MTB_SUMMARY.out.mqc.map { _meta, tsv -> tsv }      // channel: [ path(mtb_summary_mqc.tsv) ]
}

def sampleHeader() {
    return [ 'sample', 'mtbc_fraction', 'mean_depth', 'breadth', 'coverage_qc' ].join('\t')
}

//
// A sample with no coverage result (meta is null after the remainder join) keeps its ID and MTBC fraction
//
def sampleRow(id, fraction, meta) {
    def depth   = meta ? String.format('%.2f', meta.mtb_mean_depth as double) : ''
    def breadth = meta ? String.format('%.2f', meta.mtb_breadth as double) : ''
    return [ id, fraction != null ? fraction : '', depth, breadth, meta?.mtb_qc ?: '' ].join('\t')
}
