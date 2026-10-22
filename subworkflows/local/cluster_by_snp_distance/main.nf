//
// SNP distances and single-linkage clusters for PASS samples of this batch and historical consensus sequences
//

include { MAKE_SNP_ALIGNMENT } from '../../../modules/local/make_snp_alignment/main'
include { SNPSITES    } from '../../../modules/nf-core/snpsites/main'
include { SNPDISTS    } from '../../../modules/nf-core/snpdists/main'
include { SNP_CLUSTER } from '../../../modules/local/snp_cluster/main'

workflow CLUSTER_BY_SNP_DISTANCE {

    take:
    ch_consensus // channel: [ val(meta), path(fa) ] meta carries mtb_qc
    ch_fasta     // channel: [ val(meta), path(fasta) ] value channel, H37Rv
    historical   // string: directory of historical consensus FASTA files (*.fa, *.fasta); may be null

    main:

    //
    // Only PASS samples are clustered; FAIL samples stay in the report without a cluster
    //
    def ch_batch = ch_consensus
        .filter { meta, _fa -> meta.mtb_qc == 'PASS' }
        .map { meta, fa -> [ meta.id, fa ] }

    def ch_historical = historical
        ? channel.fromPath("${historical}/*.{fa,fasta}", checkIfExists: true).map { fa -> [ fa.baseName, fa ] }
        : channel.empty()

    def ch_samples = ch_batch.mix(ch_historical)

    //
    // All sequences share the H37Rv coordinates, so concatenating them is the alignment.
    // H37Rv itself goes in too, so snp-sites always finds variable sites, even when every sample is identical;
    // SNP_CLUSTER leaves it out again.
    //
    // Fewer than two samples: nothing to compare, so M5 is skipped
    def ch_alignment_input = ch_samples
        .toList()
        .filter { samples -> samples.size() >= 2 }
        .map { samples -> [ samples ] }
        .combine(ch_fasta)
        .map { samples, _meta, fasta ->
            def sequences = (samples + [ [ referenceName(), fasta ] ]).sort { entry -> entry[0] }
            [ [ id: 'snp' ], sequences.collect { entry -> entry[0] }, sequences.collect { entry -> entry[1] } ]
        }

    MAKE_SNP_ALIGNMENT(ch_alignment_input)

    SNPSITES(MAKE_SNP_ALIGNMENT.out.fasta.map { _meta, alignment -> alignment })

    SNPDISTS(SNPSITES.out.fasta.map { alignment -> [ [ id: 'snp' ], alignment ] })

    SNP_CLUSTER(SNPDISTS.out.tsv)

    emit:
    dists    = SNPDISTS.out.tsv    // channel: [ val(meta), path(tsv) ] includes the H37Rv row
    clusters = SNP_CLUSTER.out.tsv // channel: [ val(meta), path(clusters.tsv) ] one column per threshold
}

//
// Name of the H37Rv record in the alignment; passed to SNP_CLUSTER as --exclude
//
def referenceName() {
    return '_reference_H37Rv'
}
