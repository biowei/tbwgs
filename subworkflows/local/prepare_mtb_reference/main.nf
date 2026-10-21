//
// Prepare one frozen reference, index, mask and tbdb directory for the whole run
//

include { UNTAR as UNTAR_TBDB } from '../../../modules/nf-core/untar/main'
include { SAMTOOLS_FAIDX      } from '../../../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX           } from '../../../modules/nf-core/bwa/index/main'

workflow PREPARE_MTB_REFERENCE {

    take:
    tbdb      // string: tbdb archive (.tar.gz) or directory; may be null when reference and mask are given
    reference // string: H37Rv FASTA; null to use genome.fasta from tbdb
    mask      // string: mask BED; null to use mask.bed from tbdb

    main:

    //
    // tbdb directory, extracted once
    //
    def ch_tbdb = channel.empty()
    if (tbdb && (tbdb.endsWith('.tar.gz') || tbdb.endsWith('.tgz'))) {
        UNTAR_TBDB(channel.of([ [ id: 'tbdb' ], file(tbdb, checkIfExists: true) ]))
        ch_tbdb = UNTAR_TBDB.out.untar.map { _meta, db -> db }.first()
    } else if (tbdb) {
        ch_tbdb = channel.value(file(tbdb, checkIfExists: true))
    }

    //
    // Reference FASTA and mask: explicit parameter first, otherwise the copy in tbdb
    //
    def ch_fasta = reference
        ? channel.value([ [ id: 'H37Rv' ], file(reference, checkIfExists: true) ])
        : ch_tbdb.map { db -> [ [ id: 'H37Rv' ], db.resolve('genome.fasta') ] }

    def ch_mask = mask
        ? channel.value(file(mask, checkIfExists: true))
        : ch_tbdb.map { db -> db.resolve('mask.bed') }

    SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }, false)

    BWA_INDEX(ch_fasta)

    emit:
    fasta = ch_fasta                        // channel: [ val(meta), path(fasta) ]
    fai   = SAMTOOLS_FAIDX.out.fai.first()  // channel: [ val(meta), path(fai) ]
    index = BWA_INDEX.out.index.first()     // channel: [ val(meta), path(bwa) ]
    mask  = ch_mask                         // channel: path(bed)
    tbdb  = ch_tbdb                         // channel: path(tbdb_dir)
}
