//
// Drug resistance and lineage with TB-Profiler against the frozen tbdb
//

include { TBPROFILER_PROFILE } from '../../../modules/local/tbprofiler_profile/main'

workflow PROFILE_DRUG_RESISTANCE {

    take:
    ch_reads // channel: [ val(meta), [ trimmed_1, trimmed_2 ] ]
    ch_tbdb  // channel: path(tbdb_dir) value channel

    main:

    //
    // Every sample is profiled straight from the trimmed reads, without waiting for the coverage gate
    //
    TBPROFILER_PROFILE(ch_reads, ch_tbdb)

    emit:
    json = TBPROFILER_PROFILE.out.json // channel: [ val(meta), path(results.json) ] resistance, lineage and per-gene depth
    txt  = TBPROFILER_PROFILE.out.txt  // channel: [ val(meta), path(results.txt) ]
    csv  = TBPROFILER_PROFILE.out.csv  // channel: [ val(meta), path(results.csv) ]
}
