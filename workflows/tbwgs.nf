/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPROCESS_READS       } from '../subworkflows/local/preprocess_reads'
include { PREPARE_MTB_REFERENCE  } from '../subworkflows/local/prepare_mtb_reference'
include { MAP_READS_TO_H37RV     } from '../subworkflows/local/map_reads_to_h37rv'
include { CHECK_MAPPING_COVERAGE } from '../subworkflows/local/check_mapping_coverage'
include { CALL_MTB_VARIANTS      } from '../subworkflows/local/call_mtb_variants'
include { PROFILE_DRUG_RESISTANCE } from '../subworkflows/local/profile_drug_resistance'
include { CLUSTER_BY_SNP_DISTANCE } from '../subworkflows/local/cluster_by_snp_distance'
include { SUMMARISE_MTB_RESULTS  } from '../subworkflows/local/summarise_mtb_results'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_tbwgs_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TBWGS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    //
    // SUBWORKFLOW: P - merge lanes, read QC, trimming and MTBC fraction
    //
    PREPROCESS_READS(ch_samplesheet, params.kraken2db)
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_READS.out.multiqc_files)

    //
    // SUBWORKFLOW: M0 - frozen reference, BWA index, mask and tbdb
    //
    PREPARE_MTB_REFERENCE(params.mtb_tbdb, params.mtb_reference, params.mtb_mask_bed)

    //
    // SUBWORKFLOW: M1 - map to H37Rv, mark duplicates and index
    //
    MAP_READS_TO_H37RV(
        PREPROCESS_READS.out.reads,
        PREPARE_MTB_REFERENCE.out.index,
        PREPARE_MTB_REFERENCE.out.fasta,
        PREPARE_MTB_REFERENCE.out.fai
    )

    //
    // SUBWORKFLOW: M2 - coverage gate; flags samples PASS / FAIL in meta, drops none
    //
    CHECK_MAPPING_COVERAGE(
        MAP_READS_TO_H37RV.out.bam,
        PREPARE_MTB_REFERENCE.out.fasta,
        PREPARE_MTB_REFERENCE.out.fai,
        params.mtb_min_mean_depth,
        params.mtb_min_breadth,
        params.mtb_min_site_depth
    )
    ch_multiqc_files = ch_multiqc_files.mix(CHECK_MAPPING_COVERAGE.out.multiqc_files)

    //
    // SUBWORKFLOW: M3 - haploid variant calls and masked consensus, for every sample
    //
    CALL_MTB_VARIANTS(
        CHECK_MAPPING_COVERAGE.out.bam,
        CHECK_MAPPING_COVERAGE.out.quantized_bed,
        PREPARE_MTB_REFERENCE.out.fasta,
        PREPARE_MTB_REFERENCE.out.fai,
        PREPARE_MTB_REFERENCE.out.mask
    )
    ch_multiqc_files = ch_multiqc_files.mix(CALL_MTB_VARIANTS.out.multiqc_files)

    // Percent N in each consensus, for the clustering filter and the summary
    def ch_consensus = CALL_MTB_VARIANTS.out.consensus
        .map { meta, fa -> [ meta + [ mtb_consensus_n: consensusNPercent(fa) ], fa ] }

    //
    // SUBWORKFLOW: M4 - TB-Profiler drug resistance and lineage, in parallel with M1-M3
    //
    PROFILE_DRUG_RESISTANCE(PREPROCESS_READS.out.reads, PREPARE_MTB_REFERENCE.out.tbdb)

    //
    // SUBWORKFLOW: M5 - SNP distances and clusters for PASS samples plus historical consensus sequences
    //
    def ch_clusters = channel.empty()
    if (!params.mtb_skip_clustering) {
        CLUSTER_BY_SNP_DISTANCE(
            ch_consensus,
            PREPARE_MTB_REFERENCE.out.fasta,
            params.mtb_historical_consensus,
            params.mtb_max_consensus_n
        )
        ch_clusters = CLUSTER_BY_SNP_DISTANCE.out.clusters
    }

    //
    // SUBWORKFLOW: M6 - one row per sample, FAIL samples and samples without results included
    //
    SUMMARISE_MTB_RESULTS(
        PREPROCESS_READS.out.mtbc_fraction,
        CHECK_MAPPING_COVERAGE.out.bam,
        ch_consensus.map { meta, _fa -> [ meta, meta.mtb_consensus_n ] },
        PROFILE_DRUG_RESISTANCE.out.json,
        ch_clusters,
        PREPARE_MTB_REFERENCE.out.tbdb
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARISE_MTB_RESULTS.out.multiqc_files)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'tbwgs_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'tbwgs'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Percent of N (or n) in a consensus FASTA, sequence lines only
//
def consensusNPercent(fasta) {
    def total = 0L
    def n     = 0L
    fasta.eachLine { line ->
        if (!line.startsWith('>')) {
            total += line.length()
            n     += line.count('N') + line.count('n')
        }
    }
    return total ? 100.0 * n / total : 100.0
}
