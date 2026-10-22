# biowei/tbwgs: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC, after lanes of the same sample are merged
- [fastp](#fastp) - Adapter and quality trimming
- [Kraken2](#kraken2) - Read classification, used to report the fraction of reads from the _M. tuberculosis_ complex (MTBC)
- [Alignment](#alignment) - Reads mapped to H37Rv, sorted, with duplicates marked
- [Coverage](#coverage) - Depth and breadth of coverage, and the per-sample coverage QC flag
- [Variants](#variants) - Haploid variant calls against H37Rv
- [Consensus](#consensus) - Masked consensus sequence per sample
- [TB-Profiler](#tb-profiler) - Drug resistance and lineage
- [SNP clustering](#snp-clustering) - Pairwise SNP distances and transmission clusters
- [MTB summary](#mtb-summary) - One row per sample with QC, lineage, per-drug calls and clusters
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*.fastp.html`: fastp report for each sample.
  - `*.fastp.json`: fastp statistics, also read by MultiQC.
  - `*.fastp.log`: fastp log.

</details>

[fastp](https://github.com/OpenGene/fastp) trims adapters and low-quality bases. Only the trimmed reads go to the downstream steps. The trimmed FASTQ files are not saved.

### Kraken2

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `*.kraken2.report.txt`: Kraken2 report for the trimmed reads of each sample.

</details>

[Kraken2](https://ccb.jhu.edu/software/kraken2/) classifies the trimmed reads against `--kraken2db`. The pipeline reads the MTBC line (taxid 77643) from each report and calculates the percentage of reads in the MTBC. This value is a contamination flag only. It does not remove samples.

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `*.markdup.bam`: trimmed reads mapped to H37Rv with BWA-MEM, sorted by coordinate, with duplicates marked by `samtools markdup`. The read group ID and sample name are the sample ID.
  - `*.markdup.bam.bai`: BAM index.
  - `*.markdup.metrics`: `samtools markdup` statistics.

</details>

[BWA-MEM](https://github.com/lh3/bwa) maps the reads to the reference from `--mtb_tbdb` or `--mtb_reference`. [SAMtools](https://www.htslib.org/) then runs `collate`, `fixmate`, `sort` and `markdup` in one pipe. Duplicates are marked, not removed.

### Coverage

<details markdown="1">
<summary>Output files</summary>

- `coverage/`
  - `*.mosdepth.summary.txt`: mean depth; the `total` row is used for the QC flag.
  - `*.mosdepth.global.dist.txt`: fraction of the genome at each depth or more.
  - `*.quantized.bed.gz` and `*.quantized.bed.gz.csi`: the genome split into `LOW` (below `--mtb_min_site_depth`) and `PASS` regions. The variant calling step uses the `LOW` regions to mask the consensus sequence with N.
  - `*.stats`: `samtools stats` output.

</details>

[mosdepth](https://github.com/brentp/mosdepth) measures coverage of the deduplicated BAM. Each sample is flagged `PASS` when its mean depth is at least `--mtb_min_mean_depth` (default 20×) and at least `--mtb_min_breadth` percent of the genome (default 95%) is covered at `--mtb_min_site_depth` or more (default 10×); otherwise it is flagged `FAIL`. The flag never removes a sample: `FAIL` samples still get variant calls and drug resistance results, and are only left out of SNP clustering. Both coverage reports and `samtools stats` are shown in MultiQC.

### Variants

<details markdown="1">
<summary>Output files</summary>

- `variants/`
  - `*.filtered.vcf.gz` and `*.filtered.vcf.gz.tbi`: SNPs and indels that pass the filter.
  - `*.bcftools_stats.txt`: `bcftools stats` of the filtered VCF, shown in MultiQC.

</details>

[BCFtools](https://samtools.github.io/bcftools/) calls variants on every sample, `PASS` or `FAIL`. `bcftools mpileup` uses base and mapping quality of at least 20, and `bcftools call` runs the multiallelic caller with `--ploidy 1`. Multiallelic sites are split and indels left-aligned with `bcftools norm`. `bcftools filter` keeps calls with `QUAL >= 20` and depth of at least `--mtb_min_site_depth`.

A haploid caller reports one allele per site, so mixed infections and heteroresistance are not visible here. The drug resistance step reports allele frequencies for resistance mutations.

### Consensus

<details markdown="1">
<summary>Output files</summary>

- `consensus/`
  - `*.fa`: H37Rv with the sample's SNPs applied and masked sites written as N. The sequence keeps the H37Rv length and header.
  - `*.mask.bed`: the regions written as N, which are the fixed mask from `--mtb_tbdb` or `--mtb_mask_bed` plus the sample's sites below `--mtb_min_site_depth`.

</details>

`bcftools consensus` applies SNPs only. Indels are left out so that all consensus sequences have the same length and can be compared position by position for SNP clustering.

### TB-Profiler

<details markdown="1">
<summary>Output files</summary>

- `tbprofiler/`
  - `*.results.json`: full TB-Profiler result: lineage, drug resistance type, resistance mutations with allele frequency and WHO confidence grading, the database name and commit, and the median depth of every resistance gene (`qc.target_qc`).
  - `*.results.txt` and `*.results.csv`: the same result as text and CSV reports.

</details>

[TB-Profiler](https://github.com/jodyphelan/TBProfiler) 6.7.0 runs on the trimmed reads of every sample, in parallel with mapping and variant calling, so it does not wait for the coverage QC flag. It uses its own mapping and variant calling (bwa, freebayes) against the reference in `--mtb_tbdb`, with `--no_trim` because fastp has already trimmed the reads.

TB-Profiler reports a drug as susceptible when it finds no resistance mutation, even if the gene is barely covered. Use the per-gene depth in the JSON before reading a "susceptible" call: the [MTB summary](#mtb-summary) marks such drugs as `uninterpretable`.

### SNP clustering

<details markdown="1">
<summary>Output files</summary>

- `clustering/`
  - `snp.tsv`: pairwise SNP distance matrix from snp-dists. It includes a `_reference_H37Rv` row and column, the distance of each sample to H37Rv.
  - `snp.clusters.tsv`: one row per sample and one `cluster_<t>snp` column per threshold in `--mtb_snp_thresholds`. Samples in the same cluster share an ID such as `T5-1`; a sample with no other sample within the threshold gets `-`.

</details>

Only samples with coverage QC `PASS` and at most `--mtb_max_consensus_n` percent N in the consensus (default 10%; the fixed mask alone is about 7%) are clustered, together with any consensus sequences in `--mtb_historical_consensus` (file name without extension is the sample ID). All consensus sequences are on H37Rv coordinates, so they are joined into one alignment with H37Rv itself, [snp-sites](https://github.com/sanger-pathogens/snp-sites) keeps the variable columns, and [snp-dists](https://github.com/tseemann/snp-dists) counts differences between A, C, G and T only, so a masked site (N) in either sample is not counted. Clusters are single-linkage: two samples are in the same cluster when they are linked by a chain of pairs each within the threshold (5 and 12 SNPs by default). The step is skipped with `--mtb_skip_clustering`, or when fewer than two samples are available.

### MTB summary

<details markdown="1">
<summary>Output files</summary>

- `summary/`
  - `mtb_summary.tsv`: one row per sample, also shown as a table in the MultiQC report.

</details>

Every sample in the samplesheet has a row, including samples with coverage QC `FAIL` and samples for which a step gave no result; the `notes` column says what is missing and why. The columns are:

| Column                                     | Content                                                                                                                                                                 |
| ------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                                   | sample ID                                                                                                                                                               |
| `mtbc_fraction`, `mtbc_flag`               | percent of reads that Kraken2 assigns to the _M. tuberculosis_ complex; `LOW_MTBC` when below `--mtb_min_mtbc_fraction`, `not_checked` when that parameter is not set   |
| `mean_depth`, `breadth`, `coverage_qc`     | coverage QC from mosdepth (see [Coverage](#coverage))                                                                                                                   |
| `lineage`, `main_lineage`, `mixed_lineage` | TB-Profiler lineage; `mixed_lineage` is `yes` when TB-Profiler reports more than one main lineage                                                                       |
| `drtype`                                   | TB-Profiler resistance type, such as `MDR-TB` or `Pre-XDR-TB`                                                                                                           |
| one column per drug                        | `resistant`, `not_detected`, `uninterpretable` or `no_result` (see below)                                                                                               |
| `resistance_mutations`                     | mutations graded "Assoc w R" in the WHO catalogue, as `drug: gene change (allele frequency)`                                                                            |
| `interim_resistance_mutations`             | mutations graded "Assoc w R - Interim"                                                                                                                                  |
| `uncertain_mutations`                      | other mutations TB-Profiler reports in resistance genes; they do not make a drug resistant                                                                              |
| `low_coverage_genes`                       | resistance genes with a median depth below `--mtb_min_site_depth`                                                                                                       |
| `consensus_n_percent`                      | percent N in the consensus sequence; above `--mtb_max_consensus_n` the sample is not clustered                                                                         |
| `cluster_<t>snp`                           | cluster ID at each threshold in `--mtb_snp_thresholds`, `-` for a sample with no other sample within the threshold, `not_clustered` for a sample that was not clustered |
| `tbprofiler_version`, `tbdb_version`       | TB-Profiler version and tbdb name and commit                                                                                                                            |
| `notes`                                    | missing results and why                                                                                                                                                 |

A drug is `resistant` when TB-Profiler finds a mutation graded "Assoc w R" or "Assoc w R - Interim" for it. Without such a mutation, the drug is `uninterpretable` when a gene that carries at least one such graded mutation for the drug in tbdb (for example `rpoB` for rifampicin, `katG` and `inhA` for isoniazid) has a median depth below `--mtb_min_site_depth`, and `not_detected` otherwise. Other genes linked to the drug, such as genes with only "Uncertain significance" mutations, can appear in `low_coverage_genes` but do not make a drug `uninterpretable`. `no_result` means TB-Profiler gave no result for the sample.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
