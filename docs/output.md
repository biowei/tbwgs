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
