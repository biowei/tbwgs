# biowei/tbwgs: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of biowei/tbwgs, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- MTB parameter groups in `nextflow_schema.json`: reference and database, QC thresholds, SNP clustering.
- Paired-end only samplesheet schema (`sample`, `fastq_1`, `fastq_2`).
- `test` profile uses reads simulated from H37Rv, hosted in [biowei/tbwgs-test-datasets](https://github.com/biowei/tbwgs-test-datasets).
- Read preprocessing subworkflow `PREPROCESS_READS`: lane merging (`cat`), FastQC, fastp and Kraken2, with the percentage of MTBC reads per sample.
- Reference subworkflow `PREPARE_MTB_REFERENCE`: H37Rv FASTA and mask from the frozen tbdb (or `--mtb_reference` and `--mtb_mask_bed`), samtools faidx and BWA index.
- Mapping subworkflow `MAP_READS_TO_H37RV`: BWA-MEM with a read group per sample, `samtools` collate, fixmate, sort and markdup in one step, and a BAM index. BAMs are saved in `alignment/`.
- `--kraken2db` is now required. The pipeline stops early unless `--mtb_tbdb`, or both `--mtb_reference` and `--mtb_mask_bed`, are given.

### `Fixed`

### `Dependencies`

| Dependency | New version |
| ---------- | ----------- |
| `bwa`      | 0.7.19      |
| `fastp`    | 1.4.0       |
| `kraken2`  | 2.1.6       |
| `samtools` | 1.24        |

### `Deprecated`
