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
- Coverage subworkflow `CHECK_MAPPING_COVERAGE`: mosdepth and samtools stats. Each sample gets `mtb_mean_depth`, `mtb_breadth` and `mtb_qc` (`PASS` / `FAIL`) in its meta from `--mtb_min_mean_depth`, `--mtb_min_breadth` and `--mtb_min_site_depth`. No sample is dropped. Reports and a LOW / PASS coverage BED are saved in `coverage/`.
- Variant subworkflow `CALL_MTB_VARIANTS`: `bcftools mpileup` and haploid `call`, `norm`, `filter` (`QUAL >= 20`, depth at least `--mtb_min_site_depth`) and `stats` for every sample, then a SNP-only `bcftools consensus` masked with N over the fixed mask and the sample's low-coverage sites. A local module `MERGE_MASK_BED` builds the per-sample mask. VCFs are saved in `variants/` and consensus sequences in `consensus/`.
- Drug resistance subworkflow `PROFILE_DRUG_RESISTANCE` with a local `TBPROFILER_PROFILE` module (TB-Profiler 6.7.0, adapted from the nf-core 3.0.8 module) that reads the frozen tbdb given by `--mtb_tbdb`. Results are saved in `tbprofiler/`.
- `--mtb_tbdb` is now required. `--mtb_tbdb_name` is removed: the database name is the name of the tbdb directory.
- `test` profile uses the frozen tbdb archive from biowei/tbwgs-test-datasets instead of the raw tbdb reference and mask.
- `--kraken2db` is now required.

### `Fixed`

### `Dependencies`

| Dependency    | New version |
| ------------- | ----------- |
| `bcftools`    | 1.23.1      |
| `bwa`         | 0.7.19      |
| `fastp`       | 1.4.0       |
| `kraken2`     | 2.1.6       |
| `mosdepth`    | 0.3.14      |
| `samtools`    | 1.24        |
| `tb-profiler` | 6.7.0       |

### `Deprecated`
