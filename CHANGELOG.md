# biowei/tbwgs: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2026-10-22]

Initial release of biowei/tbwgs, created with the [nf-core](https://nf-co.re/) template. It is tested with TB-Profiler 6.7.0 and tbdb `who_v2+` at commit `c4167762`. This release has passed a functional check, not a formal laboratory validation (see the README).

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
- Clustering subworkflow `CLUSTER_BY_SNP_DISTANCE`: consensus sequences of `PASS` samples and of `--mtb_historical_consensus` are joined into one alignment with H37Rv (local module `MAKE_SNP_ALIGNMENT`), then snp-sites, snp-dists and single-linkage clustering at each threshold in `--mtb_snp_thresholds` (local module `SNP_CLUSTER`, script `bin/snp_cluster.py`). Results are saved in `clustering/`. A second pipeline test lowers the QC thresholds so both test samples are clustered.
- Summary subworkflow `SUMMARISE_MTB_RESULTS`: one row per sample in `summary/mtb_summary.tsv` and in MultiQC, built by a local module `MTB_SUMMARY` (script `bin/mtb_summary.py`). Samples with coverage QC `FAIL` or without a TB-Profiler or clustering result keep their row. A drug is `uninterpretable` when no resistance mutation is found and a gene linked to it has a median depth below `--mtb_min_site_depth`.
- `--mtb_tbdb` is now required. `--mtb_tbdb_name` is removed: the database name is the name of the tbdb directory.
- `test` profile uses the frozen tbdb archive from biowei/tbwgs-test-datasets instead of the raw tbdb reference and mask.
- `--kraken2db` is now required.
- Thresholds set from a run on 8 public CRyPTIC isolates ([biowei/tbtamr-validation](https://github.com/biowei/tbtamr-validation)):
  - A drug is `uninterpretable` only when a low-coverage gene has a mutation graded "Assoc w R" or "Assoc w R - Interim" for that drug in tbdb. `low_coverage_genes` still lists every resistance gene with low coverage.
  - `--mtb_min_mtbc_fraction` defaults to 90. It only flags samples.
  - New `--mtb_max_consensus_n` (default 10): consensus sequences with a higher percentage of N are not clustered. The summary has a new `consensus_n_percent` column.
- README: frozen tool and database versions, samplesheet format and validation status.

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
| `snp-dists`   | 0.8.2       |
| `snp-sites`   | 2.5.1       |
| `tb-profiler` | 6.7.0       |

### `Deprecated`
