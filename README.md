# biowei/tbwgs

**A _Mycobacterium tuberculosis_ WGS pipeline for public health laboratories, with its validation evidence in the same repository.**

This repository has two parts:

1. **Validation evidence.** Two evaluations on public CRyPTIC isolates with known phenotypes, plus the data, scripts and tables to build every number again.
2. **The pipeline.** A Nextflow (nf-core template) pipeline that turns paired-end Illumina reads into QC flags, lineage, per-drug resistance calls and SNP clusters.

## Validation status

> [!IMPORTANT]
> Release 1.0.0 has passed a functional check, a comparison with tbtAMR on 8 isolates and a 51-isolate run against phenotype, all on public data. It is not validated for clinical or public health reporting. Each laboratory must validate it on its own samples before it reports results.

| Evaluation                  | Data                                                                                     | Main result                                                                                                                       | Full report                                                        |
| --------------------------- | ---------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| Accuracy against phenotype  | 51 CRyPTIC isolates, run on Seqera Platform with AWS Batch (about $1.26)                 | Sensitivity: rifampicin 100% (23/23), isoniazid 96% (26/27), ethambutol 100% (16/16). All 11 QC, lineage and cluster checks pass. | [formal_validation_51.md](docs/validation/formal_validation_51.md) |
| Agreement with another tool | 8 CRyPTIC isolates × 15 drugs, compared with [tbtAMR](https://github.com/MDU-PHL/tbtamr) | 119 of 120 calls agree (99.2%). Every resistant call and its mutation agree.                                                      | [tbtamr_comparison.md](docs/validation/tbtamr_comparison.md)       |
| Functional check            | Simulated H37Rv reads in [nf-test](https://www.nf-test.com)                              | Every step runs. The tests run in GitHub Actions on each code change.                                                             | [tests/](tests/)                                                   |

### 51 CRyPTIC isolates against phenotype

tbwgs 1.0.0 was run on 51 public CRyPTIC isolates with HIGH-quality phenotypes. The set covers pan-susceptible, mono-resistant, RR, MDR, pre-XDR isolates, and isolates resistant to amikacin, bedaquiline or linezolid. It also has edge cases: low coverage, mixed infections, a low-MTBC sample and same-patient pairs.

| Drug         | Sensitivity (95% CI) | Specificity (95% CI) |
| ------------ | -------------------- | -------------------- |
| Rifampicin   | 100% (86–100), 23/23 | 80% (58–92), 16/20   |
| Isoniazid    | 96% (82–99), 26/27   | 88% (66–97), 15/17   |
| Ethambutol   | 100% (81–100), 16/16 | 92% (74–98), 22/24   |
| Levofloxacin | 88% (53–98), 7/8     | 100% (87–100), 25/25 |
| Moxifloxacin | 100% (61–100), 6/6   | 96% (80–99), 23/24   |
| Bedaquiline  | 25% (5–70), 1/4      | 100% (89–100), 31/31 |
| Linezolid    | 29% (8–64), 2/7      | 100% (89–100), 30/30 |

- **QC, lineage and clustering:** all 11 expected checks pass. Low-coverage isolates fail coverage QC, both suspected mixed infections are flagged, same-patient pairs cluster at the expected SNP threshold, and the lineage matches CRyPTIC for 49 of 49 isolates.
- **False positives** are mostly borderline `rpoB`, `inhA` or `embB` mutations, or a mixed infection.
- **New and repurposed drugs:** the missed isolates carry no WHO catalogue mutation. The gap is in the catalogue, not in a pipeline step.
- All 12 drugs, every disagreement with its likely cause, and the scripts are in the [full report](docs/validation/formal_validation_51.md).

![tbwgs 1.0.0 on 51 CRyPTIC isolates](docs/validation/formal_validation_51/formal_validation_51.png)

### 8 CRyPTIC isolates compared with tbtAMR

- **119 of 120** sample × drug calls agree with tbtAMR. The one difference is linezolid in ERR4813225. tbwgs calls it `uninterpretable` because the `rrl` gene has a median depth of 7×.
- **23 of 24** calls agree with the phenotype for rifampicin, isoniazid and ethambutol, for both tools. The one difference is the borderline `rpoB` p.His445Leu mutation.
- This comparison set the default QC thresholds and changed the rule for `uninterpretable` calls ([#9](https://github.com/biowei/tbwgs/pull/9)).

![tbwgs vs tbtamr on 8 CRyPTIC isolates](docs/validation/tbtamr_comparison/tbwgs_vs_tbtamr.png)

### Open findings

1. **LOW_MTBC flag:** 18 of 51 isolates get the flag, but only 1 was expected to. The other 17 have 82–90% MTBC reads, just below the default `--mtb_min_mtbc_fraction 90`. The flag does not stop any step, so their calls are not affected. The default stays at 90 until the Kraken2 reports are checked.
2. **Read integrity:** reads fetched over https were corrupted twice without an error. For runs that matter, stage the reads with an MD5 check first.
3. **Interim mutations** count as resistant. This gives one rifampicin false positive.

### How the evidence is kept

- **Frozen versions:** TB-Profiler 6.7.0 and one tbdb commit, with its sha256 (see [Frozen tool and database versions](#frozen-tool-and-database-versions)). A change to either needs a new validation.
- **Reproducible numbers:** each report has its inputs, the truth tables, the scripts and the output tables. One command builds the tables and the figure again.
- **Related repositories:** [biowei/tbtamr-validation](https://github.com/biowei/tbtamr-validation) has the tbtAMR results and raw outputs. [biowei/tbwgs-test-datasets](https://github.com/biowei/tbwgs-test-datasets) has the frozen tbdb archive and the test data.

### Before you report results

These runs use public data only. Before you use the pipeline for reporting, validate it in your laboratory: use more samples with phenotypes, include local samples, run it on your own compute and write a validation report.

## What the pipeline does

tbwgs takes paired-end Illumina reads and does read QC, mapping to H37Rv, variant calling, drug resistance prediction with TB-Profiler and SNP-distance clustering. It gives one summary row per sample, a MultiQC report and the per-step result files.

1. Merge lanes of the same sample ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))
4. Fraction of _M. tuberculosis_ complex reads ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
5. Prepare the H37Rv reference and BWA index from the frozen tbdb ([`SAMtools`](https://www.htslib.org/), [`BWA`](https://github.com/lh3/bwa))
6. Map reads to H37Rv, sort and mark duplicates ([`BWA-MEM`](https://github.com/lh3/bwa), [`SAMtools`](https://www.htslib.org/))
7. Coverage QC flag per sample, without dropping samples ([`mosdepth`](https://github.com/brentp/mosdepth), [`SAMtools`](https://www.htslib.org/))
8. Haploid variant calls and masked consensus sequence ([`BCFtools`](https://samtools.github.io/bcftools/))
9. Drug resistance and lineage from the trimmed reads ([`TB-Profiler`](https://github.com/jodyphelan/TBProfiler))
10. Pairwise SNP distances and single-linkage clusters of `PASS` samples and historical consensus sequences ([`snp-sites`](https://github.com/sanger-pathogens/snp-sites), [`snp-dists`](https://github.com/tseemann/snp-dists))
11. One summary row per sample with coverage QC, lineage, per-drug calls (resistant, not detected or uninterpretable) and clusters
12. Present QC ([`MultiQC`](http://multiqc.info/))

## How to run it

[![GitHub Actions CI Status](https://github.com/biowei/tbwgs/actions/workflows/nf-test.yml/badge.svg)](https://github.com/biowei/tbwgs/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/biowei/tbwgs/actions/workflows/linting.yml/badge.svg)](https://github.com/biowei/tbwgs/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/biowei/tbwgs)

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
TB001,TB001_S1_L001_R1_001.fastq.gz,TB001_S1_L001_R2_001.fastq.gz
TB001,TB001_S1_L002_R1_001.fastq.gz,TB001_S1_L002_R2_001.fastq.gz
```

Each row is one pair of FASTQ files. Only paired-end reads are supported. Rows with the same `sample` are merged as lanes of one sample.

Now, you can run the pipeline using:

```bash
nextflow run biowei/tbwgs \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --mtb_tbdb <TBDB_ARCHIVE> \
   --kraken2db <KRAKEN2DB> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

### Frozen tool and database versions

Drug resistance calls depend on the TB-Profiler version and on the mutation database. Release 1.0.0 is tested with these versions only:

| Component             | Version                                                                                                                                                                                                                                                       |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| TB-Profiler           | 6.7.0 (container set in `modules/local/tbprofiler_profile`)                                                                                                                                                                                                   |
| tbdb (`--mtb_tbdb`)   | `who_v2+` branch at commit [`c4167762`](https://github.com/jodyphelan/tbdb/tree/c4167762ed361c679507e3f4ebe0cd277eae6102), archive `who_v2plus_c4167762.tar.gz` in [biowei/tbwgs-test-datasets](https://github.com/biowei/tbwgs-test-datasets/tree/main/tbdb) |
| H37Rv reference, mask | `genome.fasta` and `mask.bed` from the same tbdb commit                                                                                                                                                                                                       |

The tbdb archive sha256 is `28da4bcd3f0af7f7db6bce62885b50c441289e4240c5d43f008d49e9480c501f`. If you change TB-Profiler or tbdb, you must validate the pipeline again before you report results.

## Credits

biowei/tbwgs was originally written by Chris.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use biowei/tbwgs for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
