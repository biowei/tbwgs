# biowei/tbwgs

[![GitHub Actions CI Status](https://github.com/biowei/tbwgs/actions/workflows/nf-test.yml/badge.svg)](https://github.com/biowei/tbwgs/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/biowei/tbwgs/actions/workflows/linting.yml/badge.svg)](https://github.com/biowei/tbwgs/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/biowei/tbwgs)

## Introduction

**biowei/tbwgs** is a bioinformatics pipeline for _Mycobacterium tuberculosis_ whole-genome sequencing in a public health laboratory. It takes paired-end Illumina reads and does read QC, mapping to H37Rv, variant calling, drug resistance prediction with TB-Profiler and SNP-distance clustering. It gives one summary row per sample, a MultiQC report and the per-step result files.

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

## Usage

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

## Frozen tool and database versions

Drug resistance calls depend on the TB-Profiler version and on the mutation database. Release 1.0.0 is tested with these versions only:

| Component             | Version                                                                                                                                                                                                                                                       |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| TB-Profiler           | 6.7.0 (container set in `modules/local/tbprofiler_profile`)                                                                                                                                                                                                   |
| tbdb (`--mtb_tbdb`)   | `who_v2+` branch at commit [`c4167762`](https://github.com/jodyphelan/tbdb/tree/c4167762ed361c679507e3f4ebe0cd277eae6102), archive `who_v2plus_c4167762.tar.gz` in [biowei/tbwgs-test-datasets](https://github.com/biowei/tbwgs-test-datasets/tree/main/tbdb) |
| H37Rv reference, mask | `genome.fasta` and `mask.bed` from the same tbdb commit                                                                                                                                                                                                       |

The tbdb archive sha256 is `28da4bcd3f0af7f7db6bce62885b50c441289e4240c5d43f008d49e9480c501f`. If you change TB-Profiler or tbdb, you must validate the pipeline again before you report results.

## Validation status

> [!IMPORTANT]
> Release 1.0.0 has passed a functional check only. It is not validated for clinical or public health reporting.

- The nf-test tests use simulated H37Rv reads. They show that every step runs.
- The default QC thresholds were set from a run on 8 public CRyPTIC isolates (see [biowei/tbtamr-validation](https://github.com/biowei/tbtamr-validation)). For rifampicin, isoniazid and ethambutol, 23 of 24 calls agreed with the phenotype. The one difference is the borderline `rpoB` p.His445Leu mutation.
- 8 isolates are not enough for a formal laboratory validation. Before you use the pipeline for reporting, validate it in your laboratory: use more samples with phenotypes, include local samples, run it on your own compute and write a validation report.

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
