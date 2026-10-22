#!/bin/bash
# How tbwgs was run on the 8 validation isolates (cloud container, 4 CPU, 15 GB RAM).
# ENA was not reachable from the container, so reads came from the NCBI SRA open-data mirror.
set -euo pipefail

SAMPLES="ERR4813225 ERR8975345 ERR3287775 ERR4809348 ERR4831755 ERR3287635 ERR8975009 ERR4811536"

# 1. Reads: SRA open-data bucket -> paired FASTQ (read counts match data/sample_manifest.csv)
for s in $SAMPLES; do
    curl -sSL -o "$s.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$s/$s"
    fasterq-dump --split-files "./$s.sra" -O .
    pigz "${s}_1.fastq" "${s}_2.fastq"
done

# 2. Kraken2 database for the MTBC fraction
curl -sSL https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz | tar xz -C k2_standard_08gb

# 3. Samplesheet
{
    echo "sample,fastq_1,fastq_2"
    for s in $SAMPLES; do echo "$s,$PWD/${s}_1.fastq.gz,$PWD/${s}_2.fastq.gz"; done
} > samplesheet.csv

# 4. tbwgs (biowei/tbwgs dev @ 51ded8c), conda profile, default thresholds
nextflow run biowei/tbwgs -r 51ded8c -profile conda \
    --input samplesheet.csv \
    --mtb_tbdb who_v2plus_c4167762.tar.gz \
    --kraken2db k2_standard_08gb \
    --outdir tbwgs_results

# 5. Summary with the rule merged in biowei/tbwgs PR #9 (mtb_summary.py 1.1.0),
#    from the same TB-Profiler results; samples_qc.tsv adds the consensus N percent
mtb_summary.py --samples samples_qc.tsv --genes tbdb/genes.bed --mutations tbdb/mutations.json \
    --tbprofiler tbwgs_results/tbprofiler/*.results.json \
    --clusters tbwgs_results/clustering/snp.clusters.tsv \
    --min-mtbc-fraction 90 --max-consensus-n 10 --prefix mtb_summary_v1.1.0
