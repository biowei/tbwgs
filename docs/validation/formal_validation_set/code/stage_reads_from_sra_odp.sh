#!/usr/bin/env bash
# Stage paired-end reads from the AWS SRA open-data bucket into your own S3 bucket
# as gzipped FASTQ, because the tbwgs samplesheet accepts only .fastq.gz files.
#
# Run on an EC2 instance or a Batch job in us-east-1 (no egress charge for the
# SRA open-data bucket). Needs: sra-tools >= 3.0, pigz, aws CLI.
#
# Usage:
#   stage_reads_from_sra_odp.sh <run_list.txt> s3://<your-bucket>/<prefix> [threads]
#
# run_list.txt has one run accession per line (the `run_accession` column of
# samples.csv). Output: s3://<bucket>/<prefix>/<RUN>/<RUN>_1.fastq.gz and _2.fastq.gz
set -euo pipefail

runs=$1
dest=${2%/}
threads=${3:-8}
tmp=${TMPDIR:-/tmp}/stage_sra
mkdir -p "$tmp"

while read -r run; do
    [[ -z "$run" || "$run" == \#* ]] && continue
    if aws s3 ls "$dest/$run/${run}_2.fastq.gz" > /dev/null 2>&1; then
        echo "skip $run (already staged)"
        continue
    fi
    echo "stage $run"
    aws s3 cp --no-sign-request "s3://sra-pub-run-odp/sra/$run/$run" "$tmp/$run.sra"
    fasterq-dump --split-files --threads "$threads" --outdir "$tmp" "$tmp/$run.sra"
    rm -f "$tmp/$run.sra"
    # A paired-end run must give exactly _1 and _2. Stop if it does not.
    if [[ ! -s "$tmp/${run}_1.fastq" || ! -s "$tmp/${run}_2.fastq" ]]; then
        echo "ERROR: $run did not give two read files" >&2
        exit 1
    fi
    rm -f "$tmp/$run.fastq"   # unpaired reads, if any
    pigz -p "$threads" "$tmp/${run}_1.fastq" "$tmp/${run}_2.fastq"
    aws s3 cp "$tmp/${run}_1.fastq.gz" "$dest/$run/${run}_1.fastq.gz"
    aws s3 cp "$tmp/${run}_2.fastq.gz" "$dest/$run/${run}_2.fastq.gz"
    rm -f "$tmp/${run}"_[12].fastq.gz
done < "$runs"
