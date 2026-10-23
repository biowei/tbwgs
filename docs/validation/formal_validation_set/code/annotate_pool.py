#!/usr/bin/env python3
"""Join the CRyPTIC table with the VCF summary, and keep the pool and the first 8 isolates.

Usage: annotate_pool.py <pool.table.tsv> <vcfsum.tsv> <cryptic_pool_annotated.tsv>
"""
import sys
import pandas as pd

table, vcfsum, out = sys.argv[1:4]
d = pd.read_csv(table, sep="\t", dtype={"site": str, "subj": str})
v = pd.read_csv(vcfsum, sep="\t")
d["vcf"] = d.VCF.str.split("/").str[-1]
d = d.merge(v, on="vcf", how="left")
cols = ["UNIQUEID", "ENA_SAMPLE", "ENA_RUN", "n_runs", "site", "subj", "prev8", "mbases", "instrument",
        "bioproject", "lineage", "main_lineage", "mean_depth", "n_het", "n_pass_snp", "uniqueid", "VCF"]
cols += [c for c in d.columns if c.endswith(("_BINARY_PHENOTYPE", "_MIC", "_PHENOTYPE_QUALITY"))]
d[d.uniqueid.notna() | d.prev8][cols].to_csv(out, sep="\t", index=False)
