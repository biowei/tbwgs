#!/usr/bin/env python3
"""Build the candidate pool from the CRyPTIC reuse table and the AWS SRA metadata.

Usage:
  build_pool.py <CRyPTIC_reuse_table_20240917.csv> <mtb_sra.parquet> <first8_manifest.csv> <out_prefix>

Writes:
  <out_prefix>.table.tsv  every CRyPTIC isolate with run count, site, subject, SRA Mbases
  <out_prefix>.vcf_urls   VCF URLs for the candidate pool (download, then run
                          summarise_cryptic_vcfs.py on them)

Pool: for each stratum, up to 8 isolates per site (one run, >= 30x by SRA Mbases,
not in the first 8), plus every isolate from a patient with more than one isolate,
plus all isolates with 10-22x by SRA Mbases.
"""
import sys
import numpy as np
import pandas as pd

reuse_csv, sra_parquet, first8, out = sys.argv[1:5]
VCF_BASE = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/"
GENOME_MB = 4.41

d = pd.read_csv(reuse_csv)
s = pd.read_parquet(sra_parquet)[["acc", "mbases", "instrument", "bioproject"]].drop_duplicates("acc")
d["site"] = d.UNIQUEID.str.split(".").str[1]
d["subj"] = d.UNIQUEID.str.split(".").str[3]
d["n_runs"] = d.ENA_RUN.str.count(r"\.") + 1
d = d.merge(s, left_on="ENA_RUN", right_on="acc", how="left").drop(columns="acc")
d["prev8"] = d.ENA_RUN.isin(set(pd.read_csv(first8).run_accession))

drugs = ["AMI", "BDQ", "CFZ", "DLM", "EMB", "ETH", "INH", "KAN", "LEV", "LZD", "MXF", "RIF", "RFB"]
P = {x: d[f"{x}_BINARY_PHENOTYPE"].where(d[f"{x}_PHENOTYPE_QUALITY"] == "HIGH") for x in drugs}
base = (d.n_runs == 1) & d.mbases.notna() & (d.mbases >= GENOME_MB * 30) & ~d.prev8
strata = {
    "PANS": base & (P["INH"] == "S") & (P["RIF"] == "S") & (P["EMB"] == "S") & (P["LEV"] == "S") & (P["MXF"] == "S")
    & np.logical_and.reduce([d[f"{x}_BINARY_PHENOTYPE"].eq("S") for x in drugs]),
    "HR": base & (P["INH"] == "R") & (P["RIF"] == "S") & (P["RFB"] != "R") & (P["LEV"] == "S"),
    "RR": base & (P["RIF"] == "R") & (P["INH"] == "S"),
    "MDR": base & (P["INH"] == "R") & (P["RIF"] == "R") & (P["LEV"] == "S") & (P["MXF"] == "S"),
    "PREXDR": base & (P["INH"] == "R") & (P["RIF"] == "R") & ((P["LEV"] == "R") | (P["MXF"] == "R")),
    "BDQ_R": base & (P["BDQ"] == "R"),
    "LZD_R": base & (P["LZD"] == "R"),
    "AMI_R": base & (P["AMI"] == "R") & (P["KAN"] == "R"),
}
pool = set()
for k, m in strata.items():
    pick = d[m].groupby("site", group_keys=False).apply(lambda g: g.sample(min(len(g), 8), random_state=1))
    pool |= set(pick.UNIQUEID)
multi = d.groupby(["site", "subj"]).UNIQUEID.transform("count") > 1
pool |= set(d[multi & (d.n_runs == 1)].UNIQUEID)
low = d[(d.n_runs == 1) & d.mbases.between(GENOME_MB * 10, GENOME_MB * 22) & P["INH"].notna() & P["RIF"].notna()]
pool |= set(low.UNIQUEID)

d.to_csv(out + ".table.tsv", sep="\t", index=False)
urls = d[d.UNIQUEID.isin(pool)].VCF.str.replace("../reproducibility/", VCF_BASE, regex=False)
urls.to_csv(out + ".vcf_urls", index=False, header=False)
print("pool", len(pool))
