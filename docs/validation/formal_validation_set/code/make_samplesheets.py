#!/usr/bin/env python3
"""Write the tbwgs samplesheets, the phenotype truth table and the expected QC/cluster table.

Usage: make_samplesheets.py <selected_samples.tsv> <ena_filereport.tsv> <outdir> [s3_prefix]
"""
import os, sys
import pandas as pd

sel_tsv, ena_tsv, outdir = sys.argv[1:4]
s3 = (sys.argv[4] if len(sys.argv) > 4 else "s3://YOUR-BUCKET/tbwgs-validation/reads").rstrip("/")
sel = pd.read_csv(sel_tsv, sep="\t")
ena = pd.read_csv(ena_tsv, sep="\t").set_index("run_accession")

# tbwgs mtb_summary.tsv column for each CRyPTIC drug (rifabutin has no tbwgs column)
DRUG = {"RIF": "rifampicin", "INH": "isoniazid", "EMB": "ethambutol", "LEV": "levofloxacin",
        "MXF": "moxifloxacin", "BDQ": "bedaquiline", "LZD": "linezolid", "CFZ": "clofazimine",
        "DLM": "delamanid", "AMI": "amikacin", "KAN": "kanamycin", "ETH": "ethionamide", "RFB": "rifabutin"}

ss_ena, ss_s3, truth, exp = [], [], [], []
for _, r in sel.iterrows():
    run = r.ENA_RUN
    f1, f2 = sorted(ena.loc[run, "fastq_ftp"].split(";"))
    assert f1.endswith("_1.fastq.gz") and f2.endswith("_2.fastq.gz"), run
    ss_ena.append({"sample": run, "fastq_1": "https://" + f1, "fastq_2": "https://" + f2})
    ss_s3.append({"sample": run, "fastq_1": f"{s3}/{run}/{run}_1.fastq.gz", "fastq_2": f"{s3}/{run}/{run}_2.fastq.gz"})

    t = {"sample": run, "ena_sample": r.ENA_SAMPLE, "cryptic_uniqueid": r.UNIQUEID, "stratum": r.stratum,
         "site": f"{int(r.site):02d}"}
    for code, name in DRUG.items():
        q = r[f"{code}_PHENOTYPE_QUALITY"]
        b = r[f"{code}_BINARY_PHENOTYPE"]
        # Truth is the binary phenotype only when CRyPTIC rates the MIC HIGH quality.
        t[name] = b if q == "HIGH" and b in ("R", "S") else ""
        t[f"{name}_mic"] = "" if pd.isna(r[f"{code}_MIC"]) else r[f"{code}_MIC"]
        t[f"{name}_quality"] = "" if pd.isna(q) else q
    truth.append(t)

    e = {"sample": run, "stratum": r.stratum,
         "cryptic_lineage": "" if pd.isna(r.lineage) else r.lineage,
         "cryptic_vcf_depth": "" if pd.isna(r.mean_depth) else r.mean_depth,
         "cryptic_vcf_low_frs_sites": "" if pd.isna(r.n_het) else int(r.n_het),
         "pair_id": "" if pd.isna(r.get("pair_id")) else r.pair_id,
         "cryptic_snp_distance": "" if pd.isna(r.get("cryptic_snp_distance")) else int(r.cryptic_snp_distance)}
    st = r.stratum
    if st == "low_coverage":
        e["expected"] = "coverage_qc fail (mean depth below 20x); calls in low-depth genes uninterpretable"
    elif st == "suspected_mixed":
        e["expected"] = "mixed_lineage flagged (two lineages at intermediate read fractions in the CRyPTIC VCF)"
    elif st == "pair_le5":
        e["expected"] = "same cluster_5snp as its pair (same patient, CRyPTIC distance <= 2)"
    elif st == "pair_6to12":
        e["expected"] = "same cluster_12snp as its pair; cluster_5snp may differ (CRyPTIC distance 11)"
    elif st == "pair_gt100":
        e["expected"] = "not clustered with its pair (same patient, different lineage)"
    elif st == "overlap_rpoB_H445L":
        e["expected"] = "rifampicin resistant (rpoB p.His445Leu) vs phenotype S at MIC 0.5: known discordance"
    elif st == "overlap_low_mtbc":
        e["expected"] = "mtbc_flag set (MTBC fraction about 75% in the first run); calls as in the first run"
    else:
        e["expected"] = "QC pass; compare drug calls with truth"
    exp.append(e)

os.makedirs(outdir, exist_ok=True)
pd.DataFrame(ss_ena).to_csv(f"{outdir}/samplesheet_ena.csv", index=False)
pd.DataFrame(ss_s3).to_csv(f"{outdir}/samplesheet_s3.csv", index=False)
pd.DataFrame(truth).to_csv(f"{outdir}/truth_phenotypes.tsv", sep="\t", index=False)
pd.DataFrame(exp).to_csv(f"{outdir}/expected_qc_and_clusters.tsv", sep="\t", index=False)
with open(f"{outdir}/run_list.txt", "w") as fh:
    fh.write("\n".join(sel.ENA_RUN) + "\n")
print(len(ss_ena), "samples")
