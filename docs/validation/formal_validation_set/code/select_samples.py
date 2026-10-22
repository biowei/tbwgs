#!/usr/bin/env python3
"""Select the tbwgs formal validation set from the annotated CRyPTIC pool.

Input:  tables/cryptic_pool_annotated.tsv  (one row per CRyPTIC isolate in the candidate pool)
        tables/same_patient_snp_distance.tsv
Output: tables/selected_samples.tsv        (one row per selected isolate, with its stratum)

The selection is deterministic (fixed seed). See README.md for the rules.
"""
import sys
import numpy as np
import pandas as pd

pool_tsv, pairs_tsv, out_tsv = sys.argv[1:4]
SEED = 20261022
d = pd.read_csv(pool_tsv, sep="\t", dtype={"site": str})
pairs = pd.read_csv(pairs_tsv, sep="\t")

DRUGS = ["AMI", "BDQ", "CFZ", "DLM", "EMB", "ETH", "INH", "KAN", "LEV", "LZD", "MXF", "RIF", "RFB"]
hq = {x: d[f"{x}_BINARY_PHENOTYPE"].where(d[f"{x}_PHENOTYPE_QUALITY"] == "HIGH") for x in DRUGS}

# Core samples: one run, at least 40x in the CRyPTIC VCF, no sign of mixture, lineage called,
# and HIGH quality INH, RIF and EMB phenotypes.
core = (
    (d.n_runs == 1) & (d.mean_depth >= 40) & (d.n_het <= 150) & d.lineage.notna()
    & hq["INH"].notna() & hq["RIF"].notna() & hq["EMB"].isin(["R", "S"])
    & ~d.prev8
)
all_s = np.logical_and.reduce([d[f"{x}_BINARY_PHENOTYPE"].eq("S") for x in DRUGS])
STRATA = [
    # name, size, rule
    ("pan_susceptible", 8, all_s & (hq["INH"] == "S") & (hq["RIF"] == "S") & (hq["LEV"] == "S")),
    ("INH_mono", 5, (hq["INH"] == "R") & (hq["RIF"] == "S") & (hq["RFB"] != "R") & (hq["LEV"] == "S")),
    ("RIF_R_INH_S", 4, (hq["RIF"] == "R") & (hq["INH"] == "S")),
    ("MDR_FQ_S", 6, (hq["INH"] == "R") & (hq["RIF"] == "R") & (hq["LEV"] == "S") & (hq["MXF"] == "S")),
    ("MDR_FQ_R", 6, (hq["INH"] == "R") & (hq["RIF"] == "R") & ((hq["LEV"] == "R") | (hq["MXF"] == "R"))),
    ("BDQ_R", 3, hq["BDQ"] == "R"),
    ("LZD_R", 3, hq["LZD"] == "R"),
    ("AMI_R", 2, (hq["AMI"] == "R") & (hq["KAN"] == "R")),
]

chosen = []
taken = set()


def pick(cands, n, stratum):
    """Pick n rows, spreading over main lineage first and then site (round robin)."""
    cands = cands[~cands.UNIQUEID.isin(taken)].sample(frac=1, random_state=SEED)
    groups = {k: g for k, g in cands.groupby("main_lineage")}
    out, used_sites = [], set()
    while len(out) < n and any(len(g) for g in groups.values()):
        for lin in sorted(groups):
            g = groups[lin]
            if not len(g) or len(out) >= n:
                continue
            # prefer a site not used yet in this stratum
            fresh = g[~g.site.isin(used_sites)]
            row = (fresh if len(fresh) else g).iloc[0]
            groups[lin] = g[g.UNIQUEID != row.UNIQUEID]
            used_sites.add(row.site)
            out.append(row.UNIQUEID)
    for u in out:
        taken.add(u)
        chosen.append({"UNIQUEID": u, "stratum": stratum})
    return out


for name, n, rule in STRATA:
    pick(d[core & rule], n, name)

# Same-patient pairs: known close pairs for SNP clustering, one 6-12 SNP pair and one
# pair of different strains from the same patient.
good = set(d[(d.n_runs == 1) & (d.mean_depth >= 40) & (d.n_het <= 150) & d.lineage.notna() & ~d.prev8].uniqueid)
p = pairs[pairs.id_a.isin(good) & pairs.id_b.isin(good)].sample(frac=1, random_state=SEED)
id2u = dict(zip(d.uniqueid, d.UNIQUEID))
for label, rule, n in [
    ("pair_le5", p.snp_distance <= 2, 2),
    ("pair_6to12", p.snp_distance.between(7, 11), 1),
    ("pair_gt100", p.snp_distance > 100, 1),
]:
    k = 0
    for _, r in p[rule].iterrows():
        a, b = id2u[r.id_a], id2u[r.id_b]
        if a in taken or b in taken:
            continue
        for u in (a, b):
            taken.add(u)
            chosen.append({"UNIQUEID": u, "stratum": label, "pair_id": f"{label}_{k + 1}",
                           "cryptic_snp_distance": r.snp_distance})
        k += 1
        if k == n:
            break

# Edge cases chosen by rule
edge = d[(d.n_runs == 1) & ~d.UNIQUEID.isin(taken) & hq["INH"].notna() & hq["RIF"].notna()]
low = edge[(edge.n_het <= 150) & edge.lineage.notna()]
for lo, hi in [(9, 12), (14, 19)]:  # below the tbwgs 20x mean depth threshold
    r = low[low.mean_depth.between(lo, hi)].sort_values("mean_depth").iloc[0]
    taken.add(r.UNIQUEID)
    chosen.append({"UNIQUEID": r.UNIQUEID, "stratum": "low_coverage"})
mixed = edge[(edge.mean_depth >= 60) & (edge.n_het >= 1500)].sort_values("n_het", ascending=False)
for _, r in mixed.head(2).iterrows():
    taken.add(r.UNIQUEID)
    chosen.append({"UNIQUEID": r.UNIQUEID, "stratum": "suspected_mixed"})

# Overlap with the first 8-isolate check (regression cases)
for run, label in [("ERR4811536", "overlap_rpoB_H445L"), ("ERR4831755", "overlap_low_mtbc")]:
    chosen.append({"UNIQUEID": d.loc[d.ENA_RUN == run, "UNIQUEID"].iloc[0], "stratum": label})

sel = pd.DataFrame(chosen).merge(d, on="UNIQUEID", how="left")
sel.to_csv(out_tsv, sep="\t", index=False)
print(sel.stratum.value_counts().to_string())
print("total", len(sel))
