#!/usr/bin/env python3
"""Compare tbwgs mtb_summary.tsv with tbtamr results and CRyPTIC phenotypes.

Usage (from the repository root):
    D=docs/validation/tbtamr_comparison
    python $D/code/compare_tbwgs_tbtamr.py $D/inputs/mtb_summary_v1.1.0.tsv \
        $D/inputs/tbtamr_raw_results.csv $D/inputs/sample_manifest.csv \
        $D/tables/per_drug_calls_v1.1.0.tsv

tbwgs calls: resistant -> R, not_detected -> S, uninterpretable -> U.
tbtamr calls: Resistant / Low-level resistant -> R, Susceptible* -> S.
"""
import csv
import sys

summary, tbtamr, manifest, out = sys.argv[1:5]

DRUGS = ["rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "levofloxacin", "moxifloxacin",
         "bedaquiline", "linezolid", "clofazimine", "delamanid", "amikacin", "kanamycin",
         "capreomycin", "streptomycin", "ethionamide"]
PHENO = {"isoniazid": "INH", "rifampicin": "RIF", "ethambutol": "EMB"}

tw = {r["sample"]: r for r in csv.DictReader(open(summary), delimiter="\t")}
ta = {r["seq_id"]: r for r in csv.DictReader(open(tbtamr))}
ph = {r["sample_id"]: r for r in csv.DictReader(open(manifest))}


def tw_call(v):
    return {"resistant": "R", "not_detected": "S", "uninterpretable": "U"}.get(v, v or "NA")


def ta_call(v):
    return "R" if v in ("Resistant", "Low-level resistant") else ("S" if v.startswith("Susceptible") else "NA")


rows = []
for s in sorted(ta):
    for d in DRUGS:
        t = tw.get(s, {})
        a = ta[s]
        p = ph.get(s, {}).get(f"{PHENO[d]}_BINARY_PHENOTYPE", "") if d in PHENO else ""
        rows.append({
            "sample": s, "drug": d,
            "tbwgs": tw_call(t.get(d, "")), "tbwgs_raw": t.get(d, ""),
            "tbtamr": ta_call(a.get(f"{d} - interpretation", "")),
            "tbtamr_raw": a.get(f"{d} - interpretation", ""),
            "tbtamr_mechanisms": a.get(f"{d} - mechanisms", ""),
            "phenotype": p,
        })
for r in rows:
    r["agree_tbtamr"] = "yes" if r["tbwgs"] == r["tbtamr"] else "no"
    r["agree_pheno"] = "" if not r["phenotype"] else ("yes" if r["tbwgs"] == r["phenotype"] else "no")

with open(out, "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

n = len(rows)
agree = sum(r["agree_tbtamr"] == "yes" for r in rows)
print(f"tbwgs vs tbtamr: {agree}/{n} agree")
for d in DRUGS:
    dr = [r for r in rows if r["drug"] == d]
    print(f"  {d}: {sum(r['agree_tbtamr'] == 'yes' for r in dr)}/{len(dr)}")
print("Discordant vs tbtamr:")
for r in rows:
    if r["agree_tbtamr"] == "no":
        print(f"  {r['sample']} {r['drug']}: tbwgs={r['tbwgs_raw']} tbtamr={r['tbtamr_raw']} [{r['tbtamr_mechanisms']}]")
print("vs phenotype:")
for d in PHENO:
    dr = [r for r in rows if r["drug"] == d]
    tp = sum(r["tbwgs"] == "R" and r["phenotype"] == "R" for r in dr)
    tn = sum(r["tbwgs"] == "S" and r["phenotype"] == "S" for r in dr)
    fp = sum(r["tbwgs"] == "R" and r["phenotype"] == "S" for r in dr)
    fn = sum(r["tbwgs"] != "R" and r["phenotype"] == "R" for r in dr)
    u = sum(r["tbwgs"] == "U" for r in dr)
    print(f"  {d}: TP={tp} TN={tn} FP={fp} FN={fn} U={u}")
print("Per-sample QC and lineage:")
for s in sorted(tw):
    t = tw[s]
    print(f"  {s}: mtbc={t['mtbc_fraction']} depth={t['mean_depth']} breadth={t['breadth']} qc={t['coverage_qc']} "
          f"lineage={t['lineage']} mixed={t['mixed_lineage']} lowcov={t['low_coverage_genes']} c5={t.get('cluster_5snp')} c12={t.get('cluster_12snp')}")
