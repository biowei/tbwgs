#!/usr/bin/env python3
"""Compare tbwgs calls on the 51-isolate formal validation set with CRyPTIC phenotypes and expected QC results.

Run from the repository root:
    python docs/validation/formal_validation_51/code/compare_formal_validation.py

Writes three tables to docs/validation/formal_validation_51/tables/:
    per_sample_drug_calls.tsv  one row per sample x drug, with the outcome (TP, FN, FP, TN, no_phenotype, qc_fail)
    per_drug_performance.tsv   sensitivity and specificity per drug, with 95% Wilson intervals
    expected_checks.tsv        each edge case and same-patient pair, with the expected and observed result
"""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRUGS = ["rifampicin", "isoniazid", "ethambutol", "levofloxacin", "moxifloxacin", "ethionamide",
         "kanamycin", "amikacin", "delamanid", "linezolid", "bedaquiline", "clofazimine"]


def read(path):
    with open(path, newline="") as handle:
        return {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}


def write(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def wilson(k, n, z=1.96):
    if n == 0:
        return "", ""
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return f"{max(0, centre - half):.3f}", f"{min(1, centre + half):.3f}"


summary = read(ROOT / "inputs/mtb_summary.tsv")
truth = read(ROOT / "inputs/truth_phenotypes.tsv")
expected = read(ROOT / "inputs/expected_qc_and_clusters.tsv")


def outcome(sample, drug):
    phenotype, call = truth[sample][drug], summary[sample][drug]
    if summary[sample]["coverage_qc"] != "PASS":
        return "qc_fail"
    if phenotype not in ("R", "S"):
        return "no_phenotype"
    if call == "uninterpretable":
        return "uninterpretable"
    resistant = call == "resistant"
    if phenotype == "R":
        return "TP" if resistant else "FN"
    return "FP" if resistant else "TN"


calls = []
for sample in sorted(summary):
    for drug in DRUGS:
        calls.append([sample, expected[sample]["stratum"], drug, truth[sample][drug],
                      truth[sample][drug + "_mic"], summary[sample][drug], outcome(sample, drug)])
write(ROOT / "tables/per_sample_drug_calls.tsv",
      ["sample", "stratum", "drug", "phenotype", "mic", "tbwgs_call", "outcome"], calls)

performance = []
for drug in DRUGS:
    counts = {k: sum(1 for c in calls if c[2] == drug and c[6] == k)
              for k in ("TP", "FN", "FP", "TN", "uninterpretable", "no_phenotype")}
    sens = counts["TP"] / (counts["TP"] + counts["FN"]) if counts["TP"] + counts["FN"] else None
    spec = counts["TN"] / (counts["TN"] + counts["FP"]) if counts["TN"] + counts["FP"] else None
    performance.append([drug, counts["TP"], counts["FN"], counts["FP"], counts["TN"],
                        counts["uninterpretable"], counts["no_phenotype"],
                        "" if sens is None else f"{sens:.3f}", *wilson(counts["TP"], counts["TP"] + counts["FN"]),
                        "" if spec is None else f"{spec:.3f}", *wilson(counts["TN"], counts["TN"] + counts["FP"])])
write(ROOT / "tables/per_drug_performance.tsv",
      ["drug", "TP", "FN", "FP", "TN", "uninterpretable", "no_phenotype",
       "sensitivity", "sens_low", "sens_high", "specificity", "spec_low", "spec_high"], performance)


def group(stratum):
    return [s for s in sorted(expected) if expected[s]["stratum"] == stratum]


def pairs():
    by_id = {}
    for s, r in expected.items():
        if r["pair_id"]:
            by_id.setdefault(r["pair_id"], []).append(s)
    return sorted(by_id.items())


checks = []
for s in group("low_coverage"):
    obs = summary[s]["coverage_qc"]
    checks.append(["low_coverage", s, "coverage_qc FAIL", f"{obs} ({summary[s]['mean_depth']}x)", obs == "FAIL"])
for s in group("suspected_mixed"):
    obs = summary[s]["mixed_lineage"]
    checks.append(["suspected_mixed", s, "mixed_lineage yes", f"{obs} ({summary[s]['lineage']})", obs == "yes"])
for s in group("overlap_low_mtbc"):
    obs = summary[s]["mtbc_flag"]
    checks.append(["overlap_low_mtbc", s, "mtbc_flag LOW_MTBC", f"{obs} ({summary[s]['mtbc_fraction']}%)",
                   obs == "LOW_MTBC"])
for s in group("overlap_rpoB_H445L"):
    obs = summary[s]["rifampicin"]
    checks.append(["overlap_rpoB_H445L", s, "rifampicin resistant", obs, obs == "resistant"])
for pair_id, (a, b) in pairs():
    level = "cluster_12snp" if pair_id.startswith("pair_6to12") else "cluster_5snp"
    ca, cb = summary[a][level], summary[b][level]
    if pair_id.startswith("pair_gt100"):
        ok = ca in ("-", "not_clustered") or ca != cb
        checks.append([pair_id, f"{a}/{b}", "not clustered", f"{ca} / {cb}", ok])
    else:
        ok = ca not in ("-", "not_clustered") and ca == cb
        checks.append([pair_id, f"{a}/{b}", f"same {level}", f"{ca} / {cb}", ok])
def same_lineage(called, cryptic):
    # CRyPTIC reports a less specific lineage (lineage1) than TB-Profiler (lineage1.2.2.2); a mixed call has two
    return any(c == cryptic or c.startswith(cryptic + ".") for c in called.split(";"))


lineage_ok = sum(1 for s, r in expected.items() if r["cryptic_lineage"]
                 and same_lineage(summary[s]["lineage"], r["cryptic_lineage"]))
lineage_n = sum(1 for r in expected.values() if r["cryptic_lineage"])
checks.append(["lineage", f"{lineage_n} samples", "same lineage as CRyPTIC", f"{lineage_ok}/{lineage_n}",
               lineage_ok == lineage_n])
write(ROOT / "tables/expected_checks.tsv", ["check", "samples", "expected", "observed", "pass"],
      [[*c[:4], "yes" if c[4] else "no"] for c in checks])

for row in performance:
    print("\t".join(map(str, row)))
print(f"expected checks passed: {sum(c[4] for c in checks)}/{len(checks)}")
