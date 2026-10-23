#!/usr/bin/env python3
"""Figure: tbwgs 1.0.0 on the 51-isolate formal validation set.

A: sensitivity and specificity per drug against CRyPTIC phenotypes, with 95% Wilson intervals.
B: MTBC fraction and mean depth per sample, with the QC thresholds.
C: outcome of every sample x drug call.

Run from the repository root, after compare_formal_validation.py:
    python docs/validation/formal_validation_51/code/plot_formal_validation.py
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "formal_validation_51.png"

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, YELLOW, NEUTRAL, CRITICAL = "#2a78d6", "#eb6834", "#eda100", "#f0efec", "#d03b3b"
PALE_BLUE = "#b9d3f2"

ABBR = {"rifampicin": "RIF", "isoniazid": "INH", "ethambutol": "EMB", "levofloxacin": "LFX",
        "moxifloxacin": "MFX", "ethionamide": "ETO", "kanamycin": "KAN", "amikacin": "AMK",
        "delamanid": "DLM", "linezolid": "LZD", "bedaquiline": "BDQ", "clofazimine": "CFZ"}
STRATA = ["pan_susceptible", "INH_mono", "RIF_R_INH_S", "MDR_FQ_S", "MDR_FQ_R", "AMI_R", "BDQ_R", "LZD_R",
          "pair_le5", "pair_6to12", "pair_gt100", "suspected_mixed", "low_coverage",
          "overlap_low_mtbc", "overlap_rpoB_H445L"]
SHORT = {"pair_le5": "pair <=5", "pair_6to12": "pair 6-12", "pair_gt100": "pair >100",
         "suspected_mixed": "mixed", "low_coverage": "low cov", "overlap_low_mtbc": "MTBC",
         "overlap_rpoB_H445L": "H445L", "RIF_R_INH_S": "RIF R INH S"}
OUTCOMES = [("TP", "Resistant, called resistant", BLUE),
            ("TN", "Susceptible, called not detected", PALE_BLUE),
            ("FN", "Resistant, missed", CRITICAL),
            ("FP", "Susceptible, called resistant", YELLOW),
            ("no_phenotype", "No HIGH-quality phenotype", NEUTRAL),
            ("qc_fail", "Coverage QC failed, not scored", "#c9c7c1")]


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


perf = read(ROOT / "tables/per_drug_performance.tsv")
calls = read(ROOT / "tables/per_sample_drug_calls.tsv")
summary = {r["sample"]: r for r in read(ROOT / "inputs/mtb_summary.tsv")}
expected = {r["sample"]: r for r in read(ROOT / "inputs/expected_qc_and_clusters.tsv")}
drugs = [r["drug"] for r in perf]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                     "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2})
fig = plt.figure(figsize=(13, 10.8), facecolor=SURFACE)
grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1.25, 1],
                        hspace=0.5, wspace=0.28, left=0.07, right=0.98, top=0.86, bottom=0.14)
fig.text(0.07, 0.965, "tbwgs 1.0.0 on 51 CRyPTIC isolates", fontsize=14, weight="bold")
fig.text(0.07, 0.938, "First-line drugs and fluoroquinolones reach 88-100% sensitivity. Bedaquiline, linezolid and clofazimine "
         "miss most phenotypic resistance.\nAll 11 expected QC and cluster checks pass. 18 samples fall below the 90% MTBC flag, but only 1 was expected to.",
         fontsize=9.5, color=INK2)

# A: forest plot of sensitivity and specificity
ax = fig.add_subplot(grid[0, 0], facecolor=SURFACE)
y = list(range(len(drugs)))[::-1]
for metric, low, high, colour, offset, label in (
        ("sensitivity", "sens_low", "sens_high", BLUE, 0.17, "Sensitivity"),
        ("specificity", "spec_low", "spec_high", ORANGE, -0.17, "Specificity")):
    for yi, r in zip(y, perf):
        if not r[metric]:
            continue
        v, lo, hi = float(r[metric]) * 100, float(r[low]) * 100, float(r[high]) * 100
        ax.plot([lo, hi], [yi + offset] * 2, color=colour, lw=1.6, solid_capstyle="round")
        ax.plot(v, yi + offset, "o", color=colour, ms=5.5, zorder=3,
                label=label if r is perf[0] else None)
for yi, r in zip(y, perf):
    tp, fn, fp, tn = (int(r[k]) for k in ("TP", "FN", "FP", "TN"))
    ax.text(105, yi + 0.17, f"{tp}/{tp + fn}", va="center", fontsize=8, color=BLUE)
    ax.text(105, yi - 0.17, f"{tn}/{tn + fp}", va="center", fontsize=8, color=ORANGE)
ax.set_yticks(y, [ABBR[d] for d in drugs])
ax.set_xlim(-2, 100)
ax.set_ylim(-0.7, len(drugs) - 0.3)
ax.set_xticks(range(0, 101, 20), [f"{v}%" for v in range(0, 101, 20)])
ax.grid(axis="x", color=GRID, lw=0.8)
ax.set_axisbelow(True)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.axhline(5.5, color=INK2, lw=0.6, ls=(0, (2, 2)))
ax.text(105, len(drugs) - 0.3, "n", fontsize=8, color=INK2, va="bottom")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, frameon=False, fontsize=8.5,
          handletextpad=0.3, borderaxespad=0)
ax.set_title("A  Agreement with phenotype (95% CI), 48 samples that pass coverage QC", loc="left",
             fontsize=10, weight="bold", pad=22)

# B: QC scatter
ax = fig.add_subplot(grid[0, 1], facecolor=SURFACE)
ax.axvspan(0, 90, color=NEUTRAL, zorder=0)
ax.axhspan(0, 20, color=NEUTRAL, zorder=0)
ax.axvline(90, color=INK2, lw=0.8, ls=(0, (3, 2)))
ax.axhline(20, color=INK2, lw=0.8, ls=(0, (3, 2)))
for s, r in summary.items():
    frac, depth = float(r["mtbc_fraction"]), float(r["mean_depth"])
    fail = r["coverage_qc"] != "PASS"
    ax.scatter(frac, depth, s=26, zorder=3, color=CRITICAL if fail else BLUE,
               edgecolor=SURFACE, lw=0.6)
    stratum = expected[s]["stratum"]
    if fail or stratum in ("overlap_low_mtbc", "suspected_mixed"):
        left = s == "ERR2510541"
        ax.annotate(s, (frac, depth), xytext=(-4 if left else 4, 3), textcoords="offset points", fontsize=7,
                    color=INK2, ha="right" if left else "left")
ax.set_xlim(72, 100)
ax.set_yscale("log")
ax.set_ylim(10, 500)
ax.set_yticks([10, 20, 50, 100, 200, 400], ["10x", "20x", "50x", "100x", "200x", "400x"])
ax.set_xlabel("Reads classified as MTBC by Kraken2 (%)")
ax.set_ylabel("Mean depth on H37Rv")
ax.text(90.3, 420, "--mtb_min_mtbc_fraction 90", fontsize=7.5, color=INK2)
ax.text(72.4, 21, "--mtb_min_mean_depth 20", fontsize=7.5, color=INK2)
low = sum(float(r["mtbc_fraction"]) < 90 for r in summary.values())
ax.text(72.4, 330, f"{low} of 51 below 90% MTBC\n(expected: 1)", fontsize=8, color=INK)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(handles=[Patch(color=BLUE, label="Coverage QC PASS"), Patch(color=CRITICAL, label="Coverage QC FAIL")],
          loc="lower left", frameon=False, fontsize=8)
ax.set_title("B  Sample QC", loc="left", fontsize=10, weight="bold", pad=22)

# C: sample x drug outcome grid
ax = fig.add_subplot(grid[1, :], facecolor=SURFACE)
order = sorted(summary, key=lambda s: (STRATA.index(expected[s]["stratum"]), s))
codes = {k: i for i, (k, _, _) in enumerate(OUTCOMES)}
by_key = {(c["sample"], c["drug"]): c["outcome"] for c in calls}
matrix = [[codes[by_key[(s, d)]] for s in order] for d in drugs]
ax.imshow(matrix, cmap=ListedColormap([c for _, _, c in OUTCOMES]), vmin=-0.5, vmax=len(OUTCOMES) - 0.5,
          aspect="auto", interpolation="nearest")
ax.set_xticks([x - 0.5 for x in range(1, len(order))], minor=True)
ax.set_yticks([y - 0.5 for y in range(1, len(drugs))], minor=True)
ax.grid(which="minor", color=SURFACE, lw=1.2)
ax.tick_params(which="minor", length=0)
ax.set_yticks(range(len(drugs)), [ABBR[d] for d in drugs])
ax.set_xticks(range(len(order)), order, rotation=90, fontsize=6.5)
ax.tick_params(axis="both", length=0)
for side in ax.spines.values():
    side.set_visible(False)
start = 0
for i, s in enumerate(order + [None]):
    if s is None or (i and expected[s]["stratum"] != expected[order[i - 1]]["stratum"]):
        stratum = expected[order[start]]["stratum"]
        ax.plot([start - 0.4, i - 0.6], [-0.9, -0.9], color=INK2, lw=1, clip_on=False)
        small = i - start <= 2
        ax.text((start + i - 1) / 2 - (0.3 if small else 0), -1.2, SHORT.get(stratum, stratum.replace("_", " ")),
                ha="left" if small else "center", va="bottom", fontsize=7, color=INK2, rotation=35 if small else 0,
                rotation_mode="anchor")
        start = i
ax.legend(handles=[Patch(color=c, label=label) for _, label, c in OUTCOMES], loc="upper center",
          bbox_to_anchor=(0.5, -0.3), ncol=3, frameon=False, fontsize=8)
ax.set_title("C  Every call, by sample (grouped by selection stratum) and drug", loc="left",
             fontsize=10, weight="bold", pad=40)

fig.savefig(OUT, dpi=150, facecolor=SURFACE)
print(OUT)
