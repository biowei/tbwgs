#!/usr/bin/env python3
"""Figure: tbwgs vs tbtamr per-drug calls, overall agreement, and the QC values behind the thresholds.

Run from the repository root:
    python docs/validation/tbtamr_comparison/code/plot_tbwgs_comparison.py
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tbwgs_vs_tbtamr.png"

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, YELLOW, NEUTRAL, CRITICAL = "#2a78d6", "#eb6834", "#eda100", "#f0efec", "#d03b3b"

SAMPLE_ORDER = ["ERR4813225", "ERR8975345", "ERR8975009", "ERR4831755",
                "ERR3287635", "ERR4811536", "ERR3287775", "ERR4809348"]
DRUGS = ["rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "levofloxacin", "moxifloxacin",
         "bedaquiline", "linezolid", "clofazimine", "delamanid", "amikacin", "kanamycin",
         "capreomycin", "streptomycin", "ethionamide"]
ABBR = {"rifampicin": "RIF", "isoniazid": "INH", "ethambutol": "EMB", "pyrazinamide": "PZA",
        "levofloxacin": "LFX", "moxifloxacin": "MFX", "bedaquiline": "BDQ", "linezolid": "LZD",
        "clofazimine": "CFZ", "delamanid": "DLM", "amikacin": "AMK", "kanamycin": "KAN",
        "capreomycin": "CAP", "streptomycin": "STR", "ethionamide": "ETO"}


def read(path, delimiter="\t"):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def index(rows):
    return {(r["sample"], r["drug"]): r for r in rows}


v100 = index(read(ROOT / "tables/per_drug_calls_v1.0.0.tsv"))
v110 = index(read(ROOT / "tables/per_drug_calls_v1.1.0.tsv"))
qc = {r["sample"]: r for r in read(ROOT / "inputs/samples_qc.tsv")}
profile = {r["sample_id"]: r["profile"] for r in read(ROOT / "inputs/sample_manifest.csv", ",")}


def agreement(calls, key):
    rows = list(calls.values())
    if key == "tbtamr":
        return sum(r["agree_tbtamr"] == "yes" for r in rows), len(rows)
    scored = [r for r in rows if r["phenotype"]]
    if key == "tbwgs_pheno":
        return sum(r["tbwgs"] == r["phenotype"] for r in scored), len(scored)
    return sum(r["tbtamr"] == r["phenotype"] for r in scored), len(scored)


plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                     "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2})
fig = plt.figure(figsize=(13, 8.6), facecolor=SURFACE)
grid = fig.add_gridspec(2, 3, height_ratios=[1.35, 1], width_ratios=[1.25, 1, 1],
                        hspace=0.45, wspace=0.55, left=0.15, right=0.98, top=0.86, bottom=0.08)
fig.suptitle("tbwgs vs tbtamr on 8 CRyPTIC isolates", x=0.02, y=0.985, ha="left", fontsize=14, weight="bold")
fig.text(0.02, 0.945, "tbwgs: TB-Profiler 6.7.0 + tbdb who_v2+ @ c4167762   ·   tbtamr 1.0.4 (WHO catalogue v2)   ·   "
         "phenotype: CRyPTIC binary DST (INH, RIF, EMB only)", color=INK2, fontsize=9)

# A. per-drug call matrix (tbwgs mtb_summary.py 1.1.0 rule), with cells changed from 1.0.0 marked
ax = fig.add_subplot(grid[0, :])
ax.set_facecolor(SURFACE)
colour = {"R": ORANGE, "S": NEUTRAL, "U": YELLOW}
for i, sample in enumerate(SAMPLE_ORDER):
    for j, drug in enumerate(DRUGS):
        new, old = v110[(sample, drug)], v100[(sample, drug)]
        ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, facecolor=colour[new["tbwgs"]], edgecolor="none"))
        label = new["tbwgs"]
        if old["tbwgs"] != new["tbwgs"]:
            label = f"{old['tbwgs']}→{new['tbwgs']}"
            ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fill=False, edgecolor=BLUE, linewidth=2))
        if new["phenotype"] and new["tbwgs"] != new["phenotype"]:
            ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fill=False, edgecolor=CRITICAL,
                                   linewidth=2, linestyle=(0, (3, 2))))
        ax.text(j + 0.5, i + 0.5, label, ha="center", va="center", fontsize=8.5,
                color=INK if label != "S" else INK2, weight="bold" if label != "S" else "normal")
ax.set_xlim(0, len(DRUGS))
ax.set_ylim(len(SAMPLE_ORDER), 0)
ax.set_xticks([j + 0.5 for j in range(len(DRUGS))], [ABBR[d] for d in DRUGS])
ax.xaxis.tick_top()
ax.set_yticks([i + 0.5 for i in range(len(SAMPLE_ORDER))],
              [f"{s}  ({profile[s][0]}/{profile[s][1]}/{profile[s][2]})" for s in SAMPLE_ORDER])
ax.tick_params(length=0)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.axvline(3, color=INK2, linewidth=1)
ax.set_title("A. tbwgs call per sample and drug (rule after PR #9). Every R and S agrees with tbtamr.\n"
             "Phenotype exists only for RIF, INH, EMB (left of the line).",
             loc="left", fontsize=10.5, weight="bold", pad=26)
ax.legend(handles=[
    Patch(facecolor=ORANGE, label="R: resistant in tbwgs and tbtamr"),
    Patch(facecolor=NEUTRAL, label="S: not detected in tbwgs, susceptible in tbtamr"),
    Patch(facecolor=YELLOW, label="U: tbwgs uninterpretable (low coverage), tbtamr susceptible"),
    Patch(facecolor="none", edgecolor=BLUE, linewidth=2, label="U→S: changed by the narrower rule in PR #9"),
    Patch(facecolor="none", edgecolor=CRITICAL, linewidth=2, linestyle="--",
          label="differs from phenotype (rpoB H445L, RIF MIC 0.5)"),
], loc="upper left", bbox_to_anchor=(0, -0.02), ncol=3, frameon=False, fontsize=8.5)

# B. overall agreement
ax = fig.add_subplot(grid[1, 0])
ax.set_facecolor(SURFACE)
bars = [("tbwgs, new rule vs tbtamr", agreement(v110, "tbtamr")),
        ("tbwgs, old rule vs tbtamr", agreement(v100, "tbtamr")),
        ("tbwgs vs phenotype", agreement(v110, "tbwgs_pheno")),
        ("tbtamr vs phenotype", agreement(v110, "tbtamr_pheno"))]
for k, (label, (hit, total)) in enumerate(bars):
    pct = 100 * hit / total
    ax.barh(k, pct, height=0.55, color=BLUE)
    ax.text(pct - 1, k, f"{hit}/{total}  ({pct:.1f}%)", va="center", ha="right", color="white", fontsize=9, weight="bold")
ax.set_yticks(range(len(bars)), [b[0] for b in bars])
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.set_xlabel("agreement (%)")
ax.grid(axis="x", color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.tick_params(length=0)
ax.set_title("B. Categorical agreement", loc="left", fontsize=10.5, weight="bold")


def qc_panel(position, column, threshold, title, xlabel, baseline=None, limits=None):
    ax = fig.add_subplot(position)
    ax.set_facecolor(SURFACE)
    values = [float(qc[s][column]) for s in SAMPLE_ORDER]
    ax.axvline(threshold, color=INK2, linewidth=1.2, linestyle="--")
    ax.text(threshold, -0.9, f"proposed {threshold:g}%", ha="center", fontsize=8, color=INK2)
    if baseline is not None:
        ax.axvline(baseline, color=GRID, linewidth=2)
        ax.text(baseline, -0.9, f"mask {baseline:g}%", ha="center", fontsize=8, color=INK2)
    for i, value in enumerate(values):
        ax.plot([limits[0], value], [i, i], color=GRID, linewidth=1, zorder=1)
        ax.scatter(value, i, s=48, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.text(value, i - 0.32, f"{value:.1f}", ha="center", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER)
    ax.set_ylim(len(SAMPLE_ORDER) - 0.4, -1.3)
    ax.set_xlim(*limits)
    ax.set_xlabel(xlabel)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(title, loc="left", fontsize=10.5, weight="bold")


qc_panel(grid[1, 1], "mtbc_fraction", 90, "C. MTBC fraction (Kraken2)", "reads classified as MTBC (%)", limits=(70, 100))
qc_panel(grid[1, 2], "consensus_n", 10, "D. N in consensus sequence", "consensus positions that are N (%)",
         baseline=7.07, limits=(6, 11))

fig.savefig(OUT, dpi=160, facecolor=SURFACE)
print(f"wrote {OUT}")
