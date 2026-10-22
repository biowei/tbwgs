#!/usr/bin/env python3
"""One row per sample: coverage QC, MTBC fraction, lineage, drug resistance and SNP clusters.

Every sample in --samples gets a row, whether or not TB-Profiler or clustering
produced a result for it. For each drug the call is:

  resistant         a mutation graded "Assoc w R" or "Assoc w R - Interim"
  not_detected      no such mutation, and every gene for the drug is covered
  uninterpretable   no such mutation, but a gene for the drug has a median
                    depth below --min-site-depth
  no_result         TB-Profiler gave no result for the sample

Mutations graded "Uncertain significance", or not in the catalogue, never make
a drug resistant; they are listed in their own column.
"""

import argparse
import csv
import json
import sys

__version__ = "1.0.0"

# First-line drugs first, then the rest of the WHO groups; drugs in the database but not here follow alphabetically
DRUG_ORDER = [
    "rifampicin", "isoniazid", "ethambutol", "pyrazinamide",
    "levofloxacin", "moxifloxacin", "bedaquiline", "linezolid", "clofazimine",
    "delamanid", "pretomanid", "amikacin", "kanamycin", "capreomycin",
    "streptomycin", "ethionamide", "prothionamide", "para-aminosalicylic_acid", "cycloserine", "rifapentine",
]
RESISTANT = ("Assoc w R", "Assoc w R - Interim")


def read_genes(path):
    """Drug -> gene names, from the tbdb genes.bed (chrom, start, end, locus tag, gene, drugs)."""
    drug_genes = {}
    with open(path) as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6 or cols[0].startswith("#"):
                continue
            for drug in cols[5].split(","):
                if drug:
                    drug_genes.setdefault(drug, []).extend(g for g in (cols[3], cols[4]) if g)
    order = {drug: i for i, drug in enumerate(DRUG_ORDER)}
    drugs = sorted(drug_genes, key=lambda d: (order.get(d, len(order)), d))
    return drugs, drug_genes


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mutation_label(variant, drug):
    change = variant.get("change") or variant.get("nucleotide_change") or ""
    return f"{drug}: {variant.get('gene_name') or variant.get('locus_tag', '')} {change} ({variant.get('freq', 0):.2f})"


def profile(result, drugs, drug_genes, min_site_depth):
    depth = {t["target"]: t["median_depth"] for t in result.get("qc", {}).get("target_qc", [])}
    resistant = {}
    graded = {"Assoc w R": [], "Assoc w R - Interim": [], "other": []}
    for variant in result.get("dr_variants", []):
        for entry in variant.get("drugs", []):
            drug = entry.get("drug", "")
            confidence = entry.get("confidence") or "Uncertain significance"
            key = confidence if confidence in RESISTANT else "other"
            graded[key].append(mutation_label(variant, drug))
            if confidence in RESISTANT:
                resistant[drug] = True

    calls = {}
    low_genes = set()
    for drug in drugs:
        low = [g for g in drug_genes[drug] if g in depth and depth[g] < min_site_depth]
        low_genes.update(low)
        if resistant.get(drug):
            calls[drug] = "resistant"
        elif low:
            calls[drug] = "uninterpretable"
        else:
            calls[drug] = "not_detected"

    main_lineage = result.get("main_lineage") or ""
    db = result.get("pipeline", {}).get("db_version", {})
    return {
        "lineage": result.get("sub_lineage") or "",
        "main_lineage": main_lineage,
        "mixed_lineage": "yes" if ";" in main_lineage else "no",
        "drtype": result.get("drtype", ""),
        **calls,
        "resistance_mutations": "; ".join(graded["Assoc w R"]),
        "interim_resistance_mutations": "; ".join(graded["Assoc w R - Interim"]),
        "uncertain_mutations": "; ".join(graded["other"]),
        "low_coverage_genes": ",".join(sorted(low_genes)),
        "tbprofiler_version": result.get("pipeline", {}).get("software_version", ""),
        "tbdb_version": f"{db.get('name', '')}@{db.get('commit', '')}" if db else "",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--samples", required=True, help="TSV: sample, mtbc_fraction, mean_depth, breadth, coverage_qc")
    parser.add_argument("--genes", required=True, help="tbdb genes.bed")
    parser.add_argument("--tbprofiler", nargs="*", default=[], help="TB-Profiler results.json files")
    parser.add_argument("--clusters", help="clusters TSV from snp_cluster.py")
    parser.add_argument("--thresholds", default="5,12", help="comma-separated SNP thresholds")
    parser.add_argument("--min-site-depth", type=float, default=10)
    parser.add_argument("--min-mtbc-fraction", type=float, help="percent; no flag is set when left out")
    parser.add_argument("--prefix", default="mtb_summary")
    args = parser.parse_args()

    drugs, drug_genes = read_genes(args.genes)
    results = {}
    for path in args.tbprofiler:
        with open(path) as handle:
            result = json.load(handle)
        results[result["id"]] = result
    clusters = {row["sample"]: row for row in read_tsv(args.clusters)} if args.clusters else {}
    cluster_columns = [f"cluster_{t.strip()}snp" for t in args.thresholds.split(",") if t.strip()]

    columns = (
        ["sample", "mtbc_fraction", "mtbc_flag", "mean_depth", "breadth", "coverage_qc",
         "lineage", "main_lineage", "mixed_lineage", "drtype"]
        + drugs
        + ["resistance_mutations", "interim_resistance_mutations", "uncertain_mutations", "low_coverage_genes"]
        + cluster_columns
        + ["tbprofiler_version", "tbdb_version", "notes"]
    )

    rows = []
    for sample in sorted(read_tsv(args.samples), key=lambda r: r["sample"]):
        sid = sample["sample"]
        notes = []
        row = {column: "" for column in columns}
        row.update({k: sample.get(k, "") for k in ("sample", "mtbc_fraction", "mean_depth", "breadth", "coverage_qc")})
        if not row["coverage_qc"]:
            row["coverage_qc"] = "no_result"
            notes.append("no coverage result")
        if args.min_mtbc_fraction is None or row["mtbc_fraction"] == "":
            row["mtbc_flag"] = "not_checked"
        else:
            row["mtbc_flag"] = "PASS" if float(row["mtbc_fraction"]) >= args.min_mtbc_fraction else "LOW_MTBC"

        if sid in results:
            row.update(profile(results[sid], drugs, drug_genes, args.min_site_depth))
        else:
            row.update({drug: "no_result" for drug in drugs})
            notes.append("no TB-Profiler result")

        if sid in clusters:
            row.update({c: clusters[sid].get(c, "") for c in cluster_columns})
        else:
            row.update({c: "not_clustered" for c in cluster_columns})
            if row["coverage_qc"] == "PASS":
                notes.append("not clustered: fewer than two sequences or clustering skipped")
            else:
                notes.append("not clustered: coverage QC not PASS")
        row["notes"] = "; ".join(notes)
        rows.append(row)

    with open(f"{args.prefix}.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # MultiQC custom content: same table with a header block
    with open(f"{args.prefix}_mqc.tsv", "w", newline="") as handle:
        handle.write("# id: 'mtb_summary'\n")
        handle.write("# section_name: 'MTB summary'\n")
        handle.write("# description: 'One row per sample: coverage QC, lineage, drug resistance calls and SNP clusters.'\n")
        handle.write("# plot_type: 'table'\n")
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
