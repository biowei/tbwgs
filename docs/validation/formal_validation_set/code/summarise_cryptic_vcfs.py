#!/usr/bin/env python3
"""Summarise CRyPTIC per-sample VCFs for sample selection.

For each VCF, write one row with:
  lineage       deepest tbdb barcode lineage with >= 75% of markers matching
  main_lineage  top-level lineage (lineage1 .. lineage7, or La1 etc.)
  mean_depth    median of DP / DPF over called sites (estimated genome depth)
  n_het         sites that fail MIN_FRS (< 90% read support) outside the mask:
                a high count suggests a mixed or contaminated sample
  n_pass_snp    PASS SNPs outside the mask
Also write the PASS SNP and missing positions per sample (JSON) for SNP distances.

Usage: summarise_cryptic_vcfs.py <vcf_dir> <barcode.bed> <genome.fasta> <out_prefix>
"""
import gzip, json, os, statistics, sys
from collections import defaultdict

vcf_dir, barcode_bed, fasta, out = sys.argv[1:5]

seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">"))

markers = defaultdict(list)  # lineage -> [(pos, allele)]
for line in open(barcode_bed):
    f = line.rstrip("\n").split("\t")
    markers[f[3]].append((int(f[2]), f[4]))
marker_pos = {p for m in markers.values() for p, _ in m}


def lineage_call(alleles):
    called = {}
    for lin, ms in markers.items():
        obs = [alleles.get(p, seq[p - 1]) for p, _ in ms]
        good = [(o, a) for o, (_, a) in zip(obs, ms) if o is not None]
        if not good:
            continue
        frac = sum(o == a for o, a in good) / len(good)
        if frac >= 0.75:
            called[lin] = frac
    # deepest lineage whose parents are all called (lineage4.2.2 needs lineage4.2 and lineage4)
    best = None
    for lin in called:
        parts = lin.split(".")
        parents = [".".join(parts[:i]) for i in range(1, len(parts))]
        if all(p in called for p in parents):
            if best is None or len(parts) > len(best.split(".")):
                best = lin
    return best, called


rows, snps = [], {}
for fn in sorted(os.listdir(vcf_dir)):
    if not fn.endswith(".vcf.gz"):
        continue
    alleles, depths, pass_snp, missing = {}, [], [], []
    n_het = 0
    name = None
    with gzip.open(os.path.join(vcf_dir, fn), "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                name = line.rstrip("\n").split("\t")[9]
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos, ref, alts, filt = int(f[1]), f[3], f[4].split(","), f[6]
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fmt.get("GT", "./.")
            masked = "mask-compass" in filt
            try:
                if fmt.get("DPF") not in (None, ".", "0", "0.0") and fmt.get("DP") not in (None, "."):
                    depths.append(float(fmt["DP"]) / float(fmt["DPF"]))
            except (ValueError, ZeroDivisionError):
                pass
            if "MIN_FRS" in filt and not masked:
                n_het += 1
            snp = len(ref) == 1 and all(len(a) == 1 for a in alts)
            if filt == "PASS" and gt not in ("./.", "0/0"):
                allele = alts[int(gt.split("/")[0]) - 1]
                if snp:
                    pass_snp.append(pos)
                    if pos in marker_pos:
                        alleles[pos] = allele
            elif filt == "PASS" and gt == "0/0":
                pass
            elif snp:
                missing.append(pos)
                if pos in marker_pos:
                    alleles[pos] = None
    lin, _ = lineage_call(alleles)
    main = lin.split(".")[0] if lin else None
    rows.append({
        "vcf": fn, "uniqueid": name, "lineage": lin, "main_lineage": main,
        "mean_depth": round(statistics.median(depths), 1) if depths else None,
        "n_het": n_het, "n_pass_snp": len(pass_snp),
    })
    snps[name] = {"pass": pass_snp, "missing": missing}

import csv
with open(out + ".tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
with open(out + ".snps.json", "w") as fh:
    json.dump(snps, fh)
