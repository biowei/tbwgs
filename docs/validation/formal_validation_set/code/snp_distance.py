#!/usr/bin/env python3
"""Pairwise SNP distance from summarise_cryptic_vcfs.py output.

Distance = PASS SNP positions that differ between the two samples, with sites
that are filtered (missing) in either sample left out. Masked sites are never PASS.

Pairs are all isolates from the same CRyPTIC site and subject that have a VCF.

Usage: snp_distance.py <snps.json> <cryptic_pool_annotated.tsv> <out.tsv>
"""
import csv, itertools, json, sys
from collections import defaultdict

snps = json.load(open(sys.argv[1]))
patients = defaultdict(list)
with open(sys.argv[2]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["uniqueid"]:
            patients[(r["site"], r["subj"])].append(r["uniqueid"])
with open(sys.argv[3], "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["id_a", "id_b", "snp_distance"])
    for a, b in (p for ids in patients.values() for p in itertools.combinations(ids, 2)):
        if a not in snps or b not in snps:
            continue
        pa, pb = set(snps[a]["pass"]), set(snps[b]["pass"])
        miss = set(snps[a]["missing"]) | set(snps[b]["missing"])
        w.writerow([a, b, len((pa ^ pb) - miss)])
