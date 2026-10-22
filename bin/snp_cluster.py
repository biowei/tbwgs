#!/usr/bin/env python3
"""Single-linkage SNP clusters from a snp-dists matrix.

Two samples join the same cluster when a chain of pairwise distances at or
below the threshold links them. One column of cluster IDs is written per
threshold; samples without a partner at that threshold get "-".
"""

import argparse
import csv
import sys


def read_matrix(path, exclude):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    names = rows[0][1:]
    dists = {}
    for row in rows[1:]:
        for name, value in zip(names, row[1:]):
            dists[(row[0], name)] = int(value)
    keep = [name for name in names if name not in exclude]
    return keep, dists


def clusters_at(names, dists, threshold):
    parent = {name: name for name in names}

    def find(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            if dists[(a, b)] <= threshold:
                parent[find(a)] = find(b)

    groups = {}
    for name in names:
        groups.setdefault(find(name), []).append(name)
    # Number clusters of two or more by their first member, so IDs are stable for the same input
    labels = {name: "-" for name in names}
    multi = sorted((sorted(members) for members in groups.values() if len(members) > 1), key=lambda m: m[0])
    for number, members in enumerate(multi, start=1):
        for name in members:
            labels[name] = f"T{threshold}-{number}"
    return labels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", help="snp-dists TSV")
    parser.add_argument("output", help="clusters TSV")
    parser.add_argument("--thresholds", default="5,12", help="comma-separated SNP thresholds")
    parser.add_argument("--exclude", action="append", default=[], help="sequence name to leave out, such as the reference")
    args = parser.parse_args()

    thresholds = [int(value) for value in args.thresholds.split(",") if value.strip()]
    names, dists = read_matrix(args.matrix, set(args.exclude))
    if not names:
        sys.exit("No samples left in the distance matrix")

    columns = {threshold: clusters_at(names, dists, threshold) for threshold in thresholds}
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample"] + [f"cluster_{threshold}snp" for threshold in thresholds])
        for name in sorted(names):
            writer.writerow([name] + [columns[threshold][name] for threshold in thresholds])


if __name__ == "__main__":
    main()
