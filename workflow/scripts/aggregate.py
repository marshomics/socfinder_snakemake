#!/usr/bin/env python3
"""
aggregate.py
============
Collect per-genome outputs into three master tables:

  master_SOCKS.tsv     genome_id <tab> protein_id   (one row per social gene)
  master_summary.tsv   per-genome counts (one row per genome)
  missing.txt          genome_ids with no SOCKS.csv (should be empty on success)

Outputs are streamed so memory stays flat regardless of genome count.
"""

import argparse
import csv
import os
import sys


def genome_ids(samples, outdir_root):
    if samples:
        with open(samples) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                yield row["id"]
    else:
        for name in sorted(os.listdir(outdir_root)):
            if os.path.isdir(os.path.join(outdir_root, name)):
                yield name


def main():
    ap = argparse.ArgumentParser(description="Aggregate per-genome SOCfinder outputs.")
    ap.add_argument("--outdir-root", required=True)
    ap.add_argument("--samples", help="sample sheet (defines genome order); optional")
    ap.add_argument("--out-socks", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-missing", required=True)
    args = ap.parse_args()

    n_total = n_done = n_genes = 0
    summary_cols = ["genome_id", "secondary_metabolites", "functional_annotation",
                    "extracellular", "total"]

    with open(args.out_socks, "w", newline="") as socks_fh, \
         open(args.out_summary, "w", newline="") as sum_fh, \
         open(args.out_missing, "w") as miss_fh:
        socks_w = csv.writer(socks_fh, delimiter="\t")
        socks_w.writerow(["genome_id", "protein_id"])
        sum_w = csv.writer(sum_fh, delimiter="\t")
        sum_w.writerow(summary_cols)

        for gid in genome_ids(args.samples, args.outdir_root):
            n_total += 1
            gdir = os.path.join(args.outdir_root, gid)
            socks_path = os.path.join(gdir, "SOCKS.csv")
            if not os.path.exists(socks_path):
                miss_fh.write(gid + "\n")
                continue
            n_done += 1

            with open(socks_path) as fh:
                next(fh, None)  # header
                for line in fh:
                    gene = line.strip().strip('"')
                    if gene:
                        socks_w.writerow([gid, gene])
                        n_genes += 1

            sum_path = os.path.join(gdir, "summary.csv")
            if os.path.exists(sum_path):
                with open(sum_path) as fh:
                    r = csv.DictReader(fh)
                    for row in r:
                        sum_w.writerow([row.get(c, "") for c in summary_cols])

    print(f"aggregated {n_done}/{n_total} genomes, {n_genes} social-gene rows", file=sys.stderr)
    if n_done != n_total:
        print(f"WARNING: {n_total - n_done} genomes missing (see {args.out_missing})", file=sys.stderr)


if __name__ == "__main__":
    main()
