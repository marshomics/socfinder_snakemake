#!/usr/bin/env python3
"""
run_batch.py
============
Process one batch of genomes inside a single cluster job.

A batch is just a slice of the sample sheet. Running ~hundreds of genomes per
cluster job (instead of one job per genome) keeps both the scheduler queue and
the Snakemake DAG small enough to handle 300k genomes. Within the job, genomes
run concurrently, one core each, via a thread pool whose size is the number of
cores the job was given.

Genome-level resume: each genome is processed by run_one_genome.py, which skips
itself if SOCKS.csv already exists. The batch sentinel is only written once every
genome in the batch has a SOCKS.csv, so a failed/partial batch can be re-run and
will only redo the genomes that are still missing.
"""

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
RUN_ONE = os.path.join(HERE, "run_one_genome.py")


def read_rows(path):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def process(row, args):
    gid = row["id"]
    outdir = os.path.join(args.outdir_root, gid)
    gram = (row.get("gram") or args.gram_default).strip().lower()
    kofam = (row.get("kofam") or "").strip()
    cmd = [
        sys.executable, RUN_ONE,
        "--id", gid,
        "--faa", row["faa"], "--fna", row["fna"], "--gff", row["gff"],
        "--gram", gram,
        "--outdir", outdir,
        "--blast-db-dir", args.blast_db_dir,
        "--social-ko", args.social_ko,
        "--antismash-types", args.antismash_types,
        "--antismash-bin", args.antismash_bin,
        "--antismash-mode", args.antismash_mode,
        "--threads", str(args.threads_per_genome),
    ]
    if args.run_kofamscan:
        cmd += ["--run-kofamscan",
                "--kofam-bin", args.kofam_bin,
                "--kofam-profiles", args.kofam_profiles,
                "--kofam-ko-list", args.kofam_ko_list,
                "--kofam-threshold-scale", str(args.kofam_threshold_scale),
                "--kofam-cpus", str(args.kofam_cpus)]
    elif kofam:
        cmd += ["--kofam", kofam]
    if args.no_local_scratch:
        cmd.append("--no-local-scratch")
    if args.scratch_base:
        cmd += ["--scratch-base", args.scratch_base]
    if args.keep_intermediates:
        cmd.append("--keep-intermediates")
    if args.force:
        cmd.append("--force")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    ok = os.path.exists(os.path.join(outdir, "SOCKS.csv"))
    return gid, ok, proc.returncode, proc.stdout, proc.stderr


def main():
    ap = argparse.ArgumentParser(description="Run SOCfinder on a batch of genomes.")
    ap.add_argument("--batch-tsv", required=True)
    ap.add_argument("--sentinel", required=True, help="file to touch when all genomes succeed")
    ap.add_argument("--jobs", type=int, default=1, help="genomes to run concurrently")
    ap.add_argument("--threads-per-genome", type=int, default=1)
    ap.add_argument("--outdir-root", required=True)
    ap.add_argument("--blast-db-dir", required=True)
    ap.add_argument("--social-ko", required=True)
    ap.add_argument("--antismash-types", required=True)
    ap.add_argument("--antismash-bin", default="antismash")
    ap.add_argument("--antismash-mode", default="default", choices=["default", "minimal"])
    ap.add_argument("--run-kofamscan", action="store_true")
    ap.add_argument("--kofam-bin", default="exec_annotation")
    ap.add_argument("--kofam-profiles", default="")
    ap.add_argument("--kofam-ko-list", default="")
    ap.add_argument("--kofam-threshold-scale", default="0.75")
    ap.add_argument("--kofam-cpus", type=int, default=1)
    ap.add_argument("--scratch-base", default="")
    ap.add_argument("--no-local-scratch", action="store_true")
    ap.add_argument("--gram-default", default="both", choices=["p", "n", "both"])
    ap.add_argument("--keep-intermediates", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    rows = list(read_rows(args.batch_tsv))
    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(process, r, args) for r in rows]
        for fut in as_completed(futs):
            gid, ok, rc, _out, err = fut.result()
            if not ok:
                failures.append(gid)
                tail = (err or "").strip().splitlines()[-3:]
                print(f"FAILED {gid} (rc={rc}): {' | '.join(tail)}", file=sys.stderr)

    done = sum(os.path.exists(os.path.join(args.outdir_root, r["id"], "SOCKS.csv")) for r in rows)
    print(f"batch: {done}/{len(rows)} genomes complete")

    if done == len(rows):
        os.makedirs(os.path.dirname(args.sentinel), exist_ok=True)
        open(args.sentinel, "w").close()
        return 0
    print(f"{len(rows) - done} genomes incomplete; sentinel not written "
          f"(re-run to retry only the missing ones)", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
