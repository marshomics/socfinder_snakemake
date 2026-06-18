#!/usr/bin/env python3
"""
run_one_genome.py
=================
Process a single genome end to end and write its social-gene list.

This replaces SOC_mine.py + SOC_parse.py for one genome, with these structural
changes for scale:

  * KOFAMscan is not run here by default; a precomputed table is passed with
    --kofam (it is the slowest step in the original pipeline). It can still be run
    in-pipeline with --run-kofamscan.
  * BLAST and antiSMASH are run single-threaded by default (--threads 1) so that
    many genomes run concurrently, one core each, instead of one genome hogging 32
    cores. Across 300k genomes this is the single biggest throughput win.
  * Heavy scratch I/O (the cleaned GFF, BLAST outputs, the antiSMASH working
    directory, KOFAMscan's temp dir) is done on node-local disk ($TMPDIR), and
    only the small result CSVs are written back to the shared output directory.
    With thousands of genomes running at once this keeps the shared filesystem
    from becoming the bottleneck. Use --no-local-scratch to work in the output
    directory instead.
  * The step is idempotent: if SOCKS.csv already exists it returns immediately, so
    a re-run of a partially finished batch only does the genomes that are missing.

Exit status is non-zero on failure and SOCKS.csv is not written, so a workflow
manager can safely retry the genome.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import socfinder_lib as lib

GRAM_DB = {"p": "blastdbP", "n": "blastdbN", "both": "blastdbBoth"}
BLAST_OUTFMT = "6 sseqid qacc qlen evalue bitscore sstart send slen"


def run(cmd, log):
    """Run a command, streaming stdout/stderr to a log file. Raise on failure."""
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)


def blastp(db, query, out, threads, log):
    run(["blastp", "-db", db, "-query", query, "-evalue", "10e-8",
         "-outfmt", BLAST_OUTFMT, "-out", out, "-num_threads", str(threads)], log)


def main():
    ap = argparse.ArgumentParser(description="Run SOCfinder on one genome.")
    ap.add_argument("--id", required=True)
    ap.add_argument("--faa", required=True, help="protein FASTA (.faa)")
    ap.add_argument("--fna", required=True, help="nucleotide FASTA (.fna)")
    ap.add_argument("--gff", required=True, help="GFF3")
    ap.add_argument("--kofam", default=None,
                    help="precomputed KOFAMscan detail-tsv table (omit if --run-kofamscan)")
    ap.add_argument("--run-kofamscan", action="store_true",
                    help="run KOFAMscan here instead of reading a precomputed table")
    ap.add_argument("--kofam-bin", default="exec_annotation")
    ap.add_argument("--kofam-profiles", help="KOFAM profiles .hal (e.g. prokaryote.hal)")
    ap.add_argument("--kofam-ko-list", help="KOFAM ko_list file")
    ap.add_argument("--kofam-threshold-scale", default="0.75",
                    help="KOFAMscan --threshold-scale (SOCfinder uses 0.75; tool default is 1.0)")
    ap.add_argument("--kofam-cpus", type=int, default=1)
    ap.add_argument("--kofam-recompute-scale", default="",
                    help="recompute significance from a precomputed scale-1.0 table at this "
                         "scale (e.g. 0.75) instead of trusting its '*' column")
    ap.add_argument("--gram", default="both", choices=["p", "n", "both"])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--blast-db-dir", required=True)
    ap.add_argument("--social-ko", required=True)
    ap.add_argument("--antismash-types", required=True)
    ap.add_argument("--antismash-bin", default="antismash")
    ap.add_argument("--antismash-mode", default="default", choices=["default", "minimal"])
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--scratch-base", default=None,
                    help="base dir for node-local scratch (default: $TMPDIR, else /tmp)")
    ap.add_argument("--no-local-scratch", action="store_true",
                    help="do all work in the output dir instead of node-local scratch")
    ap.add_argument("--force", action="store_true", help="recompute even if SOCKS.csv exists")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="keep cleaned GFF, BLAST outputs and the antiSMASH dir")
    args = ap.parse_args()

    out = args.outdir
    socks = os.path.join(out, "SOCKS.csv")
    failed = os.path.join(out, "FAILED")

    if os.path.exists(socks) and not args.force:
        print(f"SKIP {args.id} (SOCKS.csv exists)")
        return 0

    os.makedirs(out, exist_ok=True)
    if os.path.exists(failed):
        os.remove(failed)
    log = os.path.join(out, "run.log")
    open(log, "w").close()

    # Heavy intermediates live in `work`; small results are written to `out`.
    # By default `work` is a fresh directory on node-local scratch.
    if args.no_local_scratch:
        work, cleanup_work = out, False
    else:
        base = args.scratch_base or os.environ.get("TMPDIR") or tempfile.gettempdir()
        os.makedirs(base, exist_ok=True)
        work = tempfile.mkdtemp(prefix=f"socf_{args.id}_", dir=base)
        cleanup_work = True

    try:
        # --- GFF fix (vectorised, no SQLite) -------------------------------
        clean_gff = os.path.join(work, "clean.gff")
        n_removed = lib.fix_gff(args.gff, clean_gff)
        print(f"{args.id}: {n_removed} GFF records removed")

        # --- BLAST x3 (single-threaded by default) -------------------------
        blast_dir = os.path.join(work, "blast_outputs")
        os.makedirs(blast_dir, exist_ok=True)
        dbd = args.blast_db_dir
        blastp(os.path.join(dbd, GRAM_DB[args.gram]), args.faa,
               os.path.join(blast_dir, "file_PSORT.txt"), args.threads, log)
        blastp(os.path.join(dbd, "blastdbCExtra"), args.faa,
               os.path.join(blast_dir, "file_PSORT_E.txt"), args.threads, log)
        blastp(os.path.join(dbd, "blastdbCNonExtra"), args.faa,
               os.path.join(blast_dir, "file_PSORT_NE.txt"), args.threads, log)

        # --- antiSMASH ------------------------------------------------------
        adir = os.path.join(work, "anti_smash")
        if os.path.exists(adir):
            shutil.rmtree(adir)
        as_cmd = [args.antismash_bin, args.fna,
                  "--genefinding-gff3", clean_gff,
                  "--output-dir", adir,
                  "--cpus", str(args.threads)]
        if args.antismash_mode == "minimal":
            as_cmd.append("--minimal")
        run(as_cmd, log)
        # Record the accession (first token of the nucleotide FASTA header).
        with open(args.fna) as fh:
            first = fh.readline().strip()
        acc = first.split()[0].lstrip(">") if first else args.id
        with open(os.path.join(out, "accession.txt"), "w") as fh:
            fh.write(acc + "\n")

        # --- KOFAMscan (optional; otherwise a precomputed table is used) ----
        # The table itself (-o) is small and goes to the shared outdir so it is
        # kept; only its temp dir (--tmp-dir) lives on node-local scratch.
        if args.run_kofamscan:
            if not (args.kofam_profiles and args.kofam_ko_list):
                raise ValueError("--run-kofamscan needs --kofam-profiles and --kofam-ko-list")
            kofam_path = os.path.join(out, "kofam.txt")
            # Skip the expensive recompute if a good table is already here (retry-safe).
            if args.force or not (os.path.exists(kofam_path) and os.path.getsize(kofam_path) > 0):
                ktmp = os.path.join(work, "kofam_tmp")
                if os.path.exists(ktmp):
                    shutil.rmtree(ktmp)
                run([args.kofam_bin, args.faa, "-o", kofam_path, "-f", "detail-tsv",
                     "--threshold-scale", str(args.kofam_threshold_scale),
                     "--tmp-dir", ktmp, "--cpu", str(args.kofam_cpus),
                     "-p", args.kofam_profiles, "-k", args.kofam_ko_list], log)
                shutil.rmtree(ktmp, ignore_errors=True)
        else:
            kofam_path = args.kofam
            if not kofam_path:
                raise ValueError("no --kofam table given and --run-kofamscan not set")

        # Recompute significance only for a precomputed (scale-1.0) table; a table
        # we just generated is already at the requested scale.
        if args.run_kofamscan or not args.kofam_recompute_scale:
            kofam_recompute = None
        else:
            kofam_recompute = float(args.kofam_recompute_scale)

        # --- Parse each module + combine (small CSVs -> shared outdir) ------
        lib.parse_kofam(kofam_path, args.social_ko, os.path.join(out, "K_SOCK.csv"),
                        threshold_scale=kofam_recompute)
        lib.parse_blast(blast_dir, os.path.join(out, "B_SOCK.csv"))
        lib.parse_antismash(adir, args.antismash_types, os.path.join(out, "A_SOCK_filtered.csv"))
        summary = lib.combine(out, genome_id=args.id)
        print(f"{args.id}: total social genes = {summary['total']}")

        # --- Intermediates --------------------------------------------------
        if args.keep_intermediates and work != out:
            # Copy the heavy outputs from node-local scratch back to the outdir.
            shutil.copytree(adir, os.path.join(out, "anti_smash"), dirs_exist_ok=True)
            shutil.copytree(blast_dir, os.path.join(out, "blast_outputs"), dirs_exist_ok=True)
            shutil.copyfile(clean_gff, os.path.join(out, "clean.gff"))
        elif not args.keep_intermediates and work == out:
            shutil.rmtree(adir, ignore_errors=True)
            shutil.rmtree(blast_dir, ignore_errors=True)
            if os.path.exists(clean_gff):
                os.remove(clean_gff)

    except Exception as exc:  # noqa: BLE001 - we want any failure to be recoverable
        with open(failed, "w") as fh:
            fh.write(f"{type(exc).__name__}: {exc}\n")
        # Remove a half-written SOCKS so the genome is retried, not skipped.
        if os.path.exists(socks):
            os.remove(socks)
        print(f"FAILED {args.id}: {exc}", file=sys.stderr)
        return 1
    finally:
        # Always clear node-local scratch, even on failure.
        if cleanup_work and os.path.isdir(work):
            shutil.rmtree(work, ignore_errors=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
