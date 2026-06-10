#!/usr/bin/env python3
"""
make_samplesheet.py
===================
Build a samples.tsv from directories of inputs. One row per genome with columns:

    id   faa   fna   gff   kofam   gram

The genome id is taken from each protein FASTA filename (minus extension). The
matching nucleotide FASTA, GFF and precomputed KOFAMscan table are found by id in
their respective directories (the dirs may all be the same). Genomes missing any
input are reported and skipped.

Example
-------
    python make_samplesheet.py \
        --faa-dir proteins --fna-dir genomes --gff-dir gffs --kofam-dir kofam \
        --gram both --out config/samples.tsv
"""

import argparse
import glob
import os
import sys

FAA_EXT = [".faa", ".faa.gz", "_protein.faa", "_protein.faa.gz"]
FNA_EXT = [".fna", ".fna.gz", "_genomic.fna", "_genomic.fna.gz", ".fasta"]
GFF_EXT = [".gff", ".gff3", ".gff.gz", ".gff3.gz", "_genomic.gff"]
KO_EXT = [".txt", ".tsv", ".kofam", "_kofam.txt", ".detail.tsv"]


def strip_ext(name, exts):
    for e in sorted(exts, key=len, reverse=True):
        if name.endswith(e):
            return name[: -len(e)]
    return os.path.splitext(name)[0]


def find(directory, gid, exts):
    for e in exts:
        for cand in (os.path.join(directory, gid + e),):
            if os.path.exists(cand):
                return cand
    # Fall back to a glob so prefixed/suffixed names still match.
    hits = sorted(glob.glob(os.path.join(directory, f"*{gid}*")))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description="Build a SOCfinder sample sheet.")
    ap.add_argument("--faa-dir", required=True)
    ap.add_argument("--fna-dir", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--kofam-dir", default=None,
                    help="dir of precomputed KOFAMscan tables; omit if running KOFAMscan "
                         "in-pipeline (run_kofamscan: true)")
    ap.add_argument("--gram", default="both", choices=["p", "n", "both"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    faa_files = []
    for e in FAA_EXT:
        faa_files += glob.glob(os.path.join(args.faa_dir, f"*{e}"))
    faa_files = sorted(set(faa_files))

    rows, missing = [], []
    for faa in faa_files:
        gid = strip_ext(os.path.basename(faa), FAA_EXT)
        fna = find(args.fna_dir, gid, FNA_EXT)
        gff = find(args.gff_dir, gid, GFF_EXT)
        kofam = find(args.kofam_dir, gid, KO_EXT) if args.kofam_dir else None
        need = {"fna": fna, "gff": gff}
        if args.kofam_dir:
            need["kofam"] = kofam
        if not all(need.values()):
            missing.append((gid, need))
            continue
        rows.append((gid, os.path.abspath(faa), os.path.abspath(fna),
                     os.path.abspath(gff),
                     os.path.abspath(kofam) if kofam else "", args.gram))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("id\tfaa\tfna\tgff\tkofam\tgram\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"wrote {len(rows)} genomes to {args.out}")
    if missing:
        print(f"WARNING: skipped {len(missing)} genomes with missing inputs:", file=sys.stderr)
        for gid, info in missing[:20]:
            miss = [k for k, v in info.items() if not v]
            print(f"  {gid}: missing {', '.join(miss)}", file=sys.stderr)
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more", file=sys.stderr)


if __name__ == "__main__":
    main()
