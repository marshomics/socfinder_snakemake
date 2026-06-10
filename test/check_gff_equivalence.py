#!/usr/bin/env python3
"""
check_gff_equivalence.py
========================
Prove the vectorised fix_gff() keeps exactly the same feature rows as the
original SOC_mine.py GFF-fixing code.

It reimplements the original row-removal logic verbatim (rules A/B/C) and compares
the surviving features (by seqid, feature type, start, end, strand) against the
output of socfinder_lib.fix_gff. Exits non-zero on any mismatch.

Usage:
    python check_gff_equivalence.py GFF [GFF ...]
"""

import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
import socfinder_lib as lib


def original_kept(gff_path):
    """The original SOC_mine.py logic, returning the set of surviving features."""
    header_lines, data = [], []
    for line in open(gff_path):
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            header_lines.append(line.strip())
        else:
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, frame, attributes = fields
            data.append({
                "Sequence ID": seqid, "Source": source, "Feature Type": ftype,
                "Start": int(start), "End": int(end), "Score": score,
                "Strand": strand, "Frame": frame, "Attributes": attributes,
            })
    df = pd.DataFrame(data)
    if df.empty:
        return set()

    region_ends, rows_to_remove = {}, []
    # Rule A
    for idx, row in df.iterrows():
        if row["Feature Type"] == "region":
            region_ends[row["Sequence ID"]] = row["End"]
        else:
            if row["Sequence ID"] in region_ends and row["End"] > region_ends[row["Sequence ID"]]:
                rows_to_remove.append(idx)
    # Rule B (feature-type string longer than 15 chars)
    for idx, row in df.iterrows():
        if len(row["Feature Type"]) > 15:
            rows_to_remove.append(idx)
    # Rule C (overlapping CDS sharing End -> keep longest)
    cds = df[df["Feature Type"] == "CDS"].copy()
    cds.sort_values(by=["Sequence ID", "End", "Start"], inplace=True)
    for seqid, group in cds.groupby("Sequence ID"):
        for end, sub in group.groupby("End"):
            if len(sub) > 1:
                longest = sub.loc[(sub["End"] - sub["Start"]).idxmax()].name
                rows_to_remove.extend(sub.index.difference([longest]))

    kept = df.drop(index=set(rows_to_remove))
    return set(zip(kept["Sequence ID"], kept["Feature Type"], kept["Start"],
                   kept["End"], kept["Strand"]))


def new_kept(gff_path):
    with tempfile.NamedTemporaryFile("w", suffix=".gff", delete=False) as tmp:
        out = tmp.name
    lib.fix_gff(gff_path, out)
    rows = []
    for line in open(out):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6]))
    os.remove(out)
    return set(rows)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    all_ok = True
    for gff in sys.argv[1:]:
        o, n = original_kept(gff), new_kept(gff)
        ok = (o == n)
        all_ok &= ok
        print(f"{'OK  ' if ok else 'FAIL'}  {gff}: original kept {len(o)}, new kept {len(n)}, "
              f"identical={ok}")
        if not ok:
            print("  only in original:", list(o - n)[:5])
            print("  only in new     :", list(n - o)[:5])
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
