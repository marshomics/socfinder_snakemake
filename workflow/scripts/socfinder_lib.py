#!/usr/bin/env python3
"""
socfinder_lib.py
================
Shared, dependency-light reimplementation of the SOCfinder parsing logic.

Why this exists
---------------
The original SOCfinder splits work across SOC_mine.py (runs the tools) and
SOC_parse.py (turns tool output into social-gene lists). Two parts of that code
do not scale to hundreds of thousands of genomes:

  * SOC_mine.py builds a gffutils SQLite database for every genome
    (``gffutils.create_db``) and then never uses it, before re-parsing the GFF by
    hand with two ``DataFrame.iterrows()`` loops. Per genome that is a wasted disk
    write plus slow Python loops; at 300k genomes it is a real I/O + CPU sink.
  * SOC_parse.py parses antiSMASH GenBank files with ``exec()`` to build
    ``vec0, vec1, ...`` variables that it reads back out of ``locals()``. That is
    O(n^2), fragile, and can fail on large genomes.

This module rewrites both as vectorised pandas / Biopython. The goal is identical
social-gene *calls* to the original, just faster and robust. Each function mirrors
one block of the original SOC_parse.py and is documented with the original logic
it replaces.
"""

from __future__ import annotations

import os
import glob
import pandas as pd
import numpy as np

# Column order of a GFF3 feature line.
GFF_COLS = ["seqid", "source", "ftype", "start", "end", "score", "strand", "frame", "attributes"]


# ---------------------------------------------------------------------------
# 1. GFF fixing  (replaces the gffutils + iterrows block in SOC_mine.py)
# ---------------------------------------------------------------------------
def fix_gff(in_gff: str, out_gff: str) -> int:
    """Clean a GFF so antiSMASH accepts it, returning the number of rows removed.

    Replicates SOC_mine.py's three removal rules, vectorised:
      A) drop non-'region' features whose End runs past the contig's 'region' end;
      B) drop features whose feature-type string is longer than 15 characters
         (antiSMASH-incompatible junk feature types);
      C) for CDS sharing the same (seqid, End), keep only the longest.

    Differences from the original that do NOT change which rows survive:
      * no gffutils SQLite database is created (the original built one and threw
        it away);
      * the attributes column is kept verbatim instead of being round-tripped
        through a Python dict, which is both faster and safer (the dict round-trip
        crashes on any attribute value containing '=').
    """
    header_lines = []
    rows = []
    with open(in_gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                # Malformed / non-feature line: skip (original would have crashed).
                continue
            rows.append(parts)

    if not rows:
        # Nothing to clean: write headers only so downstream tools still get a file.
        with open(out_gff, "w") as out:
            for h in header_lines:
                out.write(h + "\n")
        return 0

    df = pd.DataFrame(rows, columns=GFF_COLS)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")

    remove = pd.Series(False, index=df.index)

    # Rule A: feature end beyond the contig 'region' end.
    region_rows = df[df["ftype"] == "region"]
    if not region_rows.empty:
        # Original keeps the LAST 'region' end seen per seqid.
        region_end = (
            region_rows.drop_duplicates("seqid", keep="last")
            .set_index("seqid")["end"]
        )
        mapped = df["seqid"].map(region_end)
        remove |= (df["ftype"] != "region") & mapped.notna() & (df["end"] > mapped)

    # Rule B: absurdly long feature-type strings.
    remove |= df["ftype"].str.len() > 15

    # Rule C: overlapping CDS with identical End -> keep the longest.
    cds = df[df["ftype"] == "CDS"].copy()
    if not cds.empty:
        cds["_len"] = cds["end"] - cds["start"]
        cds_sorted = cds.sort_values("_len", ascending=False, kind="mergesort")
        dup = cds_sorted.duplicated(subset=["seqid", "end"], keep="first")
        remove.loc[dup[dup].index] = True

    kept = df[~remove]
    n_removed = int(remove.sum())

    with open(out_gff, "w") as out:
        for h in header_lines:
            out.write(h + "\n")
        kept.to_csv(out, sep="\t", header=False, index=False)

    return n_removed


# ---------------------------------------------------------------------------
# 2. KOFAM module  (replaces the KOFAM block in SOC_parse.py)
# ---------------------------------------------------------------------------
def parse_kofam(kofam_path: str, social_ko_path: str, out_csv: str) -> int:
    """Write K_SOCK.csv (unique social gene names) and return its row count.

    KOFAMscan must have been run with ``-f detail-tsv``. Significant hits are the
    rows whose first column is '*' (KOFAMscan marks these when score >= threshold).
    To match SOCfinder exactly the table should have been produced with
    ``--threshold-scale 0.75`` (see the test/ command in the upstream repo).
    """
    sock = pd.read_csv(social_ko_path)
    cols = ["#", "gene name", "KO", "thrshld", "score", "E-value", "KO definition"]

    if os.path.getsize(kofam_path) == 0:
        pd.DataFrame({"0": []}).to_csv(out_csv, index=False)
        return 0

    data = pd.read_csv(
        kofam_path, sep="\t", quoting=3, header=None, names=cols,
        dtype=str, engine="python",
    )
    data = data[data["#"] == "*"]
    data = data[data["KO"].isin(sock["term"])]
    socks = pd.unique(data["gene name"])
    pd.DataFrame(socks, columns=["0"]).to_csv(out_csv, index=False)
    return len(socks)


# ---------------------------------------------------------------------------
# 3. BLAST module  (replaces the BLASTE block in SOC_parse.py)
# ---------------------------------------------------------------------------
_BLAST_COLS = ["subject_seq_id", "query_acc", "query_length", "evalue",
               "bitscore", "subject_start", "subject_end", "subject_length"]


def _read_blast(path: str) -> pd.DataFrame:
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        df = pd.read_table(path, header=None, sep=r"\s+")
        df.columns = _BLAST_COLS[: df.shape[1]]
        return df
    return pd.DataFrame(columns=_BLAST_COLS)


def parse_blast(blast_dir: str, out_csv: str) -> int:
    """Write B_SOCK.csv (extracellular social genes) and return its row count.

    Faithful port of the SOC_parse.py extracellular logic:
      data    = blast vs the gram PSORTb database (computational extracellular)
      data_e  = blast vs experimentally-verified extracellular proteins
      data_ne = blast vs experimentally-verified non-extracellular proteins

      1) add query if it has an exact (evalue==0) hit to experimental extracellular
      2) drop query if it has an exact hit to experimental non-extracellular
      3) drop query if it has a length-matched (+/-10%) hit to non-extracellular
      4) add query with a strong (evalue<1e-19) length-matched (+/-20%) hit to
         experimental extracellular
      5) add query with an exact hit to the computational PSORTb database
    """
    data = _read_blast(glob.glob(os.path.join(blast_dir, "*_PSORT.txt"))[0]
                       if glob.glob(os.path.join(blast_dir, "*_PSORT.txt")) else "")
    data_e = _read_blast(glob.glob(os.path.join(blast_dir, "*_PSORT_E.txt"))[0]
                        if glob.glob(os.path.join(blast_dir, "*_PSORT_E.txt")) else "")
    data_ne = _read_blast(glob.glob(os.path.join(blast_dir, "*_PSORT_NE.txt"))[0]
                        if glob.glob(os.path.join(blast_dir, "*_PSORT_NE.txt")) else "")

    sock = []

    # 1) exact match in experimental extracellular
    if not data_e.empty:
        sock = list(pd.unique(data_e.loc[data_e["evalue"] == 0, "query_acc"]))

    # 2) exact match in experimental non-extracellular -> remove
    if not data_ne.empty:
        removies = set(pd.unique(data_ne.loc[data_ne["evalue"] == 0, "query_acc"]))
        if removies:
            if not data.empty:
                data = data[~data["query_acc"].isin(removies)]
            if not data_e.empty:
                data_e = data_e[~data_e["query_acc"].isin(removies)]
            sock = [g for g in sock if g not in removies]

    # 3) length-matched hit in experimental non-extracellular -> remove
    if not data_ne.empty:
        cond = (data_ne["subject_length"] >= data_ne["query_length"] * 0.90) & \
               (data_ne["subject_length"] <= data_ne["query_length"] * 1.10)
        removies2 = set(pd.unique(data_ne.loc[cond, "query_acc"]))
        if removies2:
            if not data.empty:
                data = data[~data["query_acc"].isin(removies2)]
            if not data_e.empty:
                data_e = data_e[~data_e["query_acc"].isin(removies2)]
            sock = [g for g in sock if g not in removies2]

    sock_df = pd.DataFrame({"SOCK": sock})

    # 4) strong, length-matched hit in experimental extracellular -> add
    if not data_e.empty:
        cond = (data_e["evalue"] < 1e-19) & \
               (data_e["subject_length"] >= data_e["query_length"] * 0.80) & \
               (data_e["subject_length"] <= data_e["query_length"] * 1.20)
        add = pd.unique(data_e.loc[cond, "query_acc"])
        if len(add):
            sock_df = pd.concat([sock_df, pd.DataFrame({"SOCK": add})], ignore_index=True)

    # 5) exact match in computational PSORTb -> add
    if not data.empty:
        add = pd.unique(data.loc[data["evalue"] == 0, "query_acc"])
        if len(add):
            sock_df = pd.concat([sock_df, pd.DataFrame({"SOCK": add})], ignore_index=True)

    sock_df = sock_df.drop_duplicates(subset="SOCK")
    out = sock_df["SOCK"].astype(str).tolist() if not sock_df.empty else []
    pd.DataFrame(out, columns=["combined"]).to_csv(out_csv, index=False)
    return len(out)


# ---------------------------------------------------------------------------
# 4. antiSMASH module  (replaces the exec()-based GenBank parser in SOC_parse.py)
# ---------------------------------------------------------------------------
def _clean(value: str) -> str:
    return value.replace('"', "").replace(" ", "").replace("/", "")


def parse_antismash(anti_dir: str, antismash_types_path: str, out_csv: str) -> int:
    """Write A_SOCK_filtered.csv (social secondary-metabolite genes); return its rows.

    For every ``*region*.gbk`` antiSMASH file:
      * the cluster type is the product of the file's first 'protocluster' feature
        (matching the original, which keyed on the first protocluster);
      * a CDS is kept if it carries a ``gene_kind`` qualifier (the original kept
        rows whose gene_kind cell was populated) AND a ``product`` qualifier;
      * its identifier is ``protein_id``, falling back to ``locus_tag``.
    A gene is social if its cluster type is in antismash_types.csv where Social==1.

    Uses Biopython's GenBank parser instead of exec()/locals() string surgery.
    """
    from Bio import SeqIO  # imported lazily so the lib loads without Biopython

    types = pd.read_csv(antismash_types_path, encoding="latin-1")
    social_types = set(types.loc[types["Social"] == 1, "Label"].tolist())

    records = []
    files = sorted(f for f in glob.glob(os.path.join(anti_dir, "*.gbk"))
                   if "region" in os.path.basename(f))

    for region_idx, path in enumerate(files):
        for rec in SeqIO.parse(path, "genbank"):
            region_type = None
            for feat in rec.features:
                if feat.type == "protocluster":
                    prod = feat.qualifiers.get("product", [None])[0]
                    if prod is not None:
                        region_type = _clean(prod)
                    break
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                if "gene_kind" not in feat.qualifiers:
                    continue
                has_product = "product" in feat.qualifiers
                locus_tag = feat.qualifiers.get("locus_tag", [None])[0]
                protein_id = feat.qualifiers.get("protein_id", [None])[0]
                records.append({
                    "gene_kind": "gene_kind=" + _clean(feat.qualifiers["gene_kind"][0]),
                    "product": region_type if has_product else None,
                    "locus_tag": _clean(locus_tag) if locus_tag else None,
                    "protein_id": _clean(protein_id) if protein_id else None,
                    "region": region_idx,
                })

    df = pd.DataFrame(records, columns=["gene_kind", "product", "locus_tag", "protein_id", "region"])
    if not df.empty:
        df["protein_id"] = df["protein_id"].astype(object)
        df["protein_id"] = df["protein_id"].fillna(df["locus_tag"])
        df = df[df["product"].isin(social_types)]
    df.to_csv(out_csv, index=False)
    return len(df)


# ---------------------------------------------------------------------------
# 5. Combine  (replaces the COMBINE block in SOC_parse.py)
# ---------------------------------------------------------------------------
def combine(out_dir: str, genome_id: str | None = None) -> dict:
    """Merge the three module outputs into SOCKS.csv + summary.csv for one genome."""
    a_path = os.path.join(out_dir, "A_SOCK_filtered.csv")
    k_path = os.path.join(out_dir, "K_SOCK.csv")
    b_path = os.path.join(out_dir, "B_SOCK.csv")

    data_a = pd.read_csv(a_path) if os.path.exists(a_path) else pd.DataFrame(columns=["protein_id"])
    if "protein_id" in data_a.columns and not data_a.empty:
        data_a["protein_id"] = data_a["protein_id"].astype(str).str.replace('"', "", regex=False)
        data_a = pd.DataFrame(data_a["protein_id"])
    else:
        data_a = pd.DataFrame(columns=["protein_id"])

    data_k = pd.read_csv(k_path) if os.path.exists(k_path) else pd.DataFrame(columns=["0"])
    data_b = pd.read_csv(b_path) if os.path.exists(b_path) else pd.DataFrame(columns=["combined"])
    data_k = data_k.rename(columns={"0": "protein_id"})
    data_b = data_b.rename(columns={"combined": "protein_id"})

    combined = pd.concat([data_a, data_b, data_k], ignore_index=True)
    if combined.empty:
        unique_entries = pd.Series([], dtype=str, name="protein_id")
    else:
        unique_entries = combined.iloc[:, 0].drop_duplicates()
    unique_entries.to_csv(os.path.join(out_dir, "SOCKS.csv"), index=False)

    summary = pd.DataFrame({
        "genome_id": [genome_id if genome_id else os.path.basename(out_dir.rstrip("/"))],
        "secondary_metabolites": [len(data_a)],
        "functional_annotation": [len(data_k)],
        "extracellular": [len(data_b)],
        "total": [len(unique_entries)],
    })
    summary.to_csv(os.path.join(out_dir, "summary.csv"), index=False)
    return summary.iloc[0].to_dict()
