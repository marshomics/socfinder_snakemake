#!/usr/bin/env python3
"""
test_parsers.py
===============
Unit tests for the ported KOFAM / BLAST / antiSMASH parsers in socfinder_lib.
Builds small synthetic inputs whose correct answers are known by construction and
checks each parser against them. Exits non-zero on any failure.
"""

import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
import socfinder_lib as lib

FAILS = []


def check(name, got, expected):
    ok = got == expected
    print(f"{'OK  ' if ok else 'FAIL'}  {name}: got={got} expected={expected}")
    if not ok:
        FAILS.append(name)


def test_kofam(tmp):
    ko = os.path.join(tmp, "social_ko.csv")
    pd.DataFrame({"term": ["K00001", "K12345"]}).to_csv(ko, index=False)
    table = os.path.join(tmp, "kofam.txt")
    with open(table, "w") as fh:
        fh.write("#\tgene name\tKO\tthrshld\tscore\tE-value\tKO definition\n")
        fh.write("*\tgene1\tK00001\t100\t200\t1e-50\tsocial significant\n")
        fh.write("\tgene2\tK00001\t100\t50\t1e-2\tsocial not significant\n")
        fh.write("*\tgene3\tK99999\t100\t200\t1e-50\tnot social\n")
    out = os.path.join(tmp, "K_SOCK.csv")
    lib.parse_kofam(table, ko, out)
    got = set(pd.read_csv(out)["0"].astype(str))
    check("kofam significant+social", got, {"gene1"})


def test_blast(tmp):
    bdir = os.path.join(tmp, "blast_outputs")
    os.makedirs(bdir, exist_ok=True)
    # columns: sseqid qacc qlen evalue bitscore sstart send slen
    with open(os.path.join(bdir, "x_PSORT_E.txt"), "w") as fh:
        fh.write("s1\tQ1\t300\t0\t100\t1\t100\t300\n")        # rule1: exact -> add Q1
        fh.write("s2\tQ2\t300\t1e-30\t100\t1\t100\t330\n")    # rule4: strong+len -> add Q2
    with open(os.path.join(bdir, "x_PSORT_NE.txt"), "w") as fh:
        fh.write("s3\tQ3\t300\t0\t100\t1\t100\t300\n")        # rule2: exact -> remove Q3
        fh.write("s4\tQ4\t300\t1e-5\t100\t1\t100\t300\n")     # rule3: len-matched -> remove Q4
    with open(os.path.join(bdir, "x_PSORT.txt"), "w") as fh:
        fh.write("s5\tQ5\t300\t0\t100\t1\t100\t300\n")        # rule5: exact -> add Q5
        fh.write("s6\tQ1\t300\t1e-2\t100\t1\t100\t300\n")     # weak, ignored
    out = os.path.join(tmp, "B_SOCK.csv")
    lib.parse_blast(bdir, out)
    got = set(pd.read_csv(out)["combined"].astype(str))
    check("blast extracellular rules", got, {"Q1", "Q2", "Q5"})


def test_antismash(tmp):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO

    adir = os.path.join(tmp, "anti_smash")
    os.makedirs(adir, exist_ok=True)

    def region(path, product, cds):
        rec = SeqRecord(Seq("ATGC" * 250), id="c1", name="c1", description="test")
        rec.annotations["molecule_type"] = "DNA"
        feats = [SeqFeature(FeatureLocation(0, 1000), type="protocluster",
                            qualifiers={"product": [product]})]
        for i, q in enumerate(cds):
            feats.append(SeqFeature(FeatureLocation(i * 100, i * 100 + 90), type="CDS", qualifiers=q))
        rec.features = feats
        SeqIO.write(rec, path, "genbank")

    # Social cluster (terpene): two genes with gene_kind; one lacks protein_id.
    region(os.path.join(adir, "c1.region001.gbk"), "terpene", [
        {"gene_kind": ["biosynthetic"], "product": ["a"], "locus_tag": ["L1"], "protein_id": ["P1"]},
        {"gene_kind": ["biosynthetic-additional"], "product": ["b"], "locus_tag": ["L2"]},  # -> L2
        {"product": ["c"], "locus_tag": ["L3"], "protein_id": ["P3"]},  # no gene_kind -> excluded
    ])
    # Non-social cluster (arylpolyene): its gene must NOT appear.
    region(os.path.join(adir, "c1.region002.gbk"), "arylpolyene", [
        {"gene_kind": ["biosynthetic"], "product": ["d"], "locus_tag": ["L4"], "protein_id": ["P4"]},
    ])

    types = os.path.join(tmp, "antismash_types.csv")
    pd.DataFrame({"Label": ["terpene", "arylpolyene"],
                  "Description": ["t", "a"],
                  "Social": [1, None]}).to_csv(types, index=False)

    out = os.path.join(tmp, "A_SOCK_filtered.csv")
    lib.parse_antismash(adir, types, out)
    df = pd.read_csv(out)
    got = set(df["protein_id"].astype(str))
    check("antismash social genes + protein_id/locus_tag fallback", got, {"P1", "L2"})


def main():
    with tempfile.TemporaryDirectory() as tmp:
        test_kofam(tmp)
        test_blast(tmp)
        try:
            test_antismash(tmp)
        except ImportError:
            print("SKIP antismash test (Biopython not installed)")
    if FAILS:
        print("FAILURES:", FAILS)
        return 1
    print("ALL PARSER TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
