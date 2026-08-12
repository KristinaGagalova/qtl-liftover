#!/usr/bin/env python3
"""
gff_to_protein_bed.py -- build the protein_id -> source-genome BED that
miniprot_lift.py needs, plus a protein length table.

The IDs in column 4 MUST match the FASTA headers of the protein file you
feed to miniprot. Use --id-attr to pick the right attribute (ID, protein_id,
transcript_id, Name...) and --strip-prefix for GFF3 IDs like "mRNA:XYZ".

Only ONE representative transcript per gene is kept (the longest CDS), which
is what you want for anchors -- redundant isoforms just add noise.
"""

import argparse
import re
import sys
from collections import defaultdict

ATTR = re.compile(r"([^;=\s]+)=([^;]*)")


def attrs(s):
    d = dict(ATTR.findall(s))
    if not d:  # GTF style
        d = dict(re.findall(r'(\S+) "([^"]*)"', s))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--out-len", help="protein_id<TAB>aa_length")
    ap.add_argument("--feature", default="mRNA",
                    help="transcript feature type (mRNA / transcript)")
    ap.add_argument("--id-attr", default="ID")
    ap.add_argument("--parent-attr", default="Parent")
    ap.add_argument("--strip-prefix", action="store_true",
                    help="drop a leading 'type:' from IDs (Ensembl-style GFF3)")
    args = ap.parse_args()

    def clean(x):
        return x.split(":", 1)[1] if (args.strip_prefix and ":" in x) else x

    tx = {}        # tid -> (chrom, start, end, gene)
    cdslen = defaultdict(int)
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == args.feature:
                tid = clean(a.get(args.id_attr, ""))
                if not tid:
                    continue
                tx[tid] = (f[0], int(f[3]) - 1, int(f[4]),
                           clean(a.get(args.parent_attr, tid)))
            elif f[2] == "CDS":
                for par in a.get(args.parent_attr, "").split(","):
                    par = clean(par)
                    if par:
                        cdslen[par] += int(f[4]) - int(f[3]) + 1

    best = {}
    for tid, (c, s, e, gene) in tx.items():
        L = cdslen.get(tid, e - s)
        if gene not in best or L > best[gene][0]:
            best[gene] = (L, tid, c, s, e)

    with open(args.out_bed, "w") as bo:
        for gene, (L, tid, c, s, e) in sorted(best.items(), key=lambda x: (x[1][2], x[1][3])):
            bo.write(f"{c}\t{s}\t{e}\t{tid}\t{gene}\n")
    if args.out_len:
        with open(args.out_len, "w") as lo:
            for gene, (L, tid, *_) in best.items():
                lo.write(f"{tid}\t{max(L // 3 - 1, 1)}\n")
    sys.stderr.write(f"[gff_to_protein_bed] {len(best)} representative "
                     f"transcripts from {len(tx)} total\n")


if __name__ == "__main__":
    main()
