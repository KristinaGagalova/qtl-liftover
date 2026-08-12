#!/usr/bin/env python3
"""Generate a tiny synthetic two-variety dataset with a KNOWN answer:
target chrB1 = source chrA1 with a 1 kb insertion at position 20,000.
Used by run_test.sh together with the tool shims."""
import os
import random
import sys

random.seed(42)
OUT = sys.argv[1] if len(sys.argv) > 1 else "data"
os.makedirs(OUT, exist_ok=True)

QLEN = 49000          # source chrA1
INS_AT, INS_LEN = 20000, 1000

src = "".join(random.choice("ACGT") for _ in range(QLEN))
ins = "".join(random.choice("ACGT") for _ in range(INS_LEN))
tgt = src[:INS_AT] + ins + src[INS_AT:]


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


open(f"{OUT}/varietyA.fa", "w").write(f">chrA1\n{wrap(src)}\n")
open(f"{OUT}/varietyB.fa", "w").write(f">chrB1\n{wrap(tgt)}\n")

# QTL intervals: one before the insertion, one after (so it must shift by 1 kb)
open(f"{OUT}/qtl_intervals.bed", "w").write(
    "chrA1\t5000\t15000\tQTL_before_indel\n"
    "chrA1\t25000\t35000\tQTL_after_indel\n")

# genes every 500 bp, 300 bp CDS
CODON = ["ATG"] + ["GCT", "TTA", "GGA", "AAG", "CGT", "GAA", "TCA"] * 100
with open(f"{OUT}/varietyA.gff3", "w") as g, open(f"{OUT}/varietyA.pep.fa", "w") as p:
    g.write("##gff-version 3\n")
    for i, s in enumerate(range(1000, QLEN - 1000, 500)):
        e = s + 300
        g.write(f"chrA1\t.\tgene\t{s+1}\t{e}\t.\t+\t.\tID=g{i}\n")
        g.write(f"chrA1\t.\tmRNA\t{s+1}\t{e}\t.\t+\t.\tID=g{i}.t1;Parent=g{i}\n")
        g.write(f"chrA1\t.\tCDS\t{s+1}\t{e}\t.\t+\t0\tID=cds{i};Parent=g{i}.t1\n")
        aa = "M" + "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(99))
        p.write(f">g{i}.t1\n{aa}\n")

# probes every 100 bp, 100 bp long
with open(f"{OUT}/probes.fa", "w") as f:
    for i, s in enumerate(range(500, QLEN - 500, 100)):
        f.write(f">P{i}\n{src[s:s+100]}\n")

print(f"test data in {OUT}/ : chrA1={QLEN} bp, chrB1={len(tgt)} bp "
      f"({INS_LEN} bp insertion at {INS_AT})")
