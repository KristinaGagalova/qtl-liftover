#!/usr/bin/env python3
"""Self-test: build a synthetic alignment where the true base-to-base mapping
is known, then check paf_lift reproduces it on both strands."""
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_lift import Aln, lift_point, PafIndex  # noqa: E402

random.seed(7)
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


T = "".join(random.choice("ACGT") for _ in range(2000))

# Build query from T[300:900] with: 50M 5D(target-only) 70M 3I(query-only) 472M
tstart, tend = 300, 900
seg = T[tstart:tend]
ins = "GGG"
q = seg[0:50] + seg[55:125] + ins + seg[125:]
cigar = "50M5D70M3I475M"
qlen = len(q)
# target consumed: 50+5+70+475 = 600 -> matches tend-tstart
assert 50 + 5 + 70 + 475 == tend - tstart
assert 50 + 70 + 3 + 475 == qlen, (qlen,)

# expected mapping query_pos -> target_pos (None inside the insertion)
exp = {}
qi = 0
ti = tstart
for op, ln in [("M", 50), ("D", 5), ("M", 70), ("I", 3), ("M", 475)]:
    for k in range(ln):
        if op == "M":
            exp[qi] = ti
            qi += 1
            ti += 1
        elif op == "D":
            ti += 1
        else:
            exp[qi] = None
            qi += 1

fails = 0

# ---- plus strand -------------------------------------------------------
paf_plus = ["qchr", str(qlen), "0", str(qlen), "+", "tchr", str(len(T)),
            str(tstart), str(tend), str(qlen - 3), str(qlen + 2), "60",
            "tp:A:P", "cg:Z:" + cigar]
a = Aln(paf_plus)
for p in range(qlen):
    r = lift_point(a, p)
    if exp[p] is None:
        continue
    if r is None or r[0] != exp[p]:
        fails += 1
    elif T[r[0]] != q[p]:
        fails += 1
print(f"plus strand : {qlen} positions, {fails} failures")

# ---- minus strand ------------------------------------------------------
q_rc = rc(q)
paf_minus = list(paf_plus)
paf_minus[4] = "-"
b = Aln(paf_minus)
f2 = 0
for p in range(qlen):
    # position p on q_rc corresponds to position qlen-1-p on q
    orig = qlen - 1 - p
    r = lift_point(b, p)
    if exp[orig] is None:
        continue
    if r is None or r[0] != exp[orig]:
        f2 += 1
    elif T[r[0]] != rc(q_rc[p]):
        f2 += 1
print(f"minus strand: {qlen} positions, {f2} failures")

# ---- interval lifting end-to-end --------------------------------------
with tempfile.TemporaryDirectory() as d:
    paf = os.path.join(d, "a.paf")
    with open(paf, "w") as fh:
        fh.write("\t".join(paf_plus) + "\n")
    bed = os.path.join(d, "q.bed")
    with open(bed, "w") as fh:
        fh.write("qchr\t100\t500\tQTL1\n")
    out = os.path.join(d, "o.tsv")
    subprocess.run([sys.executable,
                    os.path.join(os.path.dirname(os.path.abspath(__file__)), "paf_lift.py"),
                    "bed", "--paf", paf, "--bed", bed, "--out", out,
                    "--min-alnlen", "100", "--min-step", "10"], check=True)
    print(open(out).read())

sys.exit(1 if (fails or f2) else 0)
