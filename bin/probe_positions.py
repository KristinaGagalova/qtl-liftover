#!/usr/bin/env python3
"""
probe_positions.py -- turn a probe-vs-genome PAF into a table of confident,
uniquely placed probe positions.

Expects something like:
    minimap2 -x sr -N 10 -p 0.5 -t 16 genome.fa probes.fa > probes_vs_genome.paf

A probe is kept only if:
  * best hit covers >= --min-cov of the probe and is >= --min-ident identical
  * the second-best hit has < --max-second x the best hit's matched bases
    (i.e. the placement is unambiguous -- essential, array probes love repeats)

Output: probe  chrom  pos  strand  ident  cov  n_hits
`pos` is the 0-based midpoint of the probe's footprint.
"""

import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cov", type=float, default=0.90)
    ap.add_argument("--min-ident", type=float, default=0.95)
    ap.add_argument("--max-second", type=float, default=0.90)
    ap.add_argument("--report", help="write filtering counts here")
    args = ap.parse_args()

    hits = defaultdict(list)
    with open(args.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            hits[f[0]].append(dict(
                qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                tname=f[5], ts=int(f[7]), te=int(f[8]),
                nmatch=int(f[9]), alnlen=int(f[10]), mapq=int(f[11])))

    n_tot = len(hits)
    n_lowqual = n_multi = n_keep = 0
    with open(args.out, "w") as out:
        out.write("probe\tchrom\tpos\tstrand\tident\tcov\tn_hits\n")
        for probe, hl in hits.items():
            hl.sort(key=lambda h: -h["nmatch"])
            b = hl[0]
            cov = (b["qe"] - b["qs"]) / b["qlen"]
            ident = b["nmatch"] / b["alnlen"] if b["alnlen"] else 0
            if cov < args.min_cov or ident < args.min_ident:
                n_lowqual += 1
                continue
            if len(hl) > 1 and hl[1]["nmatch"] > args.max_second * b["nmatch"]:
                n_multi += 1
                continue
            n_keep += 1
            pos = (b["ts"] + b["te"]) // 2
            out.write(f"{probe}\t{b['tname']}\t{pos}\t{b['strand']}\t"
                      f"{ident:.4f}\t{cov:.3f}\t{len(hl)}\n")

    msg = (f"probes with >=1 hit: {n_tot}\n"
           f"dropped low cov/ident: {n_lowqual}\n"
           f"dropped multi-locus  : {n_multi}\n"
           f"kept unique          : {n_keep}\n")
    print(msg, end="")
    if args.report:
        open(args.report, "w").write(msg)


if __name__ == "__main__":
    main()
