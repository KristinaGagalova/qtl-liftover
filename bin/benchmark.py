#!/usr/bin/env python3
"""
benchmark.py -- score the two liftover routes against probe-derived truth.

Ground truth = probes that map uniquely and confidently to BOTH assemblies.
Their (source_pos, target_pos) pairs are an independent, method-free answer
key for the coordinate transfer.

Two levels of evaluation
------------------------
A. Base-level (whole-genome-alignment route only, since only that route
   defines a coordinate function):
     - lift-over rate, wrong-chromosome rate
     - |lifted - true| error: median / mean / P90 / P99
     - fraction within 1 kb / 10 kb / 100 kb

B. Interval-level (BOTH routes, this is what actually matters for a QTL):
   For each QTL, the truth set is the probes inside it in the source
   assembly; their true target positions define the region the transferred
   interval SHOULD contain.
     - recall    = truth probes falling inside the transferred interval
     - precision = probes inside the transferred interval that were inside
                   the QTL in the source (i.e. how much junk was dragged in)
     - Jaccard, boundary errors, length inflation

Note on the answer key: probes were designed on one genome, so they are
ascertainment-biased toward regions that are present and unique in that
genome. Recall is therefore optimistic for the exact regions the probes
cover and says nothing about probe deserts. Read the numbers with that in
mind -- they compare the two routes fairly, but they are not an absolute
accuracy for the whole interval.
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def read_tsv(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def quant(v, q):
    if not v:
        return float("nan")
    v = sorted(v)
    i = min(int(q * (len(v) - 1) + 0.5), len(v) - 1)
    return v[i]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--probes-src", required=True, help="probe_positions.py on source")
    ap.add_argument("--probes-tgt", required=True, help="probe_positions.py on target")
    ap.add_argument("--lifted-probes", help="paf_lift.py points output for source probes")
    ap.add_argument("--qtl", required=True, help="QTL BED in source coords")
    ap.add_argument("--wga-intervals", help="paf_lift.py bed output")
    ap.add_argument("--miniprot-intervals", help="miniprot_lift.py output")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank", type=int, default=0,
                    help="pad transferred intervals by this many bp before scoring")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    src = {r["probe"]: (r["chrom"], int(r["pos"])) for r in read_tsv(args.probes_src)}
    tgt = {r["probe"]: (r["chrom"], int(r["pos"])) for r in read_tsv(args.probes_tgt)}
    truth = {p: (src[p], tgt[p]) for p in src if p in tgt}
    summary = []
    summary.append(f"probes unique in source      : {len(src)}")
    summary.append(f"probes unique in target      : {len(tgt)}")
    summary.append(f"TRUTH SET (unique in both)   : {len(truth)}")

    # ---------------- A. base-level accuracy of the WGA route -------------
    if args.lifted_probes:
        lifted = {r["name"]: r for r in read_tsv(args.lifted_probes)}
        n = wrong = 0
        errs = []
        with open(os.path.join(args.outdir, "probe_level_errors.tsv"), "w") as fo:
            fo.write("probe\tsrc_chrom\tsrc_pos\ttrue_tgt_chrom\ttrue_tgt_pos\t"
                     "lifted_chrom\tlifted_pos\tsame_chrom\terror_bp\n")
            for p, ((sc, sp), (tc, tp)) in truth.items():
                r = lifted.get(p)
                if r is None or r["tgt_chrom"] == "NA":
                    fo.write(f"{p}\t{sc}\t{sp}\t{tc}\t{tp}\tNA\tNA\tNA\tNA\n")
                    continue
                n += 1
                same = r["tgt_chrom"] == tc
                if not same:
                    wrong += 1
                    err = "NA"
                else:
                    e = abs(int(r["tgt_pos"]) - tp)
                    errs.append(e)
                    err = str(e)
                fo.write(f"{p}\t{sc}\t{sp}\t{tc}\t{tp}\t{r['tgt_chrom']}\t"
                         f"{r['tgt_pos']}\t{int(same)}\t{err}\n")

        N = max(len(truth), 1)
        summary.append("")
        summary.append("== A. Base-level accuracy (whole-genome alignment route) ==")
        summary.append(f"truth probes lifted          : {n} ({100*n/N:.1f}%)")
        summary.append(f"lifted to WRONG chromosome   : {wrong} ({100*wrong/max(n,1):.2f}%)")
        if errs:
            summary.append(f"median error (bp)            : {quant(errs,0.5)}")
            summary.append(f"mean error (bp)              : {sum(errs)/len(errs):.1f}")
            summary.append(f"P90 error (bp)               : {quant(errs,0.90)}")
            summary.append(f"P99 error (bp)               : {quant(errs,0.99)}")
            for thr in (100, 1000, 10000, 100000):
                k = sum(1 for e in errs if e <= thr)
                summary.append(f"within {thr:>7} bp             : {100*k/len(errs):.2f}%")

    # ---------------- B. interval-level, both routes ----------------------
    qtls = []
    with open(args.qtl) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            qtls.append((f[3] if len(f) > 3 else f"{f[0]}:{f[1]}-{f[2]}",
                         f[0], int(f[1]), int(f[2])))

    by_tgt_chrom = defaultdict(list)
    for p, (_, (tc, tp)) in truth.items():
        by_tgt_chrom[tc].append((tp, p))
    for c in by_tgt_chrom:
        by_tgt_chrom[c].sort()

    def score(table, label):
        rows = {r["qtl_id"]: r for r in read_tsv(table)}
        outp = os.path.join(args.outdir, f"interval_eval_{label}.tsv")
        rec_all, prec_all, jac_all = [], [], []
        n_ok = 0
        with open(outp, "w") as fo:
            fo.write("qtl_id\tsrc_chrom\tsrc_len\ttgt_chrom\ttgt_len\tlen_ratio\t"
                     "n_truth_probes\tn_recovered\trecall\tn_in_interval\tprecision\t"
                     "jaccard\tleft_boundary_err\tright_boundary_err\tflag\n")
            for qid, sc, s, e in qtls:
                inside = [p for p, ((c, sp), _) in truth.items()
                          if c == sc and s <= sp < e]
                r = rows.get(qid)
                if r is None or r["tgt_chrom"] == "NA":
                    fo.write(f"{qid}\t{sc}\t{e-s}\tNA\tNA\tNA\t{len(inside)}\t0\t0.000\t"
                             f"0\tNA\t0.000\tNA\tNA\t"
                             f"{r['flag'] if r else 'MISSING'}\n")
                    rec_all.append(0.0)
                    continue
                tc = r["tgt_chrom"]
                ts = int(r["tgt_start"]) - args.flank
                te = int(r["tgt_end"]) + args.flank
                rec_set = {p for p in inside
                           if truth[p][1][0] == tc and ts <= truth[p][1][1] < te}
                in_set = {p for tp, p in by_tgt_chrom.get(tc, []) if ts <= tp < te}
                recall = len(rec_set) / len(inside) if inside else float("nan")
                prec = len(rec_set) / len(in_set) if in_set else float("nan")
                union = len(set(inside) | in_set)
                jac = len(rec_set) / union if union else float("nan")
                # boundary error from the outermost truth probes
                truepos = sorted(truth[p][1][1] for p in inside
                                 if truth[p][1][0] == tc)
                lb = rb = "NA"
                if truepos:
                    lb = str(ts - truepos[0])       # <0 = interval starts too early
                    rb = str(te - truepos[-1])      # >0 = interval ends too late
                if inside:
                    rec_all.append(recall)
                if in_set:
                    prec_all.append(prec)
                if union:
                    jac_all.append(jac)
                n_ok += 1
                fo.write(f"{qid}\t{sc}\t{e-s}\t{tc}\t{int(r['tgt_len'])}\t"
                         f"{r['len_ratio']}\t{len(inside)}\t{len(rec_set)}\t{recall:.3f}\t"
                         f"{len(in_set)}\t{prec:.3f}\t{jac:.3f}\t{lb}\t{rb}\t{r['flag']}\n")

        summary.append("")
        summary.append(f"== B. Interval-level -- {label} ==")
        summary.append(f"QTLs transferred             : {n_ok}/{len(qtls)}")
        if rec_all:
            summary.append(f"mean probe recall            : {sum(rec_all)/len(rec_all):.3f}")
            summary.append(f"median probe recall          : {quant(rec_all,0.5):.3f}")
            summary.append(f"QTLs with recall == 1.0      : {sum(1 for x in rec_all if x >= 0.999)}")
        if prec_all:
            summary.append(f"mean probe precision         : {sum(prec_all)/len(prec_all):.3f}")
        if jac_all:
            summary.append(f"mean Jaccard                 : {sum(jac_all)/len(jac_all):.3f}")
        summary.append(f"per-QTL detail               : {outp}")

    if args.wga_intervals:
        score(args.wga_intervals, "minimap2_wga")
    if args.miniprot_intervals:
        score(args.miniprot_intervals, "miniprot")

    txt = "\n".join(summary) + "\n"
    print(txt, end="")
    open(os.path.join(args.outdir, "summary.txt"), "w").write(txt)


if __name__ == "__main__":
    main()
