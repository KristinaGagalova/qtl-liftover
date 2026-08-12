#!/usr/bin/env python3
"""
miniprot_scaffold_candidates.py -- when the target assembly has no
chromosome-scale sequences (thousands-to-millions of scaffolds), a single
"transferred interval" per QTL is a false precision: the true homologous
region is often genuinely split across several scaffolds, none of them
individually "correct". This reads the anchors.tsv already produced by
miniprot_lift.py and reports, per QTL, the ranked list of candidate
scaffolds instead of collapsing to one winner.

No rerun of miniprot needed -- this only re-aggregates the anchor table.

Input:  anchors.tsv from miniprot_lift.py --anchors-out
        columns: qtl_id  protein  src_chrom  src_mid  tgt_chrom  tgt_mid
                 identity  strand  used

Output:
  --out          one row per (QTL, candidate scaffold), ranked by anchor
                 count, with cumulative coverage of mapped anchors
  --summary-out  one row per QTL: scaffold count, concentration, a
                 CONCENTRATED / DIFFUSE call
"""

import argparse
import sys
from collections import defaultdict


def spearman(x, y):
    n = len(x)
    if n < 3:
        return float("nan")

    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True, nargs="+",
                    help="one or more anchors.tsv files (e.g. from separate runs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("--top-n", type=int, default=15,
                    help="max scaffolds reported per QTL (0 = all)")
    ap.add_argument("--min-anchors", type=int, default=2,
                    help="drop scaffolds with fewer than this many anchors")
    ap.add_argument("--concentrated-frac", type=float, default=0.40,
                    help="top scaffold's frac_of_mapped >= this -> CONCENTRATED")
    args = ap.parse_args()

    # qtl_id -> src_chrom, src_start/end (min/max src_mid seen)
    qtl_src = {}
    # (qtl_id, tgt_chrom) -> list of (src_mid, tgt_mid, strand)
    groups = defaultdict(list)
    n_rows = 0

    for path in args.anchors:
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(hdr)}
            need = ["qtl_id", "src_chrom", "src_mid", "tgt_chrom", "tgt_mid", "strand"]
            missing = [c for c in need if c not in idx]
            if missing:
                sys.exit(f"[miniprot_scaffold_candidates] {path}: missing columns {missing}")
            for line in fh:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                qid = f[idx["qtl_id"]]
                sc = f[idx["src_chrom"]]
                sm = int(f[idx["src_mid"]])
                tc = f[idx["tgt_chrom"]]
                tm = int(f[idx["tgt_mid"]])
                st = f[idx["strand"]]
                n_rows += 1

                if qid not in qtl_src:
                    qtl_src[qid] = [sc, sm, sm]
                else:
                    qtl_src[qid][1] = min(qtl_src[qid][1], sm)
                    qtl_src[qid][2] = max(qtl_src[qid][2], sm)

                groups[(qid, tc)].append((sm, tm, st))

    sys.stderr.write(f"[miniprot_scaffold_candidates] read {n_rows} anchors, "
                     f"{len(qtl_src)} QTLs, {len(groups)} (QTL, scaffold) groups\n")

    out = open(args.out, "w")
    out.write("qtl_id\tsrc_chrom\tsrc_span\trank\ttgt_scaffold\tn_anchors\t"
              "frac_of_mapped\tcum_frac\ttgt_start\ttgt_end\ttgt_span\t"
              "strand\tcollinearity_rho\n")

    summ = open(args.summary_out, "w")
    summ.write("qtl_id\tsrc_chrom\tsrc_span\tn_mapped_anchors\tn_scaffolds\t"
              "top_scaffold\ttop_scaffold_frac\ttop3_cum_frac\tcall\n")

    # per-QTL: gather its scaffold groups, rank, write
    by_qtl = defaultdict(list)
    for (qid, tc), rows in groups.items():
        by_qtl[qid].append((tc, rows))

    for qid, scaffolds in by_qtl.items():
        sc_chrom, smin, smax = qtl_src[qid]
        src_span = f"{sc_chrom}:{smin}-{smax}"
        n_mapped = sum(len(rows) for _, rows in scaffolds)

        scaffolds.sort(key=lambda x: -len(x[1]))
        kept = [(tc, rows) for tc, rows in scaffolds if len(rows) >= args.min_anchors]
        if args.top_n > 0:
            kept = kept[:args.top_n]

        cum = 0
        for rank, (tc, rows) in enumerate(kept, start=1):
            n_a = len(rows)
            cum += n_a
            frac = n_a / n_mapped if n_mapped else 0.0
            cum_frac = cum / n_mapped if n_mapped else 0.0
            tmins = [r[1] for r in rows]
            ts, te = min(tmins), max(tmins) + 1
            strand_counts = defaultdict(int)
            for r in rows:
                strand_counts[r[2]] += 1
            strand = max(strand_counts, key=strand_counts.get)
            rho = spearman([r[0] for r in rows], [r[1] for r in rows])
            rho_s = f"{rho:.3f}" if rho == rho else "nan"

            out.write(f"{qid}\t{sc_chrom}\t{src_span}\t{rank}\t{tc}\t{n_a}\t"
                      f"{frac:.3f}\t{cum_frac:.3f}\t{ts}\t{te}\t{te-ts}\t"
                      f"{strand}\t{rho_s}\n")

        n_scaf = len(scaffolds)
        if scaffolds:
            top_tc, top_rows = scaffolds[0]
            top_frac = len(top_rows) / n_mapped if n_mapped else 0.0
            top3_frac = sum(len(r[1]) for r in scaffolds[:3]) / n_mapped if n_mapped else 0.0
            call = "CONCENTRATED" if top_frac >= args.concentrated_frac else "DIFFUSE"
        else:
            top_tc, top_frac, top3_frac, call = "NA", 0.0, 0.0, "NO_ANCHORS"

        summ.write(f"{qid}\t{sc_chrom}\t{src_span}\t{n_mapped}\t{n_scaf}\t"
                  f"{top_tc}\t{top_frac:.3f}\t{top3_frac:.3f}\t{call}\n")

    out.close()
    summ.close()
    sys.stderr.write(f"[miniprot_scaffold_candidates] wrote {args.out} and "
                     f"{args.summary_out}\n")


if __name__ == "__main__":
    main()
