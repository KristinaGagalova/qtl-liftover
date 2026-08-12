#!/usr/bin/env python3
"""
miniprot_lift.py -- transfer QTL intervals via protein anchors.

Idea: the genes annotated inside the QTL in the SOURCE variety are the
anchors. Their proteins are mapped onto the TARGET assembly with miniprot;
the new interval is the span of the confidently, collinearly placed anchors.

Inputs
------
--gene-bed   source-genome positions of each protein: chrom start end protein_id
             (make it with gff_to_protein_bed.py)
--miniprot   miniprot --gff output of source proteins vs TARGET genome
--qtl        QTL intervals in source coordinates (BED)

Anchor filters
--------------
  * best hit per protein (lowest Rank / highest score)
  * Identity >= --min-ident, aligned protein fraction >= --min-cov
  * 1:1-ness: second-best hit score < --max-second * best score  (drops paralogs)
  * frameshifts / internal stops capped

Then: dominant target chromosome -> largest positional cluster ->
Spearman rho between source and target order as a collinearity QC.
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_lift import cluster  # noqa: E402

ATTR = re.compile(r"([A-Za-z_]+)=([^;]+)")


def parse_miniprot_gff(path):
    """Return {protein_id: [hit, ...]} from mRNA lines of miniprot --gff."""
    hits = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            attrs = dict(ATTR.findall(f[8]))
            tgt = attrs.get("Target", "")
            parts = tgt.split()
            if not parts:
                continue
            pid = parts[0]
            plen_aln = None
            if len(parts) >= 3:
                try:
                    plen_aln = int(parts[2]) - int(parts[1]) + 1
                except ValueError:
                    pass
            hits[pid].append(dict(
                chrom=f[0], start=int(f[3]) - 1, end=int(f[4]),
                score=float(f[5]) if f[5] not in (".", "") else 0.0,
                strand=f[6],
                rank=int(attrs.get("Rank", 1)),
                ident=float(attrs.get("Identity", 0)),
                positive=float(attrs.get("Positive", 0)),
                fs=int(attrs.get("Frameshift", 0)),
                stop=int(attrs.get("StopCodon", 0)),
                aln_len=plen_aln,
            ))
    return hits


def best_anchors(hits, plen, min_ident, min_cov, max_second, max_fs, max_stop):
    """hits -> {pid: hit} keeping only confident, near-unique placements."""
    keep, stats = {}, Counter()
    for pid, hl in hits.items():
        hl = sorted(hl, key=lambda h: (-h["score"], h["rank"]))
        b = hl[0]
        stats["total"] += 1
        if b["ident"] < min_ident:
            stats["low_identity"] += 1
            continue
        L = plen.get(pid)
        if L and b["aln_len"] and b["aln_len"] / L < min_cov:
            stats["low_coverage"] += 1
            continue
        if b["fs"] > max_fs or b["stop"] > max_stop:
            stats["frameshift_or_stop"] += 1
            continue
        if len(hl) > 1 and hl[1]["score"] > max_second * b["score"]:
            stats["multi_locus"] += 1
            continue
        keep[pid] = b
        stats["kept"] += 1
    return keep, stats


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
    ap.add_argument("--gene-bed", required=True)
    ap.add_argument("--miniprot", required=True)
    ap.add_argument("--qtl", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--anchors-out", help="per-anchor TSV (useful for plots/QC)")
    ap.add_argument("--prot-len", help="optional TSV protein_id<TAB>length for coverage filter")
    ap.add_argument("--min-ident", type=float, default=0.80)
    ap.add_argument("--min-cov", type=float, default=0.70)
    ap.add_argument("--max-second", type=float, default=0.90)
    ap.add_argument("--max-fs", type=int, default=2)
    ap.add_argument("--max-stop", type=int, default=1)
    ap.add_argument("--min-anchors", type=int, default=3)
    ap.add_argument("--max-gap", type=int, default=2_000_000)
    args = ap.parse_args()

    genes = []  # (chrom, start, end, pid)
    with open(args.gene_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            genes.append((f[0], int(f[1]), int(f[2]), f[3]))
    by_chrom = defaultdict(list)
    for g in genes:
        by_chrom[g[0]].append(g)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda g: g[1])

    plen = {}
    if args.prot_len:
        with open(args.prot_len) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    plen[f[0]] = int(f[1])

    hits = parse_miniprot_gff(args.miniprot)
    anchors, stats = best_anchors(hits, plen, args.min_ident, args.min_cov,
                                  args.max_second, args.max_fs, args.max_stop)
    sys.stderr.write("[miniprot_lift] anchor filtering: " +
                     ", ".join(f"{k}={v}" for k, v in stats.items()) + "\n")

    afh = open(args.anchors_out, "w") if args.anchors_out else None
    if afh:
        afh.write("qtl_id\tprotein\tsrc_chrom\tsrc_mid\ttgt_chrom\ttgt_mid\t"
                  "identity\tstrand\tused\n")

    out = open(args.out, "w")
    out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tsrc_len\ttgt_chrom\t"
              "tgt_start\ttgt_end\ttgt_len\tlen_ratio\tstrand\tn_genes\t"
              "n_anchor_mapped\tn_anchor_used\tfrac_dominant_chrom\t"
              "collinearity_rho\tflag\n")

    with open(args.qtl) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            qid = f[3] if len(f) > 3 else f"{chrom}:{s}-{e}"
            L = e - s
            inside = [g for g in by_chrom.get(chrom, []) if g[1] < e and g[2] > s]
            mapped = [(g, anchors[g[3]]) for g in inside if g[3] in anchors]

            if len(mapped) < args.min_anchors:
                out.write(f"{qid}\t{chrom}\t{s}\t{e}\t{L}\tNA\tNA\tNA\tNA\tNA\tNA\t"
                          f"{len(inside)}\t{len(mapped)}\t0\tNA\tNA\t"
                          f"FAIL_too_few_anchors\n")
                continue

            cnt = Counter(h["chrom"] for _, h in mapped)
            tchrom, n_dom = cnt.most_common(1)[0]
            dom = [(g, h) for g, h in mapped if h["chrom"] == tchrom]
            mids = [(h["start"] + h["end"]) // 2 for _, h in dom]
            keep = set(cluster(mids, max(args.max_gap, 2 * L)))
            used = [(g, h, m) for (g, h), m in zip(dom, mids) if m in keep]

            ts = min(min(h["start"] for _, h, _ in used), min(m for *_, m in used))
            te = max(h["end"] for _, h, _ in used)
            tl = te - ts
            rho = spearman([(g[1] + g[2]) // 2 for g, _, _ in used],
                           [m for *_, m in used])
            strand = Counter(h["strand"] for _, h, _ in used).most_common(1)[0][0]

            if afh:
                usedset = {g[3] for g, _, _ in used}
                for g, h in mapped:
                    afh.write(f"{qid}\t{g[3]}\t{g[0]}\t{(g[1]+g[2])//2}\t{h['chrom']}\t"
                              f"{(h['start']+h['end'])//2}\t{h['ident']:.3f}\t"
                              f"{h['strand']}\t{int(g[3] in usedset)}\n")

            flags = []
            frac_dom = n_dom / len(mapped)
            if frac_dom < 0.9:
                flags.append("SPLIT_ACROSS_CHROMS")
            if len(used) < args.min_anchors:
                flags.append("FEW_USED_ANCHORS")
            if rho == rho and abs(rho) < 0.8:
                flags.append("POOR_COLLINEARITY")
            if rho == rho and rho < 0:
                flags.append("INVERTED")
            if not (0.5 <= tl / L <= 2.0):
                flags.append("LENGTH_CHANGE")
            flag = ";".join(flags) if flags else "PASS"

            out.write(f"{qid}\t{chrom}\t{s}\t{e}\t{L}\t{tchrom}\t{ts}\t{te}\t{tl}\t"
                      f"{tl/L:.3f}\t{strand}\t{len(inside)}\t{len(mapped)}\t{len(used)}\t"
                      f"{frac_dom:.2f}\t{rho:.3f}\t{flag}\n")
    out.close()
    if afh:
        afh.close()


if __name__ == "__main__":
    main()
