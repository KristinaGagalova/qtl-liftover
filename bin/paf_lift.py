#!/usr/bin/env python3
"""
paf_lift.py -- lift coordinates from a SOURCE assembly (= PAF query) to a
TARGET assembly (= PAF target) using a minimap2 whole-genome alignment.

The PAF must be produced with the source genome as the QUERY, e.g.

    minimap2 -cx asm10 --cs -t 16 target.fa source.fa > source_to_target.paf

Sub-commands
------------
  points   lift single positions            (TSV: chrom  pos  [name])
  bed      lift intervals by tiling anchors (BED: chrom start end [name])

Interval lifting does NOT require the whole QTL to sit in one alignment
block: the interval is tiled with anchor points, every anchor is lifted
independently, and a robust consensus (dominant target chromosome +
largest positional cluster) defines the new interval. That is what you
want for multi-Mb QTL intervals crossing SVs and assembly gaps.
"""

import argparse
import bisect
import re
import sys
from collections import Counter, defaultdict

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


class Aln:
    __slots__ = ("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen",
                 "tstart", "tend", "nmatch", "alnlen", "mapq", "cigar", "tp")

    def __init__(self, f):
        self.qname = f[0]
        self.qlen = int(f[1])
        self.qstart = int(f[2])
        self.qend = int(f[3])
        self.strand = f[4]
        self.tname = f[5]
        self.tlen = int(f[6])
        self.tstart = int(f[7])
        self.tend = int(f[8])
        self.nmatch = int(f[9])
        self.alnlen = int(f[10])
        self.mapq = int(f[11])
        self.cigar = None
        self.tp = "P"
        for tag in f[12:]:
            if tag.startswith("cg:Z:"):
                self.cigar = tag[5:]
            elif tag.startswith("tp:A:"):
                self.tp = tag[5:]

    @property
    def ident(self):
        return self.nmatch / self.alnlen if self.alnlen else 0.0


def lift_point(aln, qpos):
    """Walk the CIGAR of one alignment and translate query pos -> target pos.

    Returns (tpos, exact) or None. `exact` is False when the query base falls
    inside an insertion (no target base exists; nearest target pos returned).
    PAF convention: I consumes query, D/N consumes target.
    """
    if aln.cigar is None:
        return None
    if not (aln.qstart <= qpos < aln.qend):
        return None

    t = aln.tstart
    fwd = aln.strand == "+"
    q = aln.qstart if fwd else aln.qend

    for m in CIGAR_RE.finditer(aln.cigar):
        ln = int(m.group(1))
        op = m.group(2)
        if op in "M=X":
            if fwd:
                if q <= qpos < q + ln:
                    return (t + (qpos - q), True)
                q += ln
            else:
                if q - ln <= qpos < q:
                    return (t + (q - 1 - qpos), True)
                q -= ln
            t += ln
        elif op == "I":
            if fwd:
                if q <= qpos < q + ln:
                    return (t, False)
                q += ln
            else:
                if q - ln <= qpos < q:
                    return (t, False)
                q -= ln
        elif op in "DN":
            t += ln
        elif op in "SH":
            continue
        else:
            continue
    return None


class PafIndex:
    def __init__(self, path, min_mapq=5, min_alnlen=1000, primary_only=True,
                 min_ident=0.0):
        self.by_q = defaultdict(list)
        n_tot = n_kept = 0
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                n_tot += 1
                a = Aln(f)
                if a.cigar is None:
                    continue
                if primary_only and a.tp not in ("P", "I"):
                    continue
                if a.mapq < min_mapq or a.alnlen < min_alnlen:
                    continue
                if a.ident < min_ident:
                    continue
                self.by_q[a.qname].append(a)
                n_kept += 1
        self.starts = {}
        for q, lst in self.by_q.items():
            lst.sort(key=lambda a: a.qstart)
            self.starts[q] = [a.qstart for a in lst]
        sys.stderr.write(f"[paf_lift] loaded {n_kept}/{n_tot} alignment blocks "
                         f"over {len(self.by_q)} source sequences\n")

    def lift(self, chrom, pos):
        """Return dict with target coords, or None."""
        lst = self.by_q.get(chrom)
        if not lst:
            return None
        i = bisect.bisect_right(self.starts[chrom], pos)
        best = None
        # scan backwards; blocks are sorted by qstart, overlaps are local
        for a in reversed(lst[:i]):
            if a.qend <= pos:
                # keep scanning a little: blocks may be nested/overlapping
                if best is not None and a.qstart < pos - 10_000_000:
                    break
                continue
            r = lift_point(a, pos)
            if r is None:
                continue
            cand = dict(tname=a.tname, tpos=r[0], exact=r[1], strand=a.strand,
                        mapq=a.mapq, alnlen=a.alnlen, ident=a.ident)
            if best is None or (cand["alnlen"], cand["mapq"]) > (best["alnlen"], best["mapq"]):
                best = cand
        return best


def cluster(positions, max_gap):
    """Largest single-linkage cluster of sorted positions."""
    positions = sorted(positions)
    best, cur = [], [positions[0]]
    for p in positions[1:]:
        if p - cur[-1] <= max_gap:
            cur.append(p)
        else:
            if len(cur) > len(best):
                best = cur
            cur = [p]
    if len(cur) > len(best):
        best = cur
    return best


def cmd_points(args):
    idx = PafIndex(args.paf, args.min_mapq, args.min_alnlen, not args.allow_secondary)
    out = open(args.out, "w")
    out.write("name\tsrc_chrom\tsrc_pos\ttgt_chrom\ttgt_pos\tstrand\texact\tmapq\n")
    n = ok = 0
    with open(args.tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            chrom, pos = f[0], int(f[1])
            name = f[2] if len(f) > 2 else f"{chrom}:{pos}"
            n += 1
            r = idx.lift(chrom, pos)
            if r is None:
                out.write(f"{name}\t{chrom}\t{pos}\tNA\tNA\tNA\tNA\tNA\n")
                continue
            ok += 1
            out.write(f"{name}\t{chrom}\t{pos}\t{r['tname']}\t{r['tpos']}\t"
                      f"{r['strand']}\t{int(r['exact'])}\t{r['mapq']}\n")
    out.close()
    sys.stderr.write(f"[paf_lift] lifted {ok}/{n} points ({100*ok/max(n,1):.1f}%)\n")


def cmd_bed(args):
    idx = PafIndex(args.paf, args.min_mapq, args.min_alnlen, not args.allow_secondary)
    out = open(args.out, "w")
    out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tsrc_len\ttgt_chrom\t"
              "tgt_start\ttgt_end\ttgt_len\tlen_ratio\tstrand\t"
              "src_transferred_bp\tfrac_transferred\tlargest_gap_bp\tn_segments\t"
              "anchor_step\tn_anchor\tn_lifted\tfrac_dominant_chrom\tfrac_in_cluster\tflag\n")
    seg_out = open(args.segments_out, "w") if args.segments_out else None
    if seg_out:
        seg_out.write("qtl_id\tsrc_chrom\tsrc_start\tsrc_end\tsrc_bp\t"
                      "tgt_chrom\ttgt_start\ttgt_end\ttgt_bp\tstrand\tn_anchor\n")

    with open(args.bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            name = f[3] if len(f) > 3 else f"{chrom}:{s}-{e}"
            L = e - s
            step = max(args.min_step, L // max(args.n_anchor - 1, 1))
            anchors = sorted(set(list(range(s, e, step)) + [max(e - 1, s)]))
            n_a = len(anchors)

            lifts = [(p, idx.lift(chrom, p)) for p in anchors]
            hits = [h for _, h in lifts if h]
            n_l = len(hits)
            if n_l == 0:
                out.write(f"{name}\t{chrom}\t{s}\t{e}\t{L}\tNA\tNA\tNA\tNA\tNA\tNA\t"
                          f"0\t0.000\t{L}\t0\t{step}\t{n_a}\t0\tNA\tNA\t"
                          f"FAIL_no_anchor_lifted\n")
                continue

            cnt = Counter(h["tname"] for h in hits)
            tchrom, n_dom = cnt.most_common(1)[0]
            dom_pos = [h["tpos"] for h in hits if h["tname"] == tchrom]
            max_gap = max(args.max_gap, 2 * L)
            keep = set(cluster(dom_pos, max_gap))

            # which anchors survived chromosome vote + clustering
            ok = [(p, h) for p, h in lifts
                  if h and h["tname"] == tchrom and h["tpos"] in keep]
            ok_pos = {p for p, _ in ok}
            n_ok = len(ok)

            ts = min(h["tpos"] for _, h in ok)
            te = max(h["tpos"] for _, h in ok) + 1
            tl = te - ts
            strand = Counter(h["strand"] for _, h in ok).most_common(1)[0][0]

            # contiguous runs of transferred anchors -> segments of the interval
            segments, run = [], []
            for p in anchors:
                if p in ok_pos:
                    run.append(p)
                elif run:
                    segments.append(run)
                    run = []
            if run:
                segments.append(run)

            # transferred bp = union of the segment spans (resolution = anchor_step)
            seg_spans = [(r[0], min(r[-1] + step, e)) for r in segments]
            transferred = sum(b - a for a, b in seg_spans)
            frac_tr = transferred / L if L else 0.0
            # longest stretch of consecutive failed anchors
            gap = cur = 0
            for p in anchors:
                cur = 0 if p in ok_pos else cur + 1
                gap = max(gap, cur)
            largest_gap = min(gap * step, L)

            if seg_out:
                lut = dict(ok)
                for run, (a, b) in zip(segments, seg_spans):
                    tp = [lut[p]["tpos"] for p in run]
                    seg_out.write(
                        f"{name}\t{chrom}\t{a}\t{b}\t{b - a}\t{tchrom}\t{min(tp)}\t"
                        f"{max(tp) + 1}\t{max(tp) + 1 - min(tp)}\t{strand}\t{len(run)}\n")

            frac_dom = n_dom / n_l
            frac_clu = len(keep) / n_l

            flags = []
            if frac_tr < args.warn_lifted:
                flags.append("PARTIAL_TRANSFER")
            if frac_dom < args.warn_dom:
                flags.append("SPLIT_ACROSS_CHROMS")
            if frac_clu < args.warn_dom:
                flags.append("FRAGMENTED_ON_CHROM")
            if len(segments) > 1:
                flags.append(f"GAPPED_x{len(segments)}")
            if not (0.5 <= tl / L <= 2.0):
                flags.append("LENGTH_CHANGE")
            if strand == "-":
                flags.append("INVERTED")
            flag = ";".join(flags) if flags else "PASS"

            out.write(f"{name}\t{chrom}\t{s}\t{e}\t{L}\t{tchrom}\t{ts}\t{te}\t{tl}\t"
                      f"{tl/L:.3f}\t{strand}\t{transferred}\t{frac_tr:.3f}\t"
                      f"{largest_gap}\t{len(segments)}\t{step}\t{n_a}\t{n_l}\t"
                      f"{frac_dom:.2f}\t{frac_clu:.2f}\t{flag}\n")
    out.close()
    if seg_out:
        seg_out.close()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--paf", required=True, help="source(query) -> target PAF with cg tag")
    common.add_argument("--out", required=True)
    common.add_argument("--min-mapq", type=int, default=5)
    common.add_argument("--min-alnlen", type=int, default=1000)
    common.add_argument("--allow-secondary", action="store_true")

    a = sub.add_parser("points", parents=[common])
    a.add_argument("--tsv", required=True, help="chrom<TAB>pos[<TAB>name], 0-based pos")
    a.set_defaults(func=cmd_points)

    b = sub.add_parser("bed", parents=[common])
    b.add_argument("--bed", required=True, help="QTL intervals in source coordinates")
    b.add_argument("--n-anchor", type=int, default=200, help="anchors per interval")
    b.add_argument("--min-step", type=int, default=1000)
    b.add_argument("--max-gap", type=int, default=1_000_000,
                   help="max gap when clustering lifted anchors on target chrom")
    b.add_argument("--segments-out",
                   help="BED of the contiguous pieces each interval transferred as")
    b.add_argument("--warn-lifted", type=float, default=0.5,
                   help="flag PARTIAL_TRANSFER below this transferred fraction")
    b.add_argument("--warn-dom", type=float, default=0.9)
    b.set_defaults(func=cmd_bed)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
