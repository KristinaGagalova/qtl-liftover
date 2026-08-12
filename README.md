# qtl-liftover

Transfer QTL intervals from one variety's assembly to another by two
independent routes, and score both against your probe sequences.

| | Route 1 — `minimap2` | Route 2 — `miniprot` |
|---|---|---|
| What is aligned | whole genome vs whole genome (DNA) | source proteins vs target genome |
| Coordinate resolution | base-pair | gene-level (anchors) |
| Works in intergenic space | yes | no |
| Tolerates divergence | up to ~5–10% (`asm20`) | much higher — codon/aa level |
| Breaks on | large SVs, TE-rich gaps, unplaced scaffolds | gene loss, paralog expansion, bad annotation |
| Needs annotation | no | yes (source GFF + proteins) |

They fail for **different** reasons, which is why running both is worth it:
agreement is strong evidence, disagreement localises a structural variant or
a gene-family expansion sitting on your QTL.

## Quick start

```bash
nextflow run . \
    --source_genome   varietyA.fa \
    --target_genome   varietyB.fa \
    --source_gff      varietyA.gff3 \
    --source_proteins varietyA.pep.fa \
    --qtl_bed         qtl_intervals.bed \
    --probes          probes.fa \
    --outdir          ./results \
    -profile singularity,slurm
```

Source = the assembly your QTLs are currently defined on. `qtl_bed` is BED:
`chrom  start  end  qtl_id`, 0-based half-open.

## Pipeline steps

| Step | Module | Tool | Description |
|---|---|---|---|
| 1 | `MINIMAP2_ALIGN` | minimap2 | whole-genome alignment, source as query |
| 2 | `PAF_LIFT_BED` | python | lift intervals through the CIGAR by tiled anchors |
| 3 | `GFF2PROTEIN_BED` | python | representative transcript per gene + protein lengths |
| 4 | `MINIPROT_INDEX` | miniprot | index the target genome |
| 5 | `MINIPROT_ALIGN` | miniprot | spliced protein-to-genome alignment |
| 6 | `MINIPROT_LIFT` | python | anchors → transferred intervals + collinearity |
| 7 | `MINIMAP2_PROBES` | minimap2 | probes vs each assembly (2 parallel jobs) |
| 8 | `PROBE_POSITIONS` | python | keep confident, single-locus probes |
| 9 | `PAF_LIFT_POINTS` | python | lift source probe positions through the WGA PAF |
| 10 | `BENCHMARK_LIFTOVER` | python | score both routes against the truth set |

## Parameters

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--source_genome` | ✅ | - | assembly the QTLs are defined on |
| `--target_genome` | ✅ | - | assembly to move them to |
| `--qtl_bed` | ✅ | - | QTL intervals, BED, source coordinates |
| `--source_gff` | route 2 | - | source annotation |
| `--source_proteins` | route 2 | - | source proteins, IDs matching the GFF |
| `--probes` | benchmark | - | probe / marker sequences |
| `--outdir` | - | `./results` | output directory |
| `--asm_preset` | - | `asm10` | `asm5` ~0.1%, `asm10` ~1%, `asm20` ~5% divergence |
| `--minimap2_extra` | - | `-I 8G -K 4G` | raise `-I` above your genome size |
| `--n_anchor` | - | `200` | anchor points tiled per interval |
| `--miniprot_max_second` | - | `0.90` | lower = stricter about paralogs |
| `--bench_flank` | - | `0` | pad transferred intervals before scoring |
| `--skip_miniprot` | - | `false` | route 1 only, no annotation needed |
| `--skip_benchmark` | - | `false` | no probes needed |
| `--max_cpus` / `--max_memory` / `--max_time` | - | `16` / `128.GB` / `48.h` | resource ceilings |

## Profiles

| Profile | Description |
|---|---|
| `docker` | run with Docker containers |
| `singularity` | run with Singularity containers |
| `conda` | run with Conda environments |
| `slurm` | submit to a generic Slurm scheduler |
| `setonix` | Pawsey Setonix (Slurm + Singularity + scratch binds) |
| `test` | synthetic dataset, minimal resources |

Combine them: `-profile singularity,slurm`.

## Check this before anything else

Protein FASTA headers must match the `--gff_id_attr` values in the GFF, or
route 2 silently finds zero anchors:

```bash
grep '^>' varietyA.pep.fa | head -3 | cut -c2- | cut -d' ' -f1
awk '$3=="mRNA"' varietyA.gff3 | head -3 | sed 's/.*ID=\([^;]*\).*/\1/'
```

If they differ, regenerate the proteins from the GFF so they cannot disagree:

```bash
gffread varietyA.gff3 -g varietyA.fa -y varietyA.pep.fa -S
```

## Output structure

```
results/
├── wga/
│   ├── <source>.source_to_target.paf
│   ├── <source>.wga.qtl_target.tsv        route 1 transferred intervals
│   └── <source>.wga.qtl_segments.bed      the pieces each interval arrived in
├── miniprot/
│   ├── <source>.genes.bed
│   ├── <source>.miniprot.gff
│   ├── <source>.miniprot.qtl_target.tsv   route 2 transferred intervals
│   └── <source>.miniprot.anchors.tsv      per-gene anchor placements
├── probes/
│   ├── probes_vs_{src,tgt}.paf
│   ├── {src,tgt}.probe_pos.tsv            confident unique probe positions
│   └── probe.lifted.tsv
├── bench/
│   ├── summary.txt                        headline metrics
│   ├── probe_level_errors.tsv             per-probe bp error
│   └── interval_eval_*.tsv                per-QTL recall / precision
└── pipeline_info/
    ├── timeline.html  report.html  trace.txt  dag.html
```

## How much of each interval actually transferred

`wga/*.qtl_target.tsv` has one row per input QTL:

| Column | Meaning |
|---|---|
| `src_len` | length of the interval you supplied |
| `tgt_len` | length of the interval in the target assembly |
| `src_transferred_bp` | how many bp of the source interval found a home |
| `frac_transferred` | that as a fraction of `src_len` — the headline number |
| `largest_gap_bp` | longest unbroken stretch that did not transfer |
| `n_segments` | 1 = interval moved intact; >1 = it arrived in pieces |
| `len_ratio` | `tgt_len / src_len` — expansion or contraction |
| `anchor_step` | resolution of the bp figures (raise `--n_anchor` to sharpen) |

`tgt_len` and `src_transferred_bp` answer different questions. `tgt_len` is
the span from first to last anchor, gaps included — what you would screen.
`src_transferred_bp` excludes the gaps and says how much of the original
interval is genuinely accounted for.

`wga/*.qtl_segments.bed` breaks each interval into its contiguous transferred
pieces so you can see where the holes are.

## How the two routes work

**Route 1** (`LIFT_WGA`). minimap2 runs with the source as the *query*, so
the PAF encodes source → target. `PAF_LIFT_BED` does not assume the QTL fits
in one alignment block — that breaks on multi-Mb intervals. It tiles each
interval with `--n_anchor` points, lifts each through the CIGAR
independently, then takes the dominant target chromosome and the largest
positional cluster. An interval landing half on chr3 and half on chr7 is
flagged, not silently averaged.

**Route 2** (`LIFT_PROTEIN`). Genes annotated inside the QTL become anchors.
`--miniprot_extra '--outn=5'` keeps secondary hits deliberately: a protein
whose second-best hit nearly ties its best is a paralog and gets dropped
(`--miniprot_max_second`), which stops a tandem array from dragging the
interval across the chromosome. `collinearity_rho` is the Spearman
correlation of source vs target anchor order — below ~0.8 means rearranged,
not merely shifted.

## The benchmark

`BENCHMARK_PROBES` maps probes **independently** to both assemblies. Probes
placing uniquely and confidently in both give (source pos, target pos) pairs
that owe nothing to either route. Uniqueness filtering is not optional: array
probes are full of near-repeats, and a probe with two equally good hits
manufactures a fake liftover error. `-N 20 -p 0.4` is passed on purpose so
secondary hits can be counted.

*Base level* (route 1 only, since only it defines a coordinate function):
lift-over rate, wrong-chromosome rate, median/P90/P99 bp error, fraction
within 100 bp / 1 kb / 10 kb / 100 kb. Two cultivars of one species should
give >95% lifted, <0.5% wrong chromosome, single-digit bp median.

*Interval level* (both routes), per QTL:

- **recall** — probes inside the QTL in A that land inside the transferred
  interval in B. Recall < 1 means part of your QTL was thrown away and the
  causal gene may sit outside the new interval.
- **precision** — probes inside the transferred interval that were inside the
  original QTL. Low precision means inflation.
- **boundary errors** — signed bp from each transferred edge to the outermost
  true probe.

If recall falls short for a QTL you care about, re-run with
`--bench_flank 50000` to find how much padding restores full coverage.

**Caveat on the answer key.** Probes were designed on one genome, so they are
biased toward sequence that is present, single-copy and conserved — they
under-sample exactly the regions where liftover is hardest. The metrics
compare the two routes fairly and catch gross failures, but read them as a
lower bound on error. Probe count per QTL is reported so you can see when a
"perfect" recall rests on four markers.

## Reading the two answers together

| Situation | What it means | What to do |
|---|---|---|
| Both routes agree within a few kb | clean, syntenic region | use route 1's coordinates |
| Same chromosome, route 2 much narrower | route 2 is gene-bounded by construction | use route 1; route 2 confirms |
| Route 1 `FRAGMENTED_ON_CHROM`, route 2 clean | SV / repeat expansion inside the QTL | use route 2's span, inspect |
| Route 1 fine, route 2 `FAIL_too_few_anchors` | gene-poor interval or annotation gap | use route 1 |
| Different chromosomes | translocation, or a paralogous region hijacked one route | inspect before trusting either |
| Route 1 `INVERTED` / route 2 `rho < 0` | segment or chromosome flipped in the target | often real — assemblies get published in opposite orientations |

When the routes disagree but overlap, report the **union**.

## Config layering gotcha

Tool arguments live in `conf/modules.config` as `ext.args`; container images
and `publishDir` live in `conf/containers.config`, overriding the containers
declared inside each module.

**A `withName` block inside a profile-included config replaces that
selector's entire entry**, silently wiping the `ext.args` set in
`conf/modules.config` — the process then runs with no arguments and no error.
So in `conf/test.config`, `conf/slurm.config` and
`conf/profiles/*.config`, per-process resources are set with `withLabel`,
never `withName`. Top-level `withName` blocks in `nextflow.config` are fine;
those merge. If you must override a single process inside a profile,
re-declare its `ext.args` there too.

`check_max()` is defined in `nextflow.config` (the parent) so that every
`includeConfig`'d file can call it. Moving it into an included file puts it
out of scope for the parent.

## Practical notes

- Pick `--asm_preset` from measured divergence. `mash dist A.fa B.fa` takes
  seconds; ANI ~99.9% → `asm5`, ~99% → `asm10`, ~95% → `asm20`. Too tight a
  preset silently loses alignment in divergent regions.
- For large plant genomes raise `-I` in `--minimap2_extra` above the genome
  size, or minimap2 splits the index and the PAF covers only part of the
  target. `MINIMAP2_ALIGN` requests 96 GB by default for this reason.
- QTL intervals are usually defined by flanking *markers*. If you have those
  marker names, mapping their probes straight to the target and taking the
  span is the most defensible transfer of all — `probes/tgt.probe_pos.tsv`
  gives you that for free as a third opinion.
- Coordinates are 0-based half-open (BED). Add 1 to the start before handing
  results to anything expecting GFF/VCF conventions.
- Container tags in `conf/containers.config` are pinned; bump them as needed.

## Tests

```bash
./test/run_test.sh          # end-to-end on synthetic data, no aligners needed
nextflow run . -stub-run -profile test   # validate the DAG only
python3 bin/test_paf_lift.py             # base-exact CIGAR lifting, both strands
```

`test/run_test.sh` builds a target that is the source plus a 1 kb insertion
at position 20,000, runs the full pipeline with shims in place of minimap2
and miniprot, and asserts the QTL before the insertion does not move while
the one after it shifts by exactly 1 kb.

## Repository layout

```
main.nf                              thin entry point
nextflow.config                      params, profiles, check_max
nfcore_custom.config                 resource requirements by process label
conf/
  containers.config                  per-process container + publishDir overrides
  modules.config                     per-process ext.args
  test.config  slurm.config
  profiles/singularity.config  profiles/pawsey_setonix.config
workflows/qtl-liftover.nf            QTL_LIFTOVER_FLOW
subworkflows/local/
  lift_wga.nf                        route 1
  lift_protein.nf                    route 2
  benchmark_probes.nf                probe-based scoring
modules/local/                       10 processes
bin/                                 python scripts (auto-added to PATH)
test/                                data generator, tool shims, run_test.sh
assets/NO_FILE                       placeholder for skipped optional inputs
```

> **Note:** if the executable bit was lost in transfer, run
> `chmod +x bin/*.py test/shims/* test/run_test.sh` before the first run —
> Nextflow requires `bin/` scripts to be executable.
