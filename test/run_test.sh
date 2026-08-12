#!/usr/bin/env bash
# Integration test: runs the whole pipeline on synthetic data with a known
# answer (target = source + 1 kb insertion at 20,000, plus a 3 kb unalignable
# hole inside the second QTL). minimap2 and miniprot are replaced by shims,
# so no aligners need to be installed.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
NF="${NEXTFLOW:-nextflow}"

python3 "$HERE/make_test_data.py" "$HERE/data"
export PATH="$HERE/shims:$PATH"

cd "$ROOT"
"$NF" run . -profile test "$@"

echo; echo "=== route 1 (minimap2 WGA) ==="
cat "$HERE"/results/wga/*.qtl_target.tsv
echo; echo "=== route 1 segments ==="
cat "$HERE"/results/wga/*.qtl_segments.bed
echo; echo "=== route 2 (miniprot) ==="
cat "$HERE"/results/miniprot/*.qtl_target.tsv
echo; echo "=== benchmark ==="
cat "$HERE"/results/bench/summary.txt

python3 - "$HERE/results" <<'PY'
import csv, glob, sys
d = sys.argv[1]
rows = {r["qtl_id"]: r for r in
        csv.DictReader(open(glob.glob(f"{d}/wga/*.qtl_target.tsv")[0]), delimiter="\t")}
b = rows["QTL_before_indel"]
a = rows["QTL_after_indel"]
# interval before the insertion must not move
assert (int(b["tgt_start"]), int(b["tgt_end"])) == (5000, 15000), b
# interval after it must shift by exactly 1 kb
assert (int(a["tgt_start"]), int(a["tgt_end"])) == (26000, 36000), a
# both fully transferred in this scenario
assert b["frac_transferred"] == "1.000", b
assert a["frac_transferred"] == "1.000", a
print("\nPASS: intervals transferred to the expected coordinates")
PY
