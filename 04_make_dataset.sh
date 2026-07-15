#!/usr/bin/env bash
# Stage 4 — turn HybriDetector hybrids into a labelled, miRBench-style v7 dataset.
#
#   bash chimeric_eclip/04_make_dataset.sh external/HybriDetector/hyb_pairs data/newset
#
# Runs miRBench's 7-step post-processing (their run_postprocess_pipeline.sh, single
# mode) after concatenating the per-sample hybrid tables:
#
#   0 filter+dedup      keep noncodingRNA_type == miRNA, rename to the v7 columns,
#                       drop duplicate (gene, noncodingRNA) sequence pairs
#   1 annotate          overlap each site against Ensembl-90 transcripts -> feature,
#                       dominant_region, regions_present, transcript coordinates
#   2 exclude families  hold out miRNA families that occur nowhere else (leftout set)
#   3 make negatives    THE key step: for each positive, sample a binding site from a
#                       *different gene cluster*, keeping each miRNA family's share of
#                       the negative class equal to its share of the positive class.
#                       That is the frequency-class-bias fix the paper is named for —
#                       a model cannot win by memorising which miRNAs are common.
#   4 train/test split  on the `test` column, which step 0 sets from chr.g == "1":
#                       chromosome 1 is the held-out test set, everything else trains
#   5 drop test col
#   6 add conservation  gene_phyloP + gene_phastCons per site, from the bigwigs
#
# The result has exactly the 18 columns of the existing *_v7.tsv files, so it drops
# straight into cnn/run_train.sh as TRAIN=/TESTS=.
#
# Env vars: NEG_RATIO (negatives per positive, default 1), PHYLOP_BW, PHASTCONS_BW.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env mirbench_pp

HYB_DIR="${1:?usage: 04_make_dataset.sh <HybriDetector hyb_pairs dir> <output dir>}"
OUT_DIR="${2:?usage: 04_make_dataset.sh <HybriDetector hyb_pairs dir> <output dir>}"
NAME="${NAME:-$(basename "$OUT_DIR")}"

for bw in "$PHYLOP_BW" "$PHASTCONS_BW"; do
    if [ ! -s "$bw" ]; then
        echo "Error: conservation bigwig not found: $bw" >&2
        echo "  hg38.phastCons100way.bw: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/" >&2
        echo "  hg38.phyloP100way.bw:    https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/" >&2
        exit 1
    fi
done

mkdir -p "$OUT_DIR/input"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
CONCAT="$OUT_DIR/input/$NAME.tsv"

# Concatenate the per-sample hybrid tables, keeping only the first header. (Upstream's
# concat_HD_output.sh does exactly this; inlined to avoid its SLURM preamble.)
if [ ! -s "$CONCAT" ]; then
    echo "== concatenating hybrid tables from $HYB_DIR"
    first=1
    for f in "$HYB_DIR"/*unified_length_all_types_unique_high_confidence.tsv; do
        echo "   + $(basename "$f")"
        if [ "$first" = 1 ]; then cat "$f" > "$CONCAT"; first=0
        else tail -n +2 "$f" >> "$CONCAT"; fi
    done
fi
echo "== $(( $(wc -l < "$CONCAT") - 1 )) hybrid rows"

echo "== running miRBench post-processing (7 steps)"
cd "$OUT_DIR"
mm_run mirbench_pp bash "$MIRBENCH_DIR/code/post_process/run_postprocess_pipeline.sh" \
    --mode single \
    -f "$CONCAT" \
    -o "$OUT_DIR" \
    -p "$PHYLOP_BW" \
    -c "$PHASTCONS_BW"

echo
echo "Done. Final datasets (v7 schema) in $OUT_DIR/step6_add_conservation/:"
ls -la "$OUT_DIR"/step6_add_conservation/*.tsv 2>/dev/null
echo
echo "Train the CNN on them with e.g.:"
echo "  TRAIN=<...train...conservation.tsv> TESTS=<...test...conservation.tsv> bash cnn/run_train.sh kfold"
