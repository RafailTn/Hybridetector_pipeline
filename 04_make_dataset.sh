#!/usr/bin/env bash
# Stage 4 — turn HybriDetector hybrids into a labelled, miRBench-style v7 dataset.
#
#   bash chimeric_eclip/04_make_dataset.sh external/HybriDetector/hyb_pairs data/newset
#
# WHICH SAMPLES GO IN: HybriDetector hardwires all of its per-sample hybrid tables into
# one shared external/HybriDetector/hyb_pairs/ directory, which ACCUMULATES across every
# stage-3 run you ever do. Globbing the whole directory would therefore silently fuse
# unrelated experiments/fractions (e.g. GSE297116 cytoplasm + chromatin + Manakov) into a
# single dataset. So pass one or more sample-name globs to pick exactly the samples for
# THIS dataset — matched against the filename prefix, i.e. the stage-3 <sample> name:
#
#   bash chimeric_eclip/04_make_dataset.sh external/HybriDetector/hyb_pairs \
#        data/cytoplasm 'HCT116_cytoplasm_*'
#
# With no globs it takes every table in the directory (back-compatible, and correct when
# hyb_pairs holds a single experiment) — but it lists them so a stray fraction is obvious.
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
#                       a model cannot win by memorising which miRNAs are common. The
#                       ratio is a FIXED 1 negative per positive (num_neg = block size in
#                       miRBench's make_neg_sets.py) — the sampler exposes no ratio knob.
#   4 train/test split  on the `test` column, which step 0 sets from chr.g == "1":
#                       chromosome 1 is the held-out test set, everything else trains
#   5 drop test col
#   6 add conservation  gene_phyloP + gene_phastCons per site, from the bigwigs
#
# The result has exactly the 18 columns of the existing *_v7.tsv files, so it drops
# straight into cnn/run_train.sh as TRAIN=/TESTS=.
#
# Env vars: PHYLOP_BW, PHASTCONS_BW.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env mirbench_pp

HYB_DIR="${1:?usage: 04_make_dataset.sh <HybriDetector hyb_pairs dir> <output dir> [sample-glob ...]}"
OUT_DIR="${2:?usage: 04_make_dataset.sh <HybriDetector hyb_pairs dir> <output dir> [sample-glob ...]}"
shift 2
# Remaining args = sample-name globs selecting which hybrid tables to include; none = all.
PATTERNS=("$@")
[ "${#PATTERNS[@]}" -gt 0 ] || PATTERNS=("*")
NAME="${NAME:-$(basename "$OUT_DIR")}"

HYB_SUFFIX="unified_length_all_types_unique_high_confidence.tsv"

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

# Resolve the selected sample globs to a deduplicated list of hybrid tables. Each glob is
# matched against "<sample>.<suffix>", so 'HCT116_cytoplasm_*' picks only that fraction and
# never a chromatin/whole-cell table sitting in the same shared directory.
FILES=()
for pat in "${PATTERNS[@]}"; do
    for f in "$HYB_DIR"/$pat.$HYB_SUFFIX; do
        [ -e "$f" ] || continue
        case " ${FILES[*]} " in *" $f "*) continue ;; esac   # dedupe overlapping globs
        FILES+=("$f")
    done
done
if [ "${#FILES[@]}" -eq 0 ]; then
    echo "Error: no hybrid tables in $HYB_DIR match: ${PATTERNS[*]}" >&2
    echo "  available samples:" >&2
    for f in "$HYB_DIR"/*."$HYB_SUFFIX"; do
        [ -e "$f" ] && echo "    $(basename "$f" ".$HYB_SUFFIX")" >&2
    done
    exit 1
fi

# Concatenate the per-sample hybrid tables, keeping only the first header. (Upstream's
# concat_HD_output.sh does exactly this; inlined to avoid its SLURM preamble.)
#
# The concatenation is cached so an identical re-run (e.g. resuming after a failed postprocess
# step) need not redo it — but the cache is keyed on a MANIFEST of the selected inputs (path +
# size + mtime), not merely on $CONCAT existing. So changing the sample globs, or editing or
# replacing the hybrid tables, reuses the same OUT_DIR but rebuilds the concatenation instead
# of silently serving the previous selection's file. Sorted so pattern order is irrelevant.
MANIFEST="$CONCAT.manifest"
new_manifest=$(for f in "${FILES[@]}"; do
    printf '%s\t%s\t%s\n' "$f" "$(file_size "$f")" "$(file_mtime "$f")"
done | sort)

if [ -s "$CONCAT" ] && [ -f "$MANIFEST" ] && [ "$(cat "$MANIFEST")" = "$new_manifest" ]; then
    echo "== reusing cached concatenation ($CONCAT): ${#FILES[@]} input(s) unchanged"
else
    echo "== concatenating ${#FILES[@]} hybrid table(s) from $HYB_DIR"
    first=1
    for f in "${FILES[@]}"; do
        echo "   + $(basename "$f")"
        if [ "$first" = 1 ]; then cat "$f" > "$CONCAT"; first=0
        else tail -n +2 "$f" >> "$CONCAT"; fi
    done
    printf '%s\n' "$new_manifest" > "$MANIFEST"   # only after a complete concat
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
