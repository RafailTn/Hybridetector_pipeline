#!/usr/bin/env bash
# Stage 3, on a SLURM cluster — with a DENSE STAR index (upstream defaults).
#
#   bash chimeric_eclip/03_run_hybridetector_slurm.sh data/raw/gse297116
#
# This is the same workflow as 03_run_hybridetector.sh, minus the compromises that
# 31 GB of RAM forced locally:
#
#   local (03_run_hybridetector.sh)   cluster (this script)
#   ------------------------------    ---------------------------------
#   STAR_SA_SPARSE_D=2  (sparse SA)   dense — STAR's default, and what miRBench used
#                                     for Manakov, so mapping is directly comparable
#   HD_MEM_* = 28/24/24               upstream 200/50/34, passed to sbatch as --mem
#   -j <nproc>, one machine           one sbatch job per rule, up to `jobs:` at once
#
# WHY THIS MATTERS: a dense hg38 index needs ~30 GB resident to build AND to align
# against. On the 31 GB laptop it was OOM-killed at 27.9 GB, so sparse was the only
# option there — and sparse is NOT alignment-neutral (6931 differing reads, 5.7x the
# 1224-read dense-vs-dense noise floor on a chr21 control). On a cluster no reason to accept that,
# so this script deliberately does not export STAR_SA_SPARSE_D.
#
# BEFORE FIRST RUN, edit chimeric_eclip/slurm/config.yaml:
#   - `partition=compute` -> your cluster's partition name
#   - the per-rule mem/time in `set-resources` if your scheduler is stricter
#
# Env vars: CORES (threads per rule), READ_LENGTH, IS_UMI, JOBS, EXTRA_SNAKE_ARGS.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env hybridetector

PP_DIR="${1:?usage: 03_run_hybridetector_slurm.sh <fastq dir>}"
PROFILE="${PROFILE:-$PIPELINE_DIR/slurm}"
JOBS="${JOBS:-20}"
READ_LENGTH="${READ_LENGTH:-auto}"
IS_UMI="${IS_UMI:-TRUE}"
CORES="${CORES:-16}"

command -v sbatch >/dev/null || { echo "Error: sbatch not found — this is the cluster runner. Use 03_run_hybridetector.sh locally." >&2; exit 1; }
grep -q "partition=compute" "$PROFILE/config.yaml" && \
    echo "WARNING: $PROFILE/config.yaml still says partition=compute — set your real partition." >&2

# Deliberately NOT set: STAR_SA_SPARSE_D, HD_MEM_*. Upstream defaults = dense + 200/50/34.
# This runner is the DENSE one, so refuse to silently reuse a leftover SPARSE genome index
# (e.g. one built earlier on a small local box) — Snakemake would not rebuild it on its own.
check_star_index_sparsity 1

mkdir -p "$HD_DIR/data" "$HD_DIR/slurm_logs"
FQS=()
while IFS= read -r fq; do FQS+=("$fq"); done < <(find "$PP_DIR" -name "*.pp.fastq.gz" | sort)
SUFFIX=".pp.fastq.gz"
if [ "${#FQS[@]}" -eq 0 ]; then
    while IFS= read -r fq; do FQS+=("$fq"); done < <(find "$PP_DIR" -maxdepth 1 -name "*.fastq.gz" | sort)
    SUFFIX=".fastq.gz"
fi
[ "${#FQS[@]}" -gt 0 ] || { echo "Error: no fastq found under $PP_DIR" >&2; exit 1; }

SAMPLES=()
for fq in "${FQS[@]}"; do
    sample="$(basename "$fq" "$SUFFIX")"
    ln -sfn "$(abspath "$fq")" "$HD_DIR/data/$sample.fastq.gz"
    SAMPLES+=("$sample")
done
echo "== ${#SAMPLES[@]} sample(s): ${SAMPLES[*]}"

# Per-sample read_length (see the long note in 03_run_hybridetector.sh): a batch can mix
# libraries of different length, and HybriDetector's min miRNA fraction is derived per sample,
# so a single global max would be too lenient for the shorter samples.
READ_LENGTHS=()
while IFS= read -r rl; do READ_LENGTHS+=("$rl"); done \
    < <(read_lengths_for eclip_pp "$READ_LENGTH" "${FQS[@]}")
if [ "$READ_LENGTH" = "auto" ]; then
    echo "== measured per-sample read_length (post-trim max):"
    for i in "${!SAMPLES[@]}"; do echo "     ${SAMPLES[$i]} = ${READ_LENGTHS[$i]}"; done
else
    echo "== read_length (manual override, all samples) = $READ_LENGTH"
fi
RL_JOINED=$(printf '"%s",' "${READ_LENGTHS[@]}" | sed 's/,$//')

CONFIG="$HD_DIR/config.json"
{
    printf '{"Sample":['
    printf '"%s",' "${SAMPLES[@]}" | sed 's/,$//'
    printf '],\n'
    printf ' "map_perc_single_genomic":[%s],\n' "$(printf '"0.85",%.0s' "${SAMPLES[@]}" | sed 's/,$//')"
    printf ' "map_perc_softclip":[%s],\n'      "$(printf '"0.75",%.0s' "${SAMPLES[@]}" | sed 's/,$//')"
    printf ' "is_umi":[%s],\n'                 "$(printf "\"$IS_UMI\",%.0s" "${SAMPLES[@]}" | sed 's/,$//')"
    printf ' "read_length":[%s],\n'            "$RL_JOINED"
    printf ' "cores":[%s],\n'                  "$(printf "\"$CORES\",%.0s" "${SAMPLES[@]}" | sed 's/,$//')"
    printf ' "ram":[%s]\n'                     "$(printf '"200",%.0s' "${SAMPLES[@]}" | sed 's/,$//')"
    printf '}\n'
} > "$CONFIG"
cat "$CONFIG"

cd "$HD_DIR"
mm_run hybridetector snakemake \
    --snakefile HybriDetector.smk \
    --configfile "$CONFIG" \
    --profile "$PROFILE" \
    --jobs "$JOBS" \
    ${EXTRA_SNAKE_ARGS:-}

echo
echo "Done. Hybrid tables:"
ls -la "$HD_DIR"/hyb_pairs/*unified_length_all_types_unique_high_confidence.tsv 2>/dev/null \
    || echo "  (none yet)"
