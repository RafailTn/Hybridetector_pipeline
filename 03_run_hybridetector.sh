#!/usr/bin/env bash
# Stage 3 — call chimeras. Runs the HybriDetector Snakemake workflow over the
# preprocessed FASTQs and produces, per sample, a table of miRNA→target hybrids.
#
#   bash chimeric_eclip/03_run_hybridetector.sh data/pp/manakov
#
# The output that matters downstream is
#   external/HybriDetector/hyb_pairs/<sample>.unified_length_all_types_unique_high_confidence.tsv
# (raw column names: seq.g, noncodingRNA_seq, chr.g, Nunique, ...). Note the *other*
# output, *_high_confidence_finalout.tsv, holds the same rows under the publication
# column names from the README — stage 4 needs the raw one, not that.
#
# When several samples are given, HybriDetector merges them and writes Merged.*.
#
# MEMORY IS THE BINDING CONSTRAINT, and it decides something scientific, not just how
# fast this goes. A dense hg38 STAR index needs ~32 GB resident to build and to align
# against. Below that you must fall back to a sparse suffix array (--genomeSAsparseD 2),
# which halves the RAM — but it is NOT alignment-neutral: on a chr21 control (164k real
# chr21 reads), dense-vs-dense differs for 1224 reads (0.75%, the RNG noise floor from
# --outMultimapperOrder Random) while dense-vs-sparse differs for 6931 (4.22%, 5.7x the
# floor). So sparse is a hardware concession that perturbs the mapping, not a free win, and
# a sparse-built dataset is not strictly comparable to miRBench's Manakov v7 (built dense).
#
# This script therefore SIZES ITSELF to the host: >=40 GB of RAM gets the dense index and
# upstream-ish budgets; anything less gets sparse, because dense would simply be OOM-killed
# (it was, here, at 27.9 GB). Snakemake's `--res mem` serialises the hungry rules so two
# can never run at once. Everything below is overridable, but if you force dense on a small
# box you will get an OOM kill, not a slowdown.
#
# Env vars: CORES, MEM_GB, READ_LENGTH, IS_UMI, STAR_SA_SPARSE_D, HD_MEM_*, EXTRA_SNAKE_ARGS.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env hybridetector

# Takes stage 2's output dir, OR — when the submitter already trimmed the reads and stage
# 2 was rightly skipped (see check_read_state.sh) — stage 1's raw dir directly.
PP_DIR="${1:?usage: 03_run_hybridetector.sh <fastq dir: stage 2 output, or stage 1 if already trimmed>}"

CORES="${CORES:-$(n_cpus)}"
READ_LENGTH="${READ_LENGTH:-auto}"       # 'auto' = measure it; see below
IS_UMI="${IS_UMI:-TRUE}"                 # TRUE only if the UMI is in the read header

# ── SIZE THE RUN TO THE HOST ───────────────────────────────────────────────
# Default to 90% of physical RAM: Snakemake's --res mem is a *budget it schedules
# against*, so claiming all of it would let the workflow commit the whole machine and
# leave nothing for the OS and the conda/python overhead of the rules themselves.
TOTAL_GB=$(mem_total_gb)
MEM_GB="${MEM_GB:-$(( TOTAL_GB * 9 / 10 ))}"
[ "$MEM_GB" -ge 8 ] || { echo "Error: MEM_GB=$MEM_GB is too small to run STAR at all." >&2; exit 1; }

# The dense/sparse fork. 40 GB is the threshold: a dense hg38 index peaked at 27.9 GB
# resident during `packing SA` here and was OOM-killed on a 31 GB box, and alignment
# then has to hold the whole index too, so ~32 GB is the floor and 40 leaves headroom.
if [ "$MEM_GB" -ge 40 ]; then
    DEFAULT_SPARSE=1                     # dense = STAR's default = what miRBench used
else
    DEFAULT_SPARSE=2                     # sparse: halves the RAM, perturbs the mapping
fi
export STAR_SA_SPARSE_D="${STAR_SA_SPARSE_D:-$DEFAULT_SPARSE}"

# HD_MEM_INDEX doubles as STAR's --limitGenomeGenerateRAM, so give the index rule the
# whole budget — `--res mem` guarantees it runs alone anyway. The other two keep upstream's
# values unless the host is smaller than that, hence the min().
min() { [ "$1" -lt "$2" ] && echo "$1" || echo "$2"; }
export HD_MEM_INDEX="${HD_MEM_INDEX:-$MEM_GB}"
export HD_MEM_ALIGN="${HD_MEM_ALIGN:-$(min "$MEM_GB" 50)}"
export HD_MEM_SMALL="${HD_MEM_SMALL:-$(min "$MEM_GB" 34)}"

if [ "$STAR_SA_SPARSE_D" = "1" ]; then
    echo "== ${TOTAL_GB} GB host -> budget ${MEM_GB} GB, DENSE STAR index (comparable to miRBench's Manakov v7)"
else
    echo "== ${TOTAL_GB} GB host -> budget ${MEM_GB} GB, SPARSE STAR index (--genomeSAsparseD 2)"
    echo "   Sparse is a memory concession and it CHANGES THE MAPPING (~5.7x the dense-vs-dense"
    echo "   noise floor on a chr21 control). Fine for a pilot; for a dataset you intend to"
    echo "   compare against Manakov v7, run this on a >=40 GB machine."
fi
# The announced mode is a promise about the index on disk: a leftover index of the other
# sparsity would be reused silently (Snakemake does not track the flag). Enforce the promise.
check_star_index_sparsity "$STAR_SA_SPARSE_D"
# Snakemake refuses to schedule a job whose resources exceed the global budget, so this
# would otherwise surface as a confusing "job needs mem=200, only NN available" abort.
for v in HD_MEM_INDEX HD_MEM_ALIGN HD_MEM_SMALL; do
    [ "${!v}" -le "$MEM_GB" ] || { echo "Error: $v=${!v} exceeds MEM_GB=$MEM_GB; nothing could ever be scheduled." >&2; exit 1; }
done

# HybriDetector reads its inputs from HD_DIR/data/<sample>.fastq.gz and is hardwired to
# run inside its own directory, so link the reads into place. Prefer stage 2's *.pp.
# files; fall back to plain *.fastq.gz for the already-trimmed case.
mkdir -p "$HD_DIR/data"
FQS=()
while IFS= read -r fq; do FQS+=("$fq"); done < <(find "$PP_DIR" -name "*.pp.fastq.gz" | sort)
SUFFIX=".pp.fastq.gz"
if [ "${#FQS[@]}" -eq 0 ]; then
    while IFS= read -r fq; do FQS+=("$fq"); done < <(find "$PP_DIR" -maxdepth 1 -name "*.fastq.gz" | sort)
    SUFFIX=".fastq.gz"
    [ "${#FQS[@]}" -gt 0 ] && echo "== no stage-2 output here; using raw *.fastq.gz (already-trimmed reads)"
fi

if [ "${#FQS[@]}" -eq 0 ]; then
    echo "Error: no *.pp.fastq.gz or *.fastq.gz found under $PP_DIR" >&2
    exit 1
fi

SAMPLES=()
for fq in "${FQS[@]}"; do
    sample="$(basename "$fq" "$SUFFIX")"
    ln -sfn "$(abspath "$fq")" "$HD_DIR/data/$sample.fastq.gz"
    SAMPLES+=("$sample")
done
echo "== ${#SAMPLES[@]} sample(s): ${SAMPLES[*]}"

# read_length is not cosmetic: alignment_small derives the minimum fraction of a read
# that must align to the miRNA from it (ceil(16/read_length*100)/100, i.e. "at least
# 16 nt of noncoding RNA"). Set it too high and shorter, flimsier miRNA alignments pass;
# too low and real chimeras are rejected. It must be the max length of the reads
# HybriDetector actually sees — i.e. *after* stage 2 stripped the UMIs and adapters —
# which is not the read length GEO advertises. So measure it rather than assume it, and
# do so PER SAMPLE: a batch can mix libraries of different length (pre-trimmed reads are
# variable-length), and one global max would over-set the threshold for the shorter ones.
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

# Build the config HybriDetector.py would have written, but without going through it:
# it hardcodes `--res mem=<ram>`, and we need mem and the STAR RAM cap decoupled.
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
    printf ' "ram":[%s]\n'                     "$(printf "\"$MEM_GB\",%.0s" "${SAMPLES[@]}" | sed 's/,$//')"
    printf '}\n'
} > "$CONFIG"

echo "== config: $CONFIG"
cat "$CONFIG"

# First run also downloads the Ensembl release-90 primary assembly (~900 MB) and builds
# ~10 STAR indices; budget several hours and ~40 GB of disk for that alone.
cd "$HD_DIR"
mm_run hybridetector snakemake \
    --snakefile HybriDetector.smk \
    --configfile "$CONFIG" \
    --use-conda --conda-frontend mamba \
    --res mem="$MEM_GB" \
    -j "$CORES" \
    -p --rerun-incomplete \
    ${EXTRA_SNAKE_ARGS:-}

echo
echo "Done. Hybrid tables:"
ls -la "$HD_DIR"/hyb_pairs/*unified_length_all_types_unique_high_confidence.tsv 2>/dev/null \
    || echo "  (none — dry run, or the workflow did not reach unify_length)"
