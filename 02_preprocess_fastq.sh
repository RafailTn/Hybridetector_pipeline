#!/usr/bin/env bash
# Stage 2 — trim raw reads into the form HybriDetector requires: adapters and barcodes
# precisely removed, UMI moved from the sequence into the read header.
#
#   bash chimeric_eclip/02_preprocess_fastq.sh data/raw/manakov data/pp/manakov
#
# This is a local, non-SLURM rewrite of
# miRBench_paper/code/preprocess_downloaded_eCLIP_for_HD/preprocess_eCLIP_for_HD.sh,
# which is a job array with hardcoded #SBATCH directives. Same three steps, same
# parameters, in the same order (YeoLab chim-eCLIP recipe):
#
#   1. umi_tools extract  — lift the 5' UMI into the read name
#   2. cutadapt           — trim the 3' adapter (tiled Illumina adapter suffixes)
#   3. cutadapt -u -10    — cut the 10 nt of 3' UMI remnant
#
# One deviation, forced: upstream passes `-f fastq` to cutadapt. That option was
# removed in cutadapt 2.0 and is a hard error on the 5.x in eclip_pp, so it is dropped
# (cutadapt infers the format anyway).
#
# PROTOCOL KNOBS — these are chim-eCLIP (Manakov 2022) defaults. A CLASH or CLEAR-CLIP
# library has a different barcode layout; override rather than assume:
#   UMI_PATTERN   5' UMI as a umi_tools barcode pattern ("" disables step 1)
#   UMI3_LEN      nt to cut from the 3' end after adapter trimming (0 disables step 3)
#   MIN_LEN       cutadapt -m
#
# !! CHECK THE READ STATE BEFORE RUNNING THIS !!
# The stage assumes RAW reads. Some series deposit reads that the submitter has ALREADY
# trimmed — GSE297116 is one: its UMIs were pruned and its adapters cut before upload.
# Running this stage on such reads is a silent catastrophe: it would strip 10 nt off the
# 5' end, which is exactly where the miRNA and its seed sit, and nothing would error.
# `bash chimeric_eclip/check_read_state.sh <fastq.gz>` tells you which case you are in.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env eclip_pp

IN_DIR="${1:?usage: 02_preprocess_fastq.sh <raw fastq dir> <output dir>}"
OUT_DIR="${2:?usage: 02_preprocess_fastq.sh <raw fastq dir> <output dir>}"

UMI_PATTERN="${UMI_PATTERN:-NNNNNNNNNN}"   # 10 nt 5' UMI
UMI3_LEN="${UMI3_LEN:-10}"                 # 10 nt 3' UMI remnant
MIN_LEN="${MIN_LEN:-18}"
THREADS="${THREADS:-$(nproc)}"

# Every 10-mer suffix of the Illumina 3' adapter, as upstream. Reads are short and the
# adapter is often only partially present, so cutadapt is given all the offsets.
ADAPTERS=(AGATCGGAAG GATCGGAAGA ATCGGAAGAG TCGGAAGAGC CGGAAGAGCA GGAAGAGCAC
          GAAGAGCACA AAGAGCACAC AGAGCACACG GAGCACACGT AGCACACGTC GCACACGTCT
          CACACGTCTG ACACGTCTGA CACGTCTGAA ACGTCTGAAC CGTCTGAACT GTCTGAACTC
          TCTGAACTCC CTGAACTCCA TGAACTCCAG GAACTCCAGT AACTCCAGTC ACTCCAGTCA)
ADAPTER_ARGS=()
for a in "${ADAPTERS[@]}"; do ADAPTER_ARGS+=(-a "$a"); done

mkdir -p "$OUT_DIR"
shopt -s nullglob
for fq in "$IN_DIR"/*.fastq.gz; do
    base="$(basename "$fq" .fastq.gz)"
    sample_dir="$OUT_DIR/$base"
    final="$sample_dir/$base.pp.fastq.gz"
    if [ -s "$final" ]; then
        echo "== $base already preprocessed, skipping"
        continue
    fi
    mkdir -p "$sample_dir/logs" "$sample_dir/temp"
    echo "== preprocessing $base"

    step1="$sample_dir/temp/$base.umi.fastq.gz"
    if [ -n "$UMI_PATTERN" ]; then
        mm_run eclip_pp umi_tools extract \
            --random-seed 1 \
            --stdin "$fq" \
            --bc-pattern "$UMI_PATTERN" \
            --log "$sample_dir/logs/$base.umi_tools.log" \
            --stdout "$step1"
    else
        step1="$fq"
    fi

    step2="$sample_dir/temp/$base.umi.adapter.fastq.gz"
    mm_run eclip_pp cutadapt \
        -O 1 --match-read-wildcards --times 3 -e 0.1 \
        --quality-cutoff 6 -m "$MIN_LEN" \
        "${ADAPTER_ARGS[@]}" \
        -j "$THREADS" \
        -o "$step2" "$step1" > "$sample_dir/logs/$base.3cutadapt.txt"

    if [ "$UMI3_LEN" -gt 0 ]; then
        mm_run eclip_pp cutadapt -u "-$UMI3_LEN" -j "$THREADS" \
            -o "$final" "$step2" > "$sample_dir/logs/$base.5cutadapt.txt"
    else
        cp "$step2" "$final"
    fi

    [ -s "$final" ] && rm -rf "$sample_dir/temp"
done

echo
echo "Done. Preprocessed FASTQs (feed these to stage 3):"
find "$OUT_DIR" -name "*.pp.fastq.gz" | sort
