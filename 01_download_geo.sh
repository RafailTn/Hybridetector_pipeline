#!/usr/bin/env bash
# Stage 1 — fetch the raw reads for the samples you actually want.
#
# Two modes. First expand a GEO series into a sample sheet:
#
#   bash chimeric_eclip/01_download_geo.sh --gse GSE297116 --out data/raw/gse297116
#
# which writes <out>/samplesheet.tsv listing every GSM with its title and SRR runs, and
# stops. Delete the rows you don't want (see the note on inputs below), then:
#
#   bash chimeric_eclip/01_download_geo.sh --samplesheet <sheet> --out data/raw/gse297116
#
# which downloads each run and concatenates the runs of one sample into one FASTQ. A
# curated sheet is kept in chimeric_eclip/samplesheets/ for each experiment we process.
#
# WHICH SAMPLES: only the IP libraries carry chimeras. An eCLIP series also ships
# size-matched `(input)` controls — no AGO2 pulldown, so no miRNA:target ligations —
# and pushing those through HybriDetector just burns hours to find nothing.
#
# Reads are single-end here; where a run comes back paired, only R1 is kept, which is
# where the chimera lives (miRBench did the same for Manakov).
#
# Env vars: THREADS (default: CPU count), KEEP_SRA=1 to keep the .sra prefetch cache.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env eclip_dl

GSE=""; SHEET=""; OUT_DIR=""
while [ $# -gt 0 ]; do
    case "$1" in
        --gse)         GSE="$2"; shift 2 ;;
        --samplesheet) SHEET="$2"; shift 2 ;;
        --out)         OUT_DIR="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done
: "${OUT_DIR:?--out is required}"
THREADS="${THREADS:-$(n_cpus)}"
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

# ── mode 1: GSE -> sample sheet ────────────────────────────────────────────
if [ -n "$GSE" ]; then
    SHEET="$OUT_DIR/samplesheet.tsv"
    echo "== expanding $GSE"
    # The GEO SOFT record is the authoritative sample list; edirect's `gds` database
    # mangles this query, so go at the record directly.
    soft="$OUT_DIR/.$GSE.soft"
    curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$GSE&targ=gsm&form=text&view=brief" -o "$soft"

    printf '# sample_name\tsrr_runs\t# GSM / title\n' > "$SHEET"
    gsm=""; title=""
    while IFS= read -r line; do
        case "$line" in
            "^SAMPLE = "*)       gsm="${line#^SAMPLE = }" ;;
            "!Sample_title = "*) title="${line#!Sample_title = }" ;;
            *"!Sample_relation = SRA:"*)
                srx="${line##*term=}"
                runs=$(curl -s "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=$srx" \
                       | awk -F, 'NR>1 && $1 ~ /^[SED]RR/ {printf "%s,", $1}' | sed 's/,$//')
                # sample_name must be filesystem- and HybriDetector-safe
                safe=$(echo "$title" | tr -cs 'A-Za-z0-9' '_' | sed 's/^_//; s/_$//')
                printf '%s\t%s\t# %s / %s\n' "$safe" "$runs" "$gsm" "$title" >> "$SHEET"
                ;;
        esac
    done < "$soft"
    rm -f "$soft"

    echo
    cat "$SHEET"
    echo
    echo "Wrote $SHEET. Delete the rows you don't want (inputs carry no chimeras — keep the IPs), then rerun with:"
    echo "  bash $0 --samplesheet $SHEET --out $OUT_DIR"
    exit 0
fi

# ── mode 2: sample sheet -> fastq ──────────────────────────────────────────
: "${SHEET:?pass either --gse or --samplesheet}"
[ -s "$SHEET" ] || { echo "sample sheet not found: $SHEET" >&2; exit 1; }

while IFS=$'\t' read -r sample runs _; do
    case "$sample" in ''|'#'*) continue ;; esac
    final="$OUT_DIR/$sample.fastq.gz"
    if [ -s "$final" ]; then
        echo "== $sample already downloaded, skipping"
        continue
    fi
    echo "== $sample  <-  $runs"

    parts=()
    IFS=',' read -ra srrs <<< "$runs"
    for srr in "${srrs[@]}"; do
        if [ ! -s "$OUT_DIR/$srr.fastq" ]; then
            mm_run eclip_dl prefetch --max-size u -O "$OUT_DIR/sra" "$srr"
            mm_run eclip_dl fasterq-dump --split-files --threads "$THREADS" \
                -O "$OUT_DIR" "$OUT_DIR/sra/$srr/$srr.sra"
            # keep R1 only; single-end runs come out unsuffixed
            if [ -f "$OUT_DIR/${srr}_1.fastq" ]; then
                mv "$OUT_DIR/${srr}_1.fastq" "$OUT_DIR/$srr.fastq"
                rm -f "$OUT_DIR/${srr}_2.fastq"
            fi
            [ "${KEEP_SRA:-0}" = "1" ] || rm -rf "$OUT_DIR/sra/$srr"
        fi
        parts+=("$OUT_DIR/$srr.fastq")
    done

    # The runs of one SRX are the same library sequenced more than once (here: a deep
    # run plus a small top-off), so they concatenate into one sample.
    cat "${parts[@]}" | mm_run eclip_dl pigz -p "$THREADS" -c > "$final"
    rm -f "${parts[@]}"
done < "$SHEET"

echo
echo "Done. FASTQs in $OUT_DIR:"
ls -la "$OUT_DIR"/*.fastq.gz
echo
echo "Next: bash $(dirname "$0")/02_preprocess_fastq.sh $OUT_DIR <pp dir>"
