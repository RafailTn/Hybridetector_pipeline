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
# PAIRED RUNS: HybriDetector is single-end (one data/<sample>.fastq.gz per sample), so only
# one mate can go downstream. `MATE` (default 1) picks which. Mate 1 is right for chim-eCLIP
# — Manakov's UMI+miRNA sit at the 5' end of R1, and miRBench used R1 — but that is a fact
# about THAT protocol, not about chimeric protocols in general. A CLASH/CLEAR-CLIP variant
# can carry the miRNA (or part of the duplex) in R2, and then mate 1 is the wrong choice and
# real signal is lost. So check your protocol's read layout before trusting the default.
#
# The unused mate is NOT deleted: it is parked, gzipped, under <out>/unused_mates/ so the
# choice stays reversible (DROP_R2=1 to discard it instead and save the disk). It lives in a
# subdirectory on purpose — stage 3 globs <dir>/*.fastq.gz at maxdepth 1 to find samples, so
# a parked mate sitting beside them would be picked up as a bogus extra sample.
#
# Env vars: THREADS (default: CPU count), KEEP_SRA=1 to keep the .sra prefetch cache,
#           MATE=1|2 (which mate carries the chimera), DROP_R2=1 to discard the unused mate.
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
MATE="${MATE:-1}"            # which mate carries the chimera (1 = chim-eCLIP / Manakov)
DROP_R2="${DROP_R2:-0}"      # 1 = discard the unused mate instead of parking it
case "$MATE" in 1|2) ;; *) echo "Error: MATE must be 1 or 2 (got '$MATE')" >&2; exit 1 ;; esac
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
            # --split-files gives <srr>_1/_2 when paired, plain <srr> when single-end (though
            # some single-end runs still come back as _1 only — hence the three-way branch).
            if [ -f "$OUT_DIR/${srr}_1.fastq" ] && [ -f "$OUT_DIR/${srr}_2.fastq" ]; then
                other=$(( MATE == 1 ? 2 : 1 ))
                echo "   !! $srr is PAIRED-END. Keeping mate $MATE as the chimeric read;"
                echo "      mate $other cannot go downstream (HybriDetector is single-end)."
                echo "      Mate 1 is correct for chim-eCLIP (UMI+miRNA at the 5' end of R1)."
                echo "      If your protocol puts the chimera in the other mate, set MATE=$other."
                if [ "$DROP_R2" = "1" ]; then
                    echo "      DROP_R2=1 -> discarding mate $other"
                    rm -f "$OUT_DIR/${srr}_${other}.fastq"
                else
                    mkdir -p "$OUT_DIR/unused_mates"
                    mm_run eclip_dl pigz -p "$THREADS" -c "$OUT_DIR/${srr}_${other}.fastq" \
                        > "$OUT_DIR/unused_mates/${srr}_${other}.fastq.gz"
                    rm -f "$OUT_DIR/${srr}_${other}.fastq"
                    echo "      mate $other kept at unused_mates/${srr}_${other}.fastq.gz (DROP_R2=1 to discard)"
                fi
                mv "$OUT_DIR/${srr}_${MATE}.fastq" "$OUT_DIR/$srr.fastq"
            elif [ -f "$OUT_DIR/${srr}_1.fastq" ]; then
                mv "$OUT_DIR/${srr}_1.fastq" "$OUT_DIR/$srr.fastq"   # single-end, _1-suffixed
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
