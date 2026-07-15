#!/usr/bin/env bash
# Are these reads raw, or did the submitter already trim them? Run BEFORE stage 2.
#
#   bash chimeric_eclip/check_read_state.sh data/raw/gse297116/HCT116_wholecell_rep1_IP.fastq.gz
#
# See check_read_state.py for what it measures and why the GEO metadata is not enough.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
require_env eclip_pp

FASTQ="${1:?usage: check_read_state.sh <fastq.gz> [n_reads]}"
N="${2:-400000}"
MIRNA_FA="${MIRNA_FA:-$HD_DIR/DBs/mirna/hsa_mature.fa}"

[ -s "$MIRNA_FA" ] || { echo "mature miRNA fasta not found: $MIRNA_FA (run setup.sh)" >&2; exit 1; }

mm_run eclip_pp python "$REPO_ROOT/chimeric_eclip/check_read_state.py" "$FASTQ" "$MIRNA_FA" "$N"
