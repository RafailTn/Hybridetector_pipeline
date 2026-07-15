#!/usr/bin/env bash
# Shared paths for the chimeric-eCLIP -> v7-dataset pipeline. Sourced by every stage.
#
# The environments live under micromamba (NOT pixi): HybriDetector is a Snakemake
# workflow whose rules each build their own conda env at runtime, which needs a real
# conda/mamba solver. The pixi env in dependencies/ is for the CNN and is untouched.

MICROMAMBA="${MICROMAMBA:-$HOME/.local/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
EXTERNAL="${EXTERNAL:-$REPO_ROOT/external}"
HD_DIR="${HD_DIR:-$EXTERNAL/HybriDetector}"
MIRBENCH_DIR="${MIRBENCH_DIR:-$EXTERNAL/miRBench_paper}"

# Reference bigwigs for the conservation step (postprocess step 6).
PHYLOP_BW="${PHYLOP_BW:-$REPO_ROOT/data/hg38.phyloP100way.bw}"
PHASTCONS_BW="${PHASTCONS_BW:-$REPO_ROOT/data/hg38.phastCons100way.bw}"

# Run a command inside one of the four micromamba envs.
mm_run() {
    local env_name="$1"; shift
    "$MICROMAMBA" run -n "$env_name" "$@"
}

require_env() {
    if ! "$MICROMAMBA" env list | grep -qE "^\s+$1\s"; then
        echo "Error: micromamba env '$1' not found. Run: bash chimeric_eclip/setup.sh" >&2
        exit 1
    fi
}
