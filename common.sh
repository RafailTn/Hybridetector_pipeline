#!/usr/bin/env bash
# Shared paths for the chimeric-eCLIP -> v7-dataset pipeline. Sourced by every stage.
#
# The environments live under micromamba (NOT pixi): HybriDetector is a Snakemake
# workflow whose rules each build their own conda env at runtime, which needs a real
# conda/mamba solver. The pixi env in dependencies/ is for the CNN and is untouched.

MICROMAMBA="${MICROMAMBA:-$HOME/.local/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

# Where external/ (the HD + miRBench clones) and data/ (bigwigs, dataset outputs) live.
# In the thesis these scripts sit in <thesis>/chimeric_eclip/ and deliberately SHARE the
# thesis-root external/ and data/ with the CNN side, so the root is the parent. But when
# the pipeline is cloned on its own — scripts flat at a repo's top level, no chimeric_eclip
# subfolder — the parent points ABOVE the clone and everything lands in the wrong place. So
# only climb to the parent when we really are inside a chimeric_eclip/ subfolder; otherwise
# the scripts' own directory is the root, and the standalone clone is self-contained.
# Override REPO_ROOT to force either.
# PIPELINE_DIR is where these scripts (and their siblings: the .py, patches/, envs/,
# slurm/) live. Sibling-file references must use it, never "$REPO_ROOT/chimeric_eclip",
# so they resolve in a flat clone too.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -z "${REPO_ROOT:-}" ]; then
    if [ "$(basename "$PIPELINE_DIR")" = "chimeric_eclip" ]; then
        REPO_ROOT="$(dirname "$PIPELINE_DIR")"
    else
        REPO_ROOT="$PIPELINE_DIR"
    fi
fi
export PIPELINE_DIR REPO_ROOT
EXTERNAL="${EXTERNAL:-$REPO_ROOT/external}"
HD_DIR="${HD_DIR:-$EXTERNAL/HybriDetector}"
MIRBENCH_DIR="${MIRBENCH_DIR:-$EXTERNAL/miRBench_paper}"

# Reference bigwigs for the conservation step (postprocess step 6).
PHYLOP_BW="${PHYLOP_BW:-$REPO_ROOT/data/hg38.phyloP100way.bw}"
PHASTCONS_BW="${PHASTCONS_BW:-$REPO_ROOT/data/hg38.phastCons100way.bw}"

# ── portable helpers (Linux + macOS) ───────────────────────────────────────
# The host-shell stages must run on macOS too, where nproc, /proc/meminfo, GNU
# `readlink -f`, `stat -c` and the bash-4 `mapfile` builtin are all absent (macOS ships
# bash 3.2). Route those through these shims; mapfile is replaced inline with read loops.

# Logical CPU count.
n_cpus() {
    if command -v nproc >/dev/null 2>&1; then nproc
    else sysctl -n hw.logicalcpu 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1
    fi
}

# Total physical RAM, whole GB.
mem_total_gb() {
    if [ -r /proc/meminfo ]; then
        awk '/^MemTotal:/ {printf "%d", $2/1024/1024}' /proc/meminfo
    else
        awk -v b="$(sysctl -n hw.memsize 2>/dev/null || echo 0)" \
            'BEGIN {printf "%d", b/1024/1024/1024}'
    fi
}

# Absolute, symlink-resolved path. GNU `readlink -f` when it works (Linux, macOS >=12.3),
# else a cd-based fallback that at least makes the path absolute.
abspath() {
    if readlink -f "$1" >/dev/null 2>&1; then
        readlink -f "$1"
    else
        printf '%s/%s\n' "$(cd "$(dirname "$1")" && pwd)" "$(basename "$1")"
    fi
}

# Size of a file in bytes (GNU stat vs BSD/macOS stat).
file_size() {
    stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null || echo 0
}

# Modification time of a file, epoch seconds (GNU stat vs BSD/macOS stat).
file_mtime() {
    stat -c %Y "$1" 2>/dev/null || stat -f %m "$1" 2>/dev/null || echo 0
}

# Run a command inside one of the four micromamba envs.
mm_run() {
    local env_name="$1"; shift
    "$MICROMAMBA" run -n "$env_name" "$@"
}

# Echo one read_length per input FASTQ, in argument order, one per line. HybriDetector derives
# the min miRNA-align fraction ceil(16/read_length*100)/100 *per sample*, so a batch that mixes
# libraries of different length must get a per-sample array, not one global max (too lenient for
# the shorter samples). 'auto' measures each file's own max insert length (seqkit stats col 8 =
# max_len); a numeric mode is echoed verbatim for every file (a manual override for all).
# Args: <seqkit_env> <mode: auto|N> <fastq>...
read_lengths_for() {
    local env_name="$1" mode="$2"; shift 2
    local fq
    for fq in "$@"; do
        if [ "$mode" = "auto" ]; then
            mm_run "$env_name" seqkit stats -T "$fq" | awk 'NR==2 {print $8+0}'
        else
            printf '%s\n' "$mode"
        fi
    done
}

require_env() {
    if ! "$MICROMAMBA" env list | grep -qE "^\s+$1\s"; then
        echo "Error: micromamba env '$1' not found. Run: bash $PIPELINE_DIR/setup.sh" >&2
        exit 1
    fi
}

# Refuse to reuse a genome STAR index whose --genomeSAsparseD differs from the one this run
# requests. STAR_gen_index takes the sparsity from the environment, but the rule declares
# only the genome FASTA as input, so Snakemake would silently REUSE an existing index built
# with the other sparsity while the run announces this one — e.g. a "DENSE" banner over a
# leftover sparse index that is NOT comparable to Manakov. STAR records what it actually used
# in genomeParameters.txt; compare and stop on mismatch. Fires only for a COMPLETE index
# (SAindex present) — an absent or partial one is rebuilt fresh at the requested sparsity.
# Arg: the requested --genomeSAsparseD (1 = dense, 2 = sparse).
check_star_index_sparsity() {
    local want="$1" dir="$HD_DIR/index/STAR"
    [ -s "$dir/SAindex" ] || return 0
    [ -f "$dir/genomeParameters.txt" ] || return 0
    local have
    have=$(awk '$1=="genomeSAsparseD"{print $2; exit}' "$dir/genomeParameters.txt")
    [ -n "$have" ] || return 0
    if [ "$have" != "$want" ]; then
        echo "Error: the genome STAR index in $dir was built with --genomeSAsparseD $have," >&2
        echo "  but this run requests $want (1=dense, 2=sparse). Snakemake does not track that" >&2
        echo "  flag, so it would REUSE the on-disk ($have) index while announcing $want, silently" >&2
        echo "  producing a dataset that does not match its own log." >&2
        echo "  -> rebuild at the requested sparsity:  rm -rf $dir" >&2
        echo "  -> or accept the existing index:       STAR_SA_SPARSE_D=$have <this command>" >&2
        exit 1
    fi
}
