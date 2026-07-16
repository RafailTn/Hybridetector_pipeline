#!/usr/bin/env bash
# One-time setup: clone the two upstream repos, apply the upstream bug fixes,
# and build the four micromamba envs. Idempotent — safe to re-run.
#
#   bash chimeric_eclip/setup.sh
#
# Afterwards, prebuild HybriDetector's 15 per-rule conda envs (optional, ~4 GB, but it
# gets all the solving out of the way before you start a multi-hour run):
#
#   bash chimeric_eclip/setup.sh --prebuild-hd-envs
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
ENVS_DIR="$PIPELINE_DIR/envs"

# ── micromamba ─────────────────────────────────────────────────────────────
if [ ! -x "$MICROMAMBA" ]; then
    echo "== installing micromamba to $MICROMAMBA"
    mkdir -p "$(dirname "$MICROMAMBA")"
    # The tarball holds bin/micromamba. Unpack to a scratch dir and move the binary to
    # wherever $MICROMAMBA actually points, rather than trying to line the archive's
    # layout up with the destination via --strip-components — that only works if
    # $MICROMAMBA happens to end in .../bin/micromamba, and silently drops the binary
    # one directory too high when it doesn't.
    _mm_tmp="$(mktemp -d "${TMPDIR:-/tmp}/mm.XXXXXX")"
    trap 'rm -rf "$_mm_tmp"' EXIT
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C "$_mm_tmp" bin/micromamba
    mv "$_mm_tmp/bin/micromamba" "$MICROMAMBA"
    chmod +x "$MICROMAMBA"
    rm -rf "$_mm_tmp"; trap - EXIT
fi

# ── upstream repos ─────────────────────────────────────────────────────────
mkdir -p "$EXTERNAL"
if [ ! -d "$HD_DIR" ]; then
    echo "== cloning HybriDetector (fix_clustering branch)"
    # The paper used this branch, not main: it carries the clustering fix and the
    # large-file handling needed for the 19 concatenated Manakov samples.
    git clone --branch fix_clustering https://github.com/ML-Bioinfo-CEITEC/HybriDetector.git "$HD_DIR"
fi
if [ ! -d "$MIRBENCH_DIR" ]; then
    echo "== cloning miRBench_paper"
    git clone https://github.com/BioGeMT/miRBench_paper.git "$MIRBENCH_DIR"
fi

# ── upstream bug fixes (portable: needed on laptop AND cluster) ────────────
# python<3.12 pin (snakemake 7 imports distutils), jq for the bsgenome post-link, the
# gzipped-genome length bug in STAR_gen_index, and FTP->HTTPS for Ensembl. It also makes
# the memory budgets and STAR sparsity env-switchable, with UPSTREAM values as defaults.
# See chimeric_eclip/patches/hybridetector-fixes.patch and the README.
if git -C "$HD_DIR" diff --quiet; then
    echo "== applying hybridetector-fixes.patch"
    git -C "$HD_DIR" apply "$PIPELINE_DIR/patches/hybridetector-fixes.patch"
else
    # A dirty tree means the patch is already applied — but it also means a *newer* patch
    # revision won't land here. To pick up patch changes, reset the clone first:
    #   git -C external/HybriDetector checkout .   (or delete it and let setup.sh re-clone)
    echo "== HybriDetector already patched (working tree dirty), skipping"
    echo "   (to apply an updated patch: git -C \"$HD_DIR\" checkout . && rerun setup.sh)"
fi
mkdir -p "$HD_DIR/data"

# Makes Nunique fall back to NA when the source library had no recoverable UMIs, so a
# missing read-support count cannot be mistaken for a real one. No-op on the canonical
# datasets, which all carry UMIs. See the "Datasets without UMIs" section of the README.
if git -C "$MIRBENCH_DIR" diff --quiet; then
    echo "== applying mirbench-nunique-na.patch"
    git -C "$MIRBENCH_DIR" apply "$PIPELINE_DIR/patches/mirbench-nunique-na.patch"
else
    echo "== miRBench_paper already patched (working tree dirty), skipping"
fi

# ── environments ───────────────────────────────────────────────────────────
for env_name in hybridetector eclip_pp eclip_dl mirbench_pp; do
    if "$MICROMAMBA" env list | grep -qE "^\s+$env_name\s"; then
        echo "== env '$env_name' exists, skipping"
    else
        echo "== creating env '$env_name'"
        "$MICROMAMBA" create -y -f "$ENVS_DIR/$env_name.yml"
    fi
done

# jq is not a declared dependency of bioconductor-data-packages, but the post-link
# script of bioconductor-bsgenome.hsapiens.ucsc.hg38 shells out to it via yq. Without
# it on PATH, the filter_and_collapse rule's env fails to build. Putting it in the
# driver env is enough: env creation inherits this PATH.
"$MICROMAMBA" install -y -n hybridetector -c conda-forge jq >/dev/null

if [ "${1:-}" = "--prebuild-hd-envs" ]; then
    echo "== prebuilding HybriDetector per-rule conda envs (~5 GB, tens of minutes)"
    # Snakemake needs a resolvable DAG to know which envs to build, so give it dummy
    # samples. Nothing is aligned; --conda-create-envs-only stops after env creation.
    # Use *two* samples: the merge_replicates rule only enters the DAG when there is
    # more than one, and it is one of the slow bsgenome envs — you do not want to
    # discover that hours into a real multi-sample run.
    printf '@r1\nACGT\n+\nIIII\n' | gzip > "$HD_DIR/data/_ENVBUILD.fastq.gz"
    printf '@r1\nACGT\n+\nIIII\n' | gzip > "$HD_DIR/data/_ENVBUILD2.fastq.gz"
    cat > "$HD_DIR/config_envbuild.json" <<'JSON'
{"Sample":["_ENVBUILD","_ENVBUILD2"],
 "map_perc_single_genomic":["0.85","0.85"], "map_perc_softclip":["0.75","0.75"],
 "is_umi":["FALSE","FALSE"], "read_length":["75","75"],
 "cores":["20","20"], "ram":["28","28"]}
JSON
    ( cd "$HD_DIR" && mm_run hybridetector snakemake \
        --snakefile HybriDetector.smk --configfile config_envbuild.json \
        --use-conda --conda-frontend mamba --conda-create-envs-only \
        --res mem=28 -j "$(n_cpus)" )
    rm -f "$HD_DIR/data/_ENVBUILD.fastq.gz" "$HD_DIR/data/_ENVBUILD2.fastq.gz"
fi

echo
echo "Setup complete. Environments:"
"$MICROMAMBA" env list
