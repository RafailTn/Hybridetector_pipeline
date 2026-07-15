#!/usr/bin/env bash
# Does --genomeSAsparseD 2 change the alignments, or only the RAM?
#
#   bash chimeric_eclip/validate_sparse_index.sh [n_reads]
#
# hybridetector-local.patch builds the hg38 index with a sparse suffix array (every 2nd
# suffix) so it fits this 31 GB host. STAR documents that flag as a memory/speed knob,
# and a sparse SA is in principle a lossless representation of the same index -- but that
# is a claim about the algorithm, not about STAR 2.5.3's implementation, and miRBench
# built Manakov's index DENSE. If sparse perturbs mapping at all, then our dataset is
# processed differently from Manakov for a reason that is purely about our RAM. That is
# worth an experiment rather than an assumption.
#
# The experiment: chr21 is small enough to index densely without any RAM problem, so
# build it BOTH ways with the same STAR 2.5.3 and the same alignment parameters
# HybriDetector uses, push real reads from the sample through both, and diff the
# alignment records read by read.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

N_READS="${1:-400000}"
THREADS="${THREADS:-6}"          # modest: the real run may still be using the machine
WORK="${WORK:-$REPO_ROOT/results/sparse_index_check}"
FASTQ="${FASTQ:-$REPO_ROOT/data/raw/gse297116/HCT116_wholecell_rep1_IP.fastq.gz}"
GENOME="$HD_DIR/DBs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

STAR=""
for d in "$HD_DIR"/.snakemake/conda/*/; do
    [ -x "$d/bin/STAR" ] && { STAR="$d/bin/STAR"; SAMTOOLS="$d/bin/samtools"; break; }
done
[ -n "$STAR" ] || { echo "no STAR found in the rule envs (run setup.sh)" >&2; exit 1; }

# ── wait for the whole-genome index build to get out of the way ────────────
echo "== waiting for the main STAR genome index to finish"
prev=""
while :; do
    sz=$(stat -c %s "$HD_DIR/index/STAR/SAindex" 2>/dev/null || echo 0)
    # done when SAindex exists and has stopped growing
    [ "$sz" != "0" ] && [ "$sz" = "$prev" ] && break
    prev="$sz"
    sleep 30
done
echo "== genome index present ($(du -h "$HD_DIR/index/STAR/SAindex" | cut -f1) SAindex); starting chr21 check"

mkdir -p "$WORK"
cd "$WORK"

# ── chr21 reference ────────────────────────────────────────────────────────
if [ ! -s chr21.fa ]; then
    echo "== extracting chr21 from the Ensembl primary assembly"
    awk '/^>/{p = ($0 ~ /^>21 /)} p' "$GENOME" > chr21.fa
    grep -c "^>" chr21.fa
fi
echo "== chr21: $(grep -v '>' chr21.fa | wc -m) bases"

# STAR's own rule: min(14, log2(len)/2 - 1). chr21 is ~46 Mb -> 11.
NBASES=$(python3 -c "
import math
n = $(grep -v '>' chr21.fa | wc -m)
print(min(14, math.floor(math.log(float(n),2)/2-1)))")
echo "== --genomeSAindexNbases = $NBASES"

# ── build the two indices: dense (default) vs sparse (what we ship) ────────
for D in 1 2; do
    if [ ! -s "idx_sparseD$D/SAindex" ]; then
        echo "== building chr21 index with --genomeSAsparseD $D"
        mkdir -p "idx_sparseD$D"
        "$STAR" --runMode genomeGenerate --runThreadN "$THREADS" \
            --genomeDir "idx_sparseD$D" --genomeFastaFiles chr21.fa \
            --genomeSAindexNbases "$NBASES" --genomeSAsparseD "$D" \
            --outFileNamePrefix "idx_sparseD${D}_" > /dev/null
    fi
done
echo "== index sizes:"
du -sh idx_sparseD1 idx_sparseD2 | sed 's/^/   /'

# ── align the same reads through both ─────────────────────────────────────
if [ ! -s reads.fastq.gz ]; then
    echo "== taking $N_READS reads from $(basename "$FASTQ")"
    # `head` closes the pipe early, so zcat takes SIGPIPE; under `pipefail` that is a
    # fatal 141 even though the output is complete and correct. Take the hit locally.
    set +o pipefail
    zcat "$FASTQ" | head -n $((N_READS * 4)) | gzip > reads.fastq.gz
    set -o pipefail
fi

# Same alignment parameters as HybriDetector's alignment_single_genomic rule, so the
# comparison reflects how reads are actually mapped in the pipeline.
for D in 1 2; do
    if [ ! -s "aln_D$D.sam" ]; then
        echo "== aligning with sparseD=$D"
        "$STAR" --runMode alignReads --runThreadN "$THREADS" \
            --genomeDir "idx_sparseD$D" --readFilesIn reads.fastq.gz \
            --readFilesCommand zcat --outFileNamePrefix "aln_D${D}_" \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.1 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterMatchNmin 0 --outFilterScoreMinOverLread 0.85 \
            --outFilterMatchNminOverLread 0.85 \
            --outSAMunmapped Within --outSAMattributes All \
            --outMultimapperOrder Random --outSAMtype SAM > /dev/null
        mv "aln_D${D}_Aligned.out.sam" "aln_D$D.sam"
    fi
done

# ── compare ───────────────────────────────────────────────────────────────
echo
echo "======================= RESULT ======================="
for D in 1 2; do
    tot=$(grep -vc "^@" "aln_D$D.sam" || true)
    # unmapped records (kept by --outSAMunmapped Within) carry RNAME "*"
    mapped=$(awk '!/^@/ && $3 != "*"' "aln_D$D.sam" | wc -l)
    printf "sparseD=%s : %s records, %s mapped\n" "$D" "$tot" "$mapped"
done

# Per-read alignment identity: read name, flag, chrom, pos, CIGAR, NH tag.
# outMultimapperOrder Random means multimapper ORDER can differ run to run, so sort.
for D in 1 2; do
    awk '!/^@/ {nh="NA"; for(i=12;i<=NF;i++) if($i ~ /^NH:i:/) nh=$i;
                print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"nh}' "aln_D$D.sam" | sort > "cmp_D$D.tsv"
done

n1=$(wc -l < cmp_D1.tsv); n2=$(wc -l < cmp_D2.tsv)
diffs=$(comm -3 cmp_D1.tsv cmp_D2.tsv | wc -l)
echo
echo "alignment records: dense=$n1  sparse=$n2"
echo "records differing between dense and sparse: $diffs"
if [ "$diffs" -eq 0 ]; then
    echo
    echo "VERDICT: IDENTICAL. --genomeSAsparseD 2 is purely a RAM/speed tradeoff here;"
    echo "         mapping is unaffected, so the sparse index does not make this dataset"
    echo "         incomparable to Manakov."
else
    pct=$(python3 -c "print(f'{100*$diffs/max($n1,1):.4f}')")
    echo
    echo "VERDICT: DIFFERS in $diffs records ($pct% of dense). The sparse index is NOT"
    echo "         alignment-neutral for this STAR version. Examples:"
    comm -3 cmp_D1.tsv cmp_D2.tsv | head -10
fi
echo "====================================================="
