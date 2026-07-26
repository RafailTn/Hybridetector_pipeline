"""Generate miRBench-style negatives for a *subset* of a series' positives.

Why this exists: `positive_overlap.py --novel-out` gives the positives of a new series that
are not already in the miRBench reference datasets. Those files are positives-only, so a
model scored on them yields recall and nothing else — no AUROC, no threshold. This pairs each
novel positive with one negative, using upstream's procedure
(external/miRBench_paper/code/make_neg_sets/make_neg_sets.py), so the novel set becomes a
balanced, evaluable dataset.

It is NOT a rewrite of upstream: the sampling, the per-family 1:1 ratio, and the
sha256(noncodingRNA_name)-derived seed are reproduced exactly, so a negative drawn here for a
given family block matches what upstream would draw given the same pool. Three deliberate
deviations, each because the input is a subset of a finished pipeline rather than its step-2
output:

  1. BLOCKS AND POOL ARE DIFFERENT FILES. Upstream draws negatives for every positive in its
     input and pools decoys from that same input. Here --novel supplies the positives that
     need negatives, while --pool supplies the decoy sites (the series' step6 train+test,
     label==1 — i.e. ALL of its positives). Pooling from the novel subset alone would shrink
     the eligible-cluster set for no reason and would exclude positives from large families
     via upstream's "not enough negative examples" path.

  2. THE CLUSTER-EXCLUSION SET IS COMPUTED ON THE POOL, NOT ON THE BLOCK. This is the
     load-bearing one. Upstream forbids family F from drawing a decoy in any cluster F is
     seen binding. If that exclusion were computed from the novel block, a cluster that F
     demonstrably binds elsewhere in the series (but not in the novel subset) would be
     eligible as an F negative — manufacturing a false negative out of a known interaction.
     Computing it on the pool is therefore strictly safer than upstream, never laxer.

  3. CONSERVATION RIDES ALONG. Upstream makes negatives at step 3 and adds gene_phyloP /
     gene_phastCons from the bigwigs at step 6; the pool here is already post-step6. A
     negative copies its site wholesale from a pool row, and conservation is a pure function
     of (chr, start, end) — so carrying the two columns over reproduces what the bigwig
     readout would return, with no bigwig needed. Every site-derived column is copied as a
     unit for the same reason; only the miRNA columns come from the block.

The v7 `test` column is already dropped from these inputs (pipeline step 5), so unlike
upstream this writes no train/test split: the novel set is an evaluation set in full.

Usage:
    make_matched_negatives.py --novel data/gse297116_positives_novel_v7.tsv \
        --pool data/step6_add_conservation/gse297116...train....tsv \
               data/step6_add_conservation/gse297116...test....tsv \
        --output data/gse297116_novel_labelled_v7.tsv
"""
import argparse
import hashlib

import pandas as pd

# Columns describing the SITE — copied as a unit from the sampled pool row, so a negative is a
# real binding site with every one of its own coordinate-derived annotations intact.
SITE_COLS = [
    "gene", "feature", "chr", "start", "end", "strand", "gene_cluster_ID",
    "dominant_region", "regions_present",
    "read_start_in_sel_tx_1based", "read_end_in_sel_tx_1based",
    "gene_phyloP", "gene_phastCons",
]
# Columns describing the miRNA — taken from the positive block, i.e. the re-pairing itself.
MIRNA_COLS = ["noncodingRNA", "noncodingRNA_name", "noncodingRNA_fam"]


def _seed_for(block):
    """Upstream's reproducible per-block seed: sha256 of the miRNA name, mod 2^32-1."""
    name = block["noncodingRNA_name"].iloc[0]
    return int(hashlib.sha256(name.encode()).hexdigest(), 16) % 4294967295


def _sub_blocks(block):
    """Upstream's 'unknown' handling: an unknown family is not a family, so split it into
    one block per distinct miRNA sequence rather than treating the whole bucket as related."""
    if block["noncodingRNA_fam"].iloc[0] == "unknown":
        for mirna in block["noncodingRNA"].unique().tolist():
            yield block[block["noncodingRNA"] == mirna]
    else:
        yield block


def _excluded_clusters(sub_block, pool):
    """Clusters this family (or, when unknown, this miRNA) is seen binding anywhere in the
    POOL. See deviation 2 in the module docstring: computing this on the pool rather than on
    the novel block is what stops a known interaction being handed back as a negative."""
    fam = sub_block["noncodingRNA_fam"].iloc[0]
    if fam == "unknown":
        seen = pool[pool["noncodingRNA"] == sub_block["noncodingRNA"].iloc[0]]
    else:
        seen = pool[pool["noncodingRNA_fam"] == fam]
    return set(seen["gene_cluster_ID"].unique().tolist())


def _negatives_for(sub_block, pool, columns):
    seed = _seed_for(sub_block)
    blocked = _excluded_clusters(sub_block, pool)

    eligible = pool[~pool["gene_cluster_ID"].isin(blocked)]
    # One decoy per cluster: without this a cluster's many observed sites would let the same
    # target region be drawn repeatedly for one family.
    eligible = eligible.sample(frac=1, random_state=seed).drop_duplicates(
        subset=["gene_cluster_ID"], keep="first"
    )

    label = sub_block["noncodingRNA_fam"].iloc[0]
    if label == "unknown":
        label = sub_block["noncodingRNA"].iloc[0]

    if len(eligible) == 0:
        print(f"  WARNING: block {label} — {len(sub_block)} positives dropped "
              f"(no eligible negative clusters)", flush=True)
        return None, None

    num_neg = len(sub_block)
    if num_neg > len(eligible):
        print(f"  WARNING: block {label} — only {len(eligible)} eligible clusters for "
              f"{num_neg} positives; {num_neg - len(eligible)} positives dropped", flush=True)
        sub_block = sub_block.sample(n=len(eligible), random_state=seed)
        num_neg = len(eligible)

    negatives = eligible.sample(n=num_neg, random_state=seed)[SITE_COLS].copy()
    for col in MIRNA_COLS:
        negatives[col] = sub_block[col].values
    negatives["label"] = 0
    # Upstream's rule, kept verbatim: Nunique is read support, which a negative has none of.
    # But if the positives carry no read support either (a library whose UMIs were discarded,
    # e.g. GSE297116), writing 0 here would make the column NA for one class and 0 for the
    # other — a perfect label separator disguised as data. Keep it NA in that case.
    negatives["Nunique"] = 0 if sub_block["Nunique"].notna().any() else pd.NA
    return sub_block, negatives[columns]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--novel", required=True, help="positives needing negatives (v7)")
    ap.add_argument("--pool", required=True, nargs="+",
                    help="v7 files supplying the decoy pool; label==1 rows are used")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    novel = pd.read_csv(args.novel, sep="\t", dtype={"chr": str})
    pool = pd.concat([pd.read_csv(p, sep="\t", dtype={"chr": str}) for p in args.pool],
                     ignore_index=True)
    pool = pool[pool["label"] == 1]

    if list(novel.columns) != list(pool.columns):
        raise SystemExit(f"schema mismatch:\n  novel: {list(novel.columns)}\n  pool:  {list(pool.columns)}")
    if (novel["label"] != 1).any():
        raise SystemExit(f"--novel must be positives only; found {(novel['label'] != 1).sum()} non-positive rows")

    columns = list(novel.columns)
    print(f"{len(novel)} novel positives; pool = {len(pool)} positives "
          f"over {pool['gene_cluster_ID'].nunique()} clusters", flush=True)

    kept, out = 0, []
    for _, block in novel.groupby("noncodingRNA_fam", sort=False, dropna=False):
        for sub_block in _sub_blocks(block):
            pos, neg = _negatives_for(sub_block, pool, columns)
            if neg is None:
                continue
            out.append(pos[columns])
            out.append(neg)
            kept += len(pos)

    result = pd.concat(out, ignore_index=True)
    result.to_csv(args.output, sep="\t", index=False)
    dropped = len(novel) - kept
    print(f"wrote {len(result)} rows -> {args.output} "
          f"({kept} positives + {kept} negatives; {dropped} positives dropped)")


if __name__ == "__main__":
    main()
