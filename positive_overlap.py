"""How many of a new series' positive interactions are already in the existing datasets?

Every chimeric-eCLIP series this pipeline builds (GSE297116/HCT116, GSE304955/SAEC, …)
produces a set of positive (miRNA, target-site) pairs. Before treating those as *new*
training signal, you want to know how many are simply re-observations of interactions the
miRBench datasets already contain — Manakov being ~1.4M positives, so most of the prior.

This is the script behind the numbers in each series' write-up. It was run inline three
times before it earned a file; the logic is unchanged from those runs.

THREE definitions of "already seen", because they answer different questions and the right
one to quote is the middle:

  1. exact coords + miRNA    same chr/start/end/strand AND same miRNA sequence. Almost
                             always an undercount: HybriDetector anchors its 50 nt window
                             on the chimeric read, so two experiments rarely produce
                             byte-identical intervals even at the same site.
  2. site overlap + miRNA    the interval overlaps a reference positive's interval (same
                             chr/strand) AND the miRNA sequence matches. THE ONE TO QUOTE —
                             it is robust to the window-placement jitter that sinks (1).
  3. site overlap, any miRNA the interval overlaps a reference positive at all, regardless
                             of which miRNA. The gap between (2) and (3) is the fraction of
                             sites that are the SAME AGO2-bound locus but a DIFFERENT miRNA
                             — i.e. novelty that lives in the assignment, not the target.

Overlap uses a half-open interval test (a.start < b.end AND a.end > b.start), grouped by
(chr, strand) and bisected on the sorted reference starts so it stays fast against Manakov.

  --novel-out is written under definition (2): a query positive is novel unless some
  reference positive overlaps its interval with the same miRNA. That is the file you feed
  onward; the merged positives (all of them) go to --merged-out.

Usage:
    positive_overlap.py --input SERIES_v7.tsv [more_v7.tsv ...] \
        [--ref REF_v7.tsv ...] \
        [--merged-out positives_merged.tsv] [--novel-out positives_novel.tsv]

--input takes any number of v7 TSVs (e.g. the two step6 train/test conservation files);
they are concatenated and filtered to label==1. --ref defaults to the three existing
miRBench datasets (Manakov train+test+leftout, Hejret train+test, Klimentova test); pass
your own to compare against a different baseline. Reference files are read as label==1 only.
"""
import argparse
import glob
import sys

import numpy as np
import polars as pl

# chr is a string (X, Y, MT); start/end are read as strings then cast, because at least one
# HybriDetector-derived row writes a coordinate in scientific notation (1.9e+07) that a
# strict i64 parse chokes on. Cast via Float64 to absorb that, then to Int64.
_KEYCOLS = ["gene", "noncodingRNA", "noncodingRNA_fam", "chr", "start", "end", "strand", "label"]
_EXACT = ["chr", "start", "end", "strand", "noncodingRNA"]

DEFAULT_REFS = [
    "data/AGO2_eCLIP_Manakov2022_train_v7.tsv",
    "data/AGO2_eCLIP_Manakov2022_test_v7.tsv",
    "data/AGO2_eCLIP_Manakov2022_leftout_v7.tsv",
    "data/AGO2_CLASH_Hejret2023_train_v7.tsv",
    "data/AGO2_CLASH_Hejret2023_test_v7.tsv",
    "data/AGO2_eCLIP_Klimentova2022_test_v7.tsv",
]


def _read_positives(path):
    """v7 TSV -> label==1 rows, with start/end as Int64. Robust to sci-notation coords."""
    df = pl.read_csv(path, separator="\t", schema_overrides={"chr": pl.String},
                     infer_schema=False)
    return df.filter(pl.col("label") == "1").with_columns(
        pl.col("start").cast(pl.Float64).cast(pl.Int64),
        pl.col("end").cast(pl.Float64).cast(pl.Int64),
    )


def _overlap_masks(query, ref):
    """(overlap_any, overlap_and_same_mirna) boolean masks over query rows.

    Half-open interval overlap grouped by (chr, strand); reference starts are sorted and
    bisected, so a query only scans references whose start is within the widest reference
    interval of its own [start, end) window.
    """
    idx = {}
    for (c, s), sub in ref.group_by(["chr", "strand"]):
        sub = sub.sort("start")
        st, en, mi = sub["start"].to_numpy(), sub["end"].to_numpy(), sub["noncodingRNA"].to_numpy()
        idx[(c, s)] = (st, en, mi, int((en - st).max()) if len(en) else 0)
    any_ov = np.zeros(len(query), bool)
    mir_ov = np.zeros(len(query), bool)
    for i, r in enumerate(query.iter_rows(named=True)):
        k = (r["chr"], r["strand"])
        if k not in idx:
            continue
        st, en, mi, maxlen = idx[k]
        lo = np.searchsorted(st, r["start"] - maxlen, "left")
        hi = np.searchsorted(st, r["end"], "left")
        if lo >= hi:
            continue
        h = en[lo:hi] > r["start"]
        if h.any():
            any_ov[i] = True
            mir_ov[i] = bool((mi[lo:hi][h] == r["noncodingRNA"]).any())
    return any_ov, mir_ov


def _exact_mask(query, ref):
    hit = query.with_row_index("_i").join(
        ref.select(_EXACT).unique(), on=_EXACT, how="semi")["_i"].to_numpy()
    m = np.zeros(len(query), bool)
    m[hit] = True
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", nargs="+", required=True,
                    help="one or more v7 TSVs for the new series (concatenated, label==1 kept)")
    ap.add_argument("--ref", nargs="+", default=None,
                    help=f"reference v7 TSVs (label==1 kept). Default: the six miRBench files.")
    ap.add_argument("--merged-out", default=None,
                    help="write all input positives here (v7 schema)")
    ap.add_argument("--novel-out", default=None,
                    help="write positives NOT overlapping a same-miRNA reference site here")
    args = ap.parse_args()

    inputs = [p for pat in args.input for p in sorted(glob.glob(pat)) or [pat]]
    refs_spec = args.ref if args.ref is not None else DEFAULT_REFS
    refs = [p for pat in refs_spec for p in sorted(glob.glob(pat)) or [pat]]

    query_full = pl.concat([_read_positives(p) for p in inputs], how="vertical_relaxed")
    N = len(query_full)
    q = query_full.select(_KEYCOLS[:-1])
    print(f"input positives : {N}  ({len(inputs)} file(s))")
    print(f"  distinct miRNA seqs {q['noncodingRNA'].n_unique()}  families {q['noncodingRNA_fam'].n_unique()}")
    if args.merged_out:
        query_full.write_csv(args.merged_out, separator="\t")
        print(f"  -> {args.merged_out}")

    # Each reference file becomes its own named column of results, keyed on its basename
    # stripped to the study token (Manakov/Hejret/Klimentova collapse across splits).
    def tag(path):
        base = path.split("/")[-1]
        for t in ("Manakov", "Hejret", "Klimentova"):
            if t in base:
                return t.lower()
        return base.replace(".tsv", "")

    by_tag = {}
    for p in refs:
        by_tag.setdefault(tag(p), []).append(p)

    res = {}
    print(f"\nreference sets : {', '.join(by_tag)}\n")
    for name, paths in by_tag.items():
        ref = pl.concat([_read_positives(p) for p in paths], how="vertical_relaxed")
        any_ov, mir_ov = _overlap_masks(q, ref)
        exact = _exact_mask(q, ref)
        res[name] = (exact, mir_ov, any_ov)
        print(f"{name:<12} n={len(ref):>8}   exact(coord+miRNA) {exact.sum():>7}"
              f"   overlap+miRNA {mir_ov.sum():>7}   overlap any {any_ov.sum():>7}")

    print()
    labels = [("exact coords + miRNA", 0),
              ("site overlap + same miRNA  <- quote this", 1),
              ("site overlap (any miRNA)", 2)]
    for lvl, ix in labels:
        union = np.zeros(N, bool)
        for name in res:
            union |= res[name][ix]
        print(f"[{lvl}]")
        for name in res:
            c = res[name][ix].sum()
            print(f"    seen in {name:<12} {c:>7}  ({100 * c / N:5.1f}%)")
        new = N - union.sum()
        print(f"    seen in ANY  {union.sum():>7}   ->  NEW {new:>7}  ({100 * new / N:5.1f}% of {N})\n")

    if args.novel_out:
        union_mir = np.zeros(N, bool)
        for name in res:
            union_mir |= res[name][1]
        query_full.filter(pl.Series(~union_mir)).write_csv(args.novel_out, separator="\t")
        print(f"wrote {args.novel_out}  ({int((~union_mir).sum())} rows, "
              f"novel under 'site overlap + same miRNA')")


if __name__ == "__main__":
    main()
