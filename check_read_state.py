"""Decide whether a chimeric-eCLIP FASTQ holds RAW reads or ALREADY-TRIMMED reads.

Run this before stage 2. GEO's protocol text is not reliable on this point: two series
using the same protocol can deposit reads in different states, and the failure is silent
— trimming already-trimmed reads strips 10 nt off the 5' end, which is exactly where the
miRNA and its seed live, and nothing errors.

Both observed cases have been MEASURED with this script (400k reads each), not assumed:

                            Manakov GSE198250      GSE297116
                            (SRR18281063 R1)       (SRR33558257)
  read lengths              151, 100% uniform      variable, 105-111
  3' adapter                80.5% of reads         ~0.001% of reads
  5' U-spike (miRNA start)  position 10            position 1
  miRNA match offset        peak at 10-11          peak at 1-2
  5' 10-mer diversity       151480/200000 uniq     (insert, not a randomer)
  VERDICT                   RAW                    ALREADY TRIMMED

Four independent signals, none of which relies on the metadata:

  1. Read lengths. Raw Illumina reads are all the SAME length. A spread of lengths, and
     especially a floor exactly at cutadapt's -m, means trimming already happened.
  2. 3' adapter. Present in a decent fraction of raw reads; absent if already trimmed.
  3. Where the 5'-U spike sits. Mature miRNAs overwhelmingly start with U (AGO loading
     prefers a 5' uridine), and in a chimera the miRNA is at the 5' end of the insert.
     So the position of the T-fraction spike IS the length of whatever precedes the
     insert: ~0-2 if trimmed, ~UMI length if not. In Manakov the spike is at position 10
     (55.9% T) sitting behind a 10 nt randomer; in GSE297116 it is at position 1.
     NOTE: do NOT test this by asking whether the first 10 nt are compositionally FLAT.
     An earlier version of this script did, on the theory that a random UMI is ~25% each
     base. Real randomers are not flat — Manakov's is C-skewed (up to 37% C, mean A/C/G/T
     spread 15.9 points), so the flatness test called "no UMI" on reads that demonstrably
     have one. It is a randomer nonetheless: 151480 distinct 10-mers in 200000 reads,
     115806 of them singletons, the commonest accounting for 0.009% of reads. Diversity,
     not flatness, is what makes a UMI a UMI.
  4. Where the miRNA sits, by sequence. Same logic as 3 but matching mature miRNA
     prefixes outright instead of leaning on the U-bias. Needs the miRNA fasta.

Signals 3 and 4 can abstain (a library with no clear U-spike, or no miRNA matches at all);
the verdict is a majority of whichever signals were informative. Signal 5 below never
votes — it reports whether the UMI is recoverable from the read names, which HybriDetector
needs (is_umi=TRUE) to deduplicate PCR duplicates.
"""
import collections
import gzip
import sys

L = 18   # miRNA prefix length used for matching; long enough to be specific
W = 15   # how far into the read we look for the insert start
UMI_GUESS = 10  # prefix length whose diversity we report (the chim-eCLIP randomer)


def read_mirna_prefixes(path):
    prefixes, seq = set(), []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    s = "".join(seq).upper().replace("U", "T")
                    if len(s) >= L:
                        prefixes.add(s[:L])
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        s = "".join(seq).upper().replace("U", "T")
        if len(s) >= L:
            prefixes.add(s[:L])
    return prefixes


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: check_read_state.py <fastq.gz> <mature_mirna.fa> [n_reads]")
    fastq, mirna_fa = sys.argv[1], sys.argv[2]
    n_reads = int(sys.argv[3]) if len(sys.argv) > 3 else 400_000

    prefixes = read_mirna_prefixes(mirna_fa)
    lens = collections.Counter()
    comp = [collections.Counter() for _ in range(W)]
    mir_at = collections.Counter()
    prefix_counts = collections.Counter()
    n_adapter = 0
    names_with_umi = 0
    total = 0
    name = None

    opener = gzip.open if fastq.endswith(".gz") else open
    with opener(fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i // 4 >= n_reads:
                break
            if i % 4 == 0:
                name = line.strip()
            elif i % 4 == 1:
                r = line.strip()
                total += 1
                lens[len(r)] += 1
                if "AGATCGGAAGAGC" in r:
                    n_adapter += 1
                for p in range(min(W, len(r))):
                    comp[p][r[p]] += 1
                if len(r) >= UMI_GUESS:
                    prefix_counts[r[:UMI_GUESS]] += 1
                for k in range(W):
                    if r[k:k + L] in prefixes:
                        mir_at[k] += 1
                        break
                # umi_tools extract writes the UMI into the read name as NAME_<UMI>
                head = name.split()[0]
                if "_" in head and head.rsplit("_", 1)[1].strip("ACGTN") == "":
                    names_with_umi += 1

    if not total:
        sys.exit("no reads read")

    modal = lens.most_common(1)[0]
    frac_modal = modal[1] / total
    frac_adapter = n_adapter / total
    t_frac = [comp[p]["T"] / max(1, sum(comp[p].values())) for p in range(W)]
    spike = max(range(W), key=lambda p: t_frac[p])
    peak = max(mir_at, key=lambda k: mir_at[k]) if mir_at else None

    # Each signal votes "raw", "trimmed", or abstains (None).
    v_len = "raw" if frac_modal > 0.95 else "trimmed"
    v_adapter = "raw" if frac_adapter > 0.005 else "trimmed"
    if t_frac[spike] < 0.35:
        v_spike = None                      # no clear 5'-U spike: uninformative
    elif spike >= 5:
        v_spike = "raw"
    elif spike <= 2:
        v_spike = "trimmed"
    else:
        v_spike = None
    if peak is None or not (peak >= 5 or peak <= 2):
        v_mir = None
    else:
        v_mir = "raw" if peak >= 5 else "trimmed"

    votes = [v for v in (v_len, v_adapter, v_spike, v_mir) if v is not None]
    n_raw = votes.count("raw")
    n_trim = votes.count("trimmed")
    verdict = "raw" if n_raw > n_trim else "trimmed" if n_trim > n_raw else "ambiguous"
    umi_in_names = names_with_umi > total * 0.5

    print(f"file            : {fastq}")
    print(f"reads scanned   : {total}")
    print()
    print("1. read lengths")
    print(f"   min={min(lens)}  max={max(lens)}  modal={modal[0]} ({100*frac_modal:.1f}% of reads)")
    print(f"   -> {'UNIFORM (raw)' if v_len == 'raw' else 'VARIABLE (already trimmed)'}")
    print()
    print("2. 3' Illumina adapter (AGATCGGAAGAGC)")
    print(f"   present in {n_adapter} reads ({100*frac_adapter:.3f}%)")
    print(f"   -> {'PRESENT (raw)' if v_adapter == 'raw' else 'ABSENT (already trimmed)'}")
    print()
    print("3. 5'-U spike: mature miRNAs start with U, and the miRNA starts the insert,")
    print("   so the T spike marks where the insert begins (= length of any 5' prefix)")
    print("   per-position T%: " + " ".join(f"{100*t:.0f}" for t in t_frac))
    print(f"   max T at position {spike} ({100*t_frac[spike]:.1f}% T)")
    if v_spike is None:
        print("   -> NO CLEAR SPIKE: uninformative, abstaining")
    elif v_spike == "raw":
        print(f"   -> insert starts at {spike}: a {spike} nt prefix (the UMI) is still there (raw)")
    else:
        print(f"   -> insert starts at {spike}: reads begin at the miRNA (already trimmed)")
    uniq = len(prefix_counts)
    n_pref = sum(prefix_counts.values())
    singl = sum(1 for v in prefix_counts.values() if v == 1)
    top = prefix_counts.most_common(1)[0] if prefix_counts else ("", 0)
    print(f"   (diversity of the first {UMI_GUESS} nt: {uniq} distinct in {n_pref} reads, "
          f"{singl} singletons, commonest {100*top[1]/max(1,n_pref):.3f}%")
    print("    — high diversity = a randomer/UMI; a fixed barcode would repeat. Note that a")
    print("    real randomer is NOT compositionally flat: synthesis skew is normal.)")
    print()
    print("4. offset at which reads start matching a mature miRNA")
    if peak is None:
        print("   no miRNA matches found (is this really a chimeric eCLIP library?)")
    else:
        for k in sorted(mir_at):
            if mir_at[k] > 0:
                bar = "#" * int(40 * mir_at[k] / mir_at[peak])
                print(f"   offset {k:2d}: {mir_at[k]:6d} {bar}")
        print(f"   -> peak at offset {peak}: {peak} nt precede the miRNA")
    if v_mir is None:
        print("   -> uninformative, abstaining")
    print()
    # NB this is NOT a raw-vs-trimmed signal, and it does not vote. In a RAW read the UMI
    # sits in the SEQUENCE (the 5' randomer); it only reaches the read NAME once stage 2's
    # `umi_tools extract` moves it there. So bare names are exactly what a raw library
    # looks like at this point, and mean nothing. What IS_UMI should be therefore depends
    # on the verdict, not on this count alone.
    print("5. where is the UMI, and what should IS_UMI be?")
    print(f"   {names_with_umi}/{total} read names carry a trailing _<UMI>")
    if umi_in_names:
        print("   -> in the NAMES: these reads have already been through umi_tools extract.")
        print("      Run HybriDetector with IS_UMI=TRUE.")
    elif verdict == "raw":
        print(f"   -> not in the names, but that is EXPECTED for raw reads: the UMI is still in")
        print(f"      the SEQUENCE (the {spike} nt 5' randomer). Stage 2 moves it to the name.")
        print("      Run stage 2, then HybriDetector with IS_UMI=TRUE. PCR dedup works.")
    else:
        print("   -> not in the names, and the reads are already trimmed, so the UMI was")
        print("      excised and DISCARDED: it is unrecoverable. IS_UMI=FALSE, no PCR dedup.")
    print()

    print("=" * 72)
    print(f"votes: raw={n_raw}  trimmed={n_trim}  (abstained={4 - len(votes)})")
    if verdict == "raw":
        print("VERDICT: RAW reads. Run stage 2 (02_preprocess_fastq.sh) before HybriDetector,")
        print("         then stage 3 with IS_UMI=TRUE (stage 2 puts the UMI in the read name).")
    elif verdict == "trimmed":
        print("VERDICT: ALREADY TRIMMED. SKIP stage 2 — feed these straight to stage 3.")
        print("         Running stage 2 would cut 10 nt off the 5' end, destroying the miRNA.")
        if not umi_in_names:
            print("         UMIs are not recoverable -> IS_UMI=FALSE, and Nunique will be NA.")
    else:
        print("VERDICT: AMBIGUOUS — inspect the signals above by hand before proceeding.")
    print("=" * 72)


if __name__ == "__main__":
    main()
