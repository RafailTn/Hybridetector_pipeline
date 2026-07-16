# Construction of a miRBench-style v7 dataset from GSE297116

Methods record for the AGO2 chimeric-eCLIP dataset built with `chimeric_eclip/`. It documents
what was run, in what order, with which parameters, and — separately — every point at which the
pipeline departs from the published upstream tools, with the reason. Written to be reportable:
each deviation is either a fix for a bug that made upstream fail outright, or a decision forced
by a property of this dataset that the Manakov data did not have.

---

## 1. Provenance

| component | source | version used |
| --- | --- | --- |
| Chimera calling | [HybriDetector](https://github.com/ML-Bioinfo-CEITEC/HybriDetector), branch **`fix_clustering`** | commit `dc7a355` ("final version used to process Manakov2022 data", 2024-08-20) |
| Labelling / negatives | [miRBench_paper](https://github.com/BioGeMT/miRBench_paper) | commit `6a52d87` |
| Reference genome | Ensembl release 90, `Homo_sapiens.GRCh38.dna.primary_assembly` | downloaded by HybriDetector's `prepare_ref` rule |
| Conservation | UCSC hg38 `phyloP100way`, `phastCons100way` | bigWig, from `hgdownload.cse.ucsc.edu` |
| Environments | micromamba, per-rule conda envs | see §7 for how to dump exact versions |

**`fix_clustering`, not `main`.** This is the branch miRBench used to process Manakov 2022 (the
commit message says so explicitly). `main` lacks the clustering fix. Using `main` would not
reproduce the reference dataset.

## 2. Sample selection

Series **GSE297116** — "Nuclear Argonaute:miRNA complexes recognize target sequences within
chromatin and silence gene expression". AGO2 miR-eCLIP (Eclipsebio) in HCT116.

Selected: **wild-type HCT116, whole-cell fraction, IP libraries only** — two biological replicates.

| sample | SRA runs | GEO |
| --- | --- | --- |
| `HCT116_wholecell_rep1_IP` | SRR33558257, SRR33558258 | GSM8984092 |
| `HCT116_wholecell_rep2_IP` | SRR33558254, SRR33558255 | GSM8984094 |

Each sample has two runs (a deep run plus a top-off of the same library); stage 1 concatenates
them into one FASTQ per sample.

Excluded, and why:
- **`(input)` libraries** — size-matched controls with no AGO2 pulldown, therefore no chimeras.
- **NLS-AGO2 samples** — engineered nuclear-localised AGO2. A different targetome; mixing it with
  wild-type would confound the dataset.
- **Cytoplasmic and chromatin fractions** — only the whole-cell fraction was requested.

## 3. Stage 0 — read-state determination (**this pipeline adds this step**)

Run **before** any trimming:

```bash
bash chimeric_eclip/check_read_state.sh data/raw/gse297116
```

**Finding: GSE297116 reads are already preprocessed.** Adapters are trimmed, and the UMIs have
been excised *and discarded* (not moved to the read header). This is unlike Manakov 2022, whose
SRA submission carries raw reads.

The call was **validated against Manakov as a positive control**: the same script was run,
unmodified, on 400,000 reads of `SRR18281063` R1 (GSM5942089, `Expt1_293T_NoGelChimeCLIP_total_rep1`
— an unenriched chimeric-eCLIP IP library, the closest analogue to ours in the series miRBench
built v7 from). It returns the opposite verdict on every signal, unanimously:

| signal | Manakov GSE198250 (SRR18281063 R1) | GSE297116 (SRR33558257) |
| --- | --- | --- |
| read lengths | 151 nt, **100.0% uniform** | **variable**, 18–111 nt (modal 110, 19.2%) |
| 3′ adapter `AGATCGGAAGAGC` | **80.5%** of reads | **0.000%** (1 read in 400,000) |
| 5′-U spike (= insert start) | position **10** (56.0% T) | position **0** (44.8% T) |
| mature-miRNA match offset | peak at **10–11** | peak at **1–2** |
| UMI in read names | absent | absent |
| **verdict** | **RAW (4–0)** | **ALREADY TRIMMED (4–0)** |

The four signals, none of which relies on GEO's metadata:

1. **Read lengths.** Raw Illumina output is a single uniform length (Manakov: 151 nt, 100% of
   reads). GSE297116 spreads 18–111 nt, with the floor sitting exactly on cutadapt's `-m 18` —
   the fingerprint of a trimmer having already run.
2. **3′ adapter.** An untrimmed library reads through into the adapter whenever the insert is
   shorter than the read — in Manakov, 80.5% of reads. In GSE297116 it is gone (1 read in 400,000).
3. **Position of the 5′-U spike.** Mature miRNAs overwhelmingly begin with U (AGO loading prefers a
   5′ uridine), and in a chimera the miRNA begins the insert. So the position of the T-fraction
   spike *is* the length of whatever precedes the insert. Manakov spikes at position 10 (56.0% T),
   behind its 10 nt UMI; GSE297116 spikes at position 0 (44.8% T) — the read *starts* at the miRNA.
4. **Where the miRNA sits, by sequence.** Matching mature-miRNA prefixes outright, independently of
   the U-bias: Manakov's miRNAs start at offset 10–11, GSE297116's at offset 1–2.

Corroborated at the archive level: `vdb-dump` shows a single BIOLOGICAL read per spot with no
technical read, and `fastq-dump --origfmt` yields bare names — there is nowhere a UMI could hide.

**Consequence: stage 2 (preprocessing) was SKIPPED, and HybriDetector was run with `is_umi=FALSE`.**

> **A test that does *not* work, recorded because it is the intuitive one.** Do not ask whether the
> first 10 nt are compositionally *flat* on the theory that a random UMI is ~25% per base. Real
> randomers are not flat: Manakov's is C-skewed (up to 37% C; mean A/C/G/T spread 15.9 points), so
> a flatness test reports "no UMI" on reads that demonstrably carry one. It *is* a randomer —
> 262,918 distinct 10-mers in 400,000 reads, 176,722 of them singletons, the commonest accounting
> for 0.008% of reads. **Diversity, not flatness, is what identifies a UMI.** An earlier version of
> this script used the flatness test; it happened to reach the right verdict for GSE297116 on the
> strength of the other signals, but it voted the wrong way on Manakov and has been replaced by
> signal 3 above.

**Why this check exists.** Running the standard preprocessing on already-trimmed reads is silently
destructive: `umi_tools extract` does not *detect* a UMI, it *assumes* one and unconditionally
removes the first 10 bases of every read. On these reads that amputates real sequence — precisely
where the miRNA seed lies in a chimeric read. There is no error and no warning; the run completes
with degraded chimera yield. The only symptom is an unexplained low hybrid count.

## 4. Stage 1 — download

```bash
bash chimeric_eclip/01_download_geo.sh \
     --samplesheet chimeric_eclip/samplesheets/GSE297116_wholecell_wt.tsv \
     --out data/raw/gse297116
```

`prefetch` + `fasterq-dump` (SRA Toolkit), then concatenation of the runs belonging to each sample.

## 5. Stage 2 — preprocessing: **SKIPPED** (see §3)

For a *raw* series this stage runs `umi_tools extract` (10 nt 5′ UMI → read header) → `cutadapt`
(tiled 3′ adapters) → `cutadapt -u -10` (3′ UMI remnant), reproducing miRBench's preprocessing.
It was not run here, because the submitter had already done it.

## 6. Stage 3 — chimera calling (HybriDetector)

```bash
IS_UMI=FALSE bash chimeric_eclip/03_run_hybridetector.sh data/raw/gse297116
```

Parameters (as written to `external/HybriDetector/config.json`):

| parameter | value | note |
| --- | --- | --- |
| `is_umi` | `FALSE` | forced by §3 — no recoverable UMIs |
| `read_length` | **measured**, not assumed | max post-trim read length via `seqkit stats`; see below |
| `map_perc_single_genomic` | `0.85` | upstream default |
| `map_perc_softclip` | `0.75` | upstream default |
| STAR suffix array | **dense** (`--genomeSAsparseD 1`) | STAR's default, and what miRBench used — see §9 |

**`read_length` is measured, not taken from GEO.** HybriDetector derives the minimum fraction of a
read that must align to the miRNA as `ceil(16 / read_length * 100) / 100` — i.e. "at least 16 nt of
non-coding RNA". The value must therefore be the length of the reads *HybriDetector actually sees*,
after trimming, which is not the read length the sequencer or GEO reports. Set it too high and
flimsy miRNA alignments pass; too low and genuine chimeras are rejected. The value used by this run
is recorded in `external/HybriDetector/config.json` and echoed in the run log as
`== measured read_length (post-trim max) = N`.

Output consumed downstream, **per sample**:

```
hyb_pairs/<sample>.unified_length_all_types_unique_high_confidence.tsv
```

Raw column names (`seq.g`, `noncodingRNA_seq`, `chr.g`, `Nunique`, …). Note the sibling
`*_high_confidence_finalout.tsv` holds the same rows under the *publication* column names from
HybriDetector's README; miRBench's `filtering.py` cannot read those. The raw file is the correct
input.

## 7. Stage 4 — labelling (miRBench post-processing)

```bash
bash chimeric_eclip/04_make_dataset.sh data/hyb_gse297116 data/gse297116
```

The two per-sample hybrid tables are concatenated (one header retained), then miRBench's seven-step
pipeline runs in `--mode single`:

| step | what it does |
| --- | --- |
| 0 | filter to `noncodingRNA_type == miRNA`, rename to v7 columns, drop duplicate `(gene, noncodingRNA)` sequence pairs — **this is where the replicate overlap is deduplicated** |
| 1 | annotate against Ensembl-90 transcripts → `feature`, `dominant_region`, `regions_present`, transcript coordinates |
| 2 | hold out miRNA families occurring nowhere else (the `leftout` set) |
| 3 | **negative sampling** (see below) |
| 4 | train/test split on the `test` column, set in step 0 from `chr.g == "1"` |
| 5 | drop the `test` column |
| 6 | add `gene_phyloP`, `gene_phastCons` from the bigWigs |

**Negative sampling (step 3) is the substance of the miRBench method** and should be described as
such in any write-up. Negatives are *not* random miRNA–site shuffles. For each positive, a binding
site is sampled from a **different gene cluster**, and each miRNA family is kept at the same share
of the negative class as it holds of the positive class. This removes the miRNA frequency-class
bias — without it, a classifier scores well simply by learning which miRNAs are common. A direct
consequence: a positive and its negative twin **share a target site**, so any feature that is a
property of the site alone (genomic coordinates, target conservation, accessibility) is constant
within a twin pair and cannot separate them.

**Train/test split is by chromosome, not at random**: chromosome 1 is held out as the test set.

Result: TSVs with the 18-column v7 schema in `data/gse297116/step6_add_conservation/`, directly
usable as `TRAIN=` / `TESTS=` for `cnn/run_train.sh`.

---

## 8. Deviations from upstream

Two categories, and the distinction matters for reporting. **A–E are upstream bugs**: the tools do
not run at all without them. **F–H are decisions** forced by properties of this dataset.

### Upstream bugs fixed (HybriDetector)

Patch: `chimeric_eclip/patches/hybridetector-fixes.patch`.

**A. `merge_replicates` is bypassed; multi-sample runs take the per-sample path.**
The Snakefile has declared this rule's output as `hyb_pairs/Merged.{type}.tsv` since commit
`cc6939a`, but `merge_replicates.R` still writes `Merged.hybrids_deduplicated_filtered_collapsed_{type}.tsv`.
The job exits 0 having written a file under the old name; Snakemake fails it with
`MissingOutputException`. **This is true of both released branches (`main` and `fix_clustering`) —
the rule has never worked.** Upstream never noticed because *miRBench never used it*: their
`concat_HD_output.sh` globs the **per-sample** `*.unified_length_all_types_unique_high_confidence.tsv`
files and concatenates them, letting post-processing step 0 deduplicate the replicate overlap.
This pipeline does the same. It is therefore both the working path and the one that reproduces
Manakov v7. (`HD_MERGE_REPLICATES=1` restores the upstream path for anyone who repairs it.)

**B. `unify_length`'s `merged` flag decoupled from the sample count.**
Upstream set `merged = "TRUE" if len(SAMPLES) > 1`, *independently* of the branch `all_inputs`
chose. With the per-sample path (A), `merged=TRUE` told `unify_length` its input came from
`merge_replicates` — so it selected a `rep` column that does not exist and wrote `Merged.*`
filenames the rule does not declare. Both decisions are now bound to a single constant
(`USE_MERGE_REPLICATES`) so they cannot drift apart.

**C. `merge_replicates.R` crashed on `is_umi=FALSE` (`Object 'Nunique' not found`).**
The no-UMI branch correctly drops the `Nunique` column (without UMIs there is no unique-read
count), then an unguarded `order(-Nunique)` sorts by it. Guarded to match the same function's
*other* sort, which already does the right thing (`-Nunique` with UMIs, `-Ndups` without). Only
reachable via `HD_MERGE_REPLICATES=1` after fix A, but fixed so that path is not doubly broken.

**D. `STAR_gen_index` died with `math domain error`.**
It computed genome length as `grep -v '>' <file> | wc -m`, but the whole-genome input is gzipped —
grep saw binary, printed nothing, `wc -m` returned 0, and `log(0)` raised. It only ever worked
because the small ncRNA indices take plain `.fa` input. Now read through `zcat`.

**E. Rule conda environments were unusable.**
Two independent problems. (i) No environment pinned Python, so conda resolved 3.12+, where
`distutils` no longer exists — and Snakemake 7.18's script runner imports it, so *every* `script:`
rule died with `ModuleNotFoundError`. Pinned `python <3.12` in all rule envs. Latent and nasty:
environment *creation* succeeds, so `--conda-create-envs-only` cannot detect it. (ii)
`bioconductor-bsgenome.hsapiens.ucsc.hg38`'s post-link script shells out to `jq` via `yq` without
declaring it, breaking three environments; `jq` added to those `env.yaml` files. (Installing `jq`
in the driver environment does *not* work — mamba runs post-link scripts with a sanitised `PATH`.)

Also: `prepare_ref` fetched Ensembl over **FTP**, which many networks block. Changed to HTTPS with
`wget -c` (resumable — it is a ~900 MB download).

### Upstream bugs fixed (miRBench)

Patch: `chimeric_eclip/patches/mirbench-nunique-na.patch`.

**F. `Nunique` written as `NA` when the library has no UMIs.**
With `is_umi=FALSE`, HybriDetector emits no `Nunique` column and miRBench's `filtering.py` crashes
on the missing key. `Nunique` was established to be a **pure passthrough** — it is used in no
filter, threshold, clustering, split, or negative-sampling decision anywhere in the post-processing.
It is written as **`NA`**, not `0`, deliberately: a fabricated count would be indistinguishable from
a measured one downstream. Verified a no-op on the canonical UMI-bearing datasets.

> Worth knowing when using the resulting data: in **Manakov v7**, every negative has `Nunique = 0`
> and every positive has `Nunique > 0`. It is a perfect label leak, not a feature. Any model must
> exclude it. In this dataset it is `NA` throughout, so the leak cannot arise.

### Decisions forced by this dataset

**G. Stage 2 skipped, `is_umi=FALSE`** — see §3. This is the single most consequential decision in
the build, and the one most likely to be got wrong by someone repeating it on a different series.

**H. `read_length` measured rather than assumed** — see §6.

Also: cutadapt's `-f fastq` flag (used by upstream's preprocessing) was **removed in cutadapt 2.0**
and is a hard error on current versions; dropped, since the format is inferred anyway. This only
affects the stage-2 path, which was not used here.

---

## 9. A note on STAR index sparsity (relevant if the run is ever repeated on a small machine)

A dense hg38 STAR index needs ~32 GB resident, both to build and to load for each alignment. On a
31 GB machine this is not merely slow but impossible: STAR was OOM-killed at 27.9 GB resident during
`packing SA` (kernel memcg log). The only way to run there is a sparse suffix array
(`--genomeSAsparseD 2`), which roughly halves index RAM.

**Sparse indexing is not alignment-neutral.** Measured on chr21 with 164,171 real chr21-mapping
reads and HybriDetector's own parameters, counting the distinct **reads** whose placement changes
(name, flag, chrom, pos, CIGAR, NH), with a dense-vs-dense control for the RNG noise floor:

| comparison | differing reads |
| --- | --- |
| dense vs dense (same index, run twice) | 1,224 (0.75%) — the noise floor, from `--outMultimapperOrder Random` |
| dense vs **sparse** | **6,931 (4.22%) — 5.7× the noise floor** |

(Per read, not per SAM record: a multimapper whose alignment changes emits several differing
records, so a naive record diff over-counts by 3–4×; `validate_sparse_index.sh` collapses to
distinct reads.)

**miRBench built Manakov's index dense.** A sparse run therefore does not produce mapping-comparable
output. **This dataset was built with a dense index**, on a machine with sufficient RAM, and is
comparable. `03_run_hybridetector.sh` sizes itself to the host (dense at ≥40 GB, sparse below) and
prints which it chose; the run log line
`== NNN GB host -> budget NNN GB, DENSE STAR index` is the record of that.

---

## 10. Numbers to record for the write-up

These are properties of *your* run; capture them rather than quoting anyone else's.

```bash
# hybrids called, per replicate
wc -l external/HybriDetector/hyb_pairs/*.unified_length_all_types_unique_high_confidence.tsv

# final dataset size and class balance
for f in data/gse297116/step6_add_conservation/*.tsv; do
    echo "$f: $(( $(wc -l < "$f") - 1 )) rows"
    awk -F'\t' 'NR>1 {c[$6]++} END {for (l in c) print "   label="l": "c[l]}' "$f"
done

# read_length actually used, and is_umi
cat external/HybriDetector/config.json

# exact tool versions in the environments that ran (record from the machine that ran them)
micromamba list -n hybridetector
for e in external/HybriDetector/.snakemake/conda/*/; do
    echo "== $e"; "$e/bin/conda" list 2>/dev/null | grep -E "^(star|samtools|cutadapt|umi_tools|viennarna|r-base)\s"
done
```

Two sanity checks worth doing explicitly:

- **`Nunique` must be `NA` throughout**, not `0`. If it is `0`, the no-UMI patch did not apply and
  the label leak described in §8F is present.
- **Row count should be a plausible fraction of Manakov's** (~780k rows from 19 samples). Two
  samples yielding a similar order of magnitude would indicate something wrong.

## 11. Known limitations

- **Two replicates only.** Manakov v7 draws on 19 samples; this dataset is correspondingly smaller,
  and rare miRNA families will be sparsely represented — relevant because step 2 holds out families
  that occur nowhere else, and because family is the grouping variable for `StratifiedGroupKFold`.
- **A different cell line** (HCT116) and a different laboratory from Manakov's, so batch effects are
  confounded with any biological difference. Treat cross-dataset comparisons accordingly.
- **`Nunique` is unavailable**, so any analysis conditioning on read support is not possible on this
  dataset (it is not possible on Manakov either without leaking the label — see §8F).
- **The magnitude of the sparse-vs-dense difference in §9 was quantified only on chr21**, which is
  enough to refute neutrality but not to size the whole-genome effect. Moot for this dataset, which
  was built dense.
