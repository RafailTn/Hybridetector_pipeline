# Building miRBench-style datasets from a new chimeric eCLIP experiment

This turns raw chimeric-eCLIP/CLASH/CLEAR-CLIP reads into a labelled v7 TSV — the same
18-column `v7` schema the miRBench datasets use (listed [below](#the-four-stages)).

It is a thin local driver over two upstream repos, cloned into `external/` (gitignored):

| Repo | Role |
| --- | --- |
| [ML-Bioinfo-CEITEC/HybriDetector](https://github.com/ML-Bioinfo-CEITEC/HybriDetector) (`fix_clustering` branch) | reads → chimeric miRNA:target hybrids |
| [BioGeMT/miRBench_paper](https://github.com/BioGeMT/miRBench_paper) | hybrids → labelled train/test datasets |

The papers: [HybriDetector / Hejret 2023](https://doi.org/10.1038/s41598-023-49757-z)
and [miRBench / Sammut, Gresova et al., *Bioinformatics* 2025](https://academic.oup.com/bioinformatics/article/41/Supplement_1/i542/8199406).
The `fix_clustering` branch is the one the miRBench paper used, not `main`.

## Where it lives — a subfolder of a larger repo, or a flat clone of its own

The scripts locate their working root from their own path, so this pipeline runs **both**
as a subfolder of a larger repo and cloned on its own:

- **As a subfolder of a larger repo** (in a `chimeric_eclip/` directory) the root is the
  *parent* — so `external/` (the upstream clones) and `data/` (the conservation bigwigs and
  the dataset outputs) are the parent-level directories, **shared** with the rest of that repo.
- **Cloned flat** — the scripts at a repo's top level, no `chimeric_eclip/` subfolder — the
  scripts' own directory is the root, so the clone is self-contained: `external/` and `data/`
  live inside it, and nothing outside is needed.

So the commands below are written for the subfolder layout (`bash chimeric_eclip/setup.sh`, …).
**In a flat clone, drop the `chimeric_eclip/` prefix** — `bash setup.sh`,
`bash 01_download_geo.sh …`, and so on. The scripts' own "next step" hints already print the
right path for however you invoked them.

Set `REPO_ROOT=/path` to force the root explicitly. The one layout that misfires is cloning
the *flat* repo into a directory named exactly `chimeric_eclip`: the scripts then assume the
subfolder layout and reach one level up for `external/`/`data/`. Name the clone
anything else, or pass `REPO_ROOT`.

## Setup (once)

```bash
bash chimeric_eclip/setup.sh --prebuild-hd-envs
```

Clones both repos, applies `patches/hybridetector-fixes.patch`, and builds four
micromamba environments. HybriDetector is a Snakemake workflow where each rule builds its
own conda env at runtime (`--use-conda`), which needs a real conda/mamba solver — hence
micromamba rather than a single pinned environment.

| Env | Used by | Contents |
| --- | --- | --- |
| `eclip_dl` | stage 1 | sra-tools, entrez-direct, GEOparse |
| `eclip_pp` | stage 2 | umi_tools, cutadapt |
| `hybridetector` | stage 3 | snakemake 7.18.2 + mamba (the 13 per-rule envs — STAR, R/Bioconductor, ViennaRNA — are built on top) |
| `mirbench_pp` | stage 4 | pandas, pyBigWig, pyranges, DECIPHER, genomic_region_annotator |

> **Updating an existing clone.** `setup.sh` applies `hybridetector-fixes.patch` only to a
> *clean* HybriDetector checkout — it detects an already-patched (dirty) working tree and
> skips, so it will **not** carry a later patch revision onto a clone you set up earlier.
> When the patch changes (e.g. the `onstart` sparsity guard was added to it), refresh the
> clone yourself: `git -C external/HybriDetector stash` (or `checkout .`) then re-run
> `setup.sh`, or just delete `external/HybriDetector` and let `setup.sh` re-clone it.

## The four stages

```bash
# 1. raw reads. --gse writes an editable sample sheet; --samplesheet downloads it.
#    Only IP libraries carry chimeras — drop the (input) controls.
bash chimeric_eclip/01_download_geo.sh --gse GSE297116 --out data/raw/myset
bash chimeric_eclip/01_download_geo.sh --samplesheet <sheet> --out data/raw/myset

# 1b. ALWAYS DO THIS. Are the reads raw, or did the submitter already trim them?
bash chimeric_eclip/check_read_state.sh data/raw/myset/<sample>.fastq.gz

# 2. trim: 5' UMI -> read header, 3' adapter, 3' UMI remnant.
#    SKIP THIS ENTIRELY if 1b says the reads are already trimmed.
bash chimeric_eclip/02_preprocess_fastq.sh data/raw/myset data/pp/myset

# 3. call chimeras (the long one: hours, and see the memory note below)
bash chimeric_eclip/03_run_hybridetector.sh data/pp/myset

# 4. hybrids -> labelled v7 dataset
bash chimeric_eclip/04_make_dataset.sh external/HybriDetector/hyb_pairs data/myset
```

## Check the read state before you trim anything

**Do not trust the GEO protocol text.** Two series can run the same protocol and deposit
reads in *different states*, and getting this wrong fails silently rather than loudly:
trimming already-trimmed reads strips 10 nt off the 5' end, which is exactly where the
miRNA and its seed sit. Nothing errors — you just lose the chimeras.

| | Manakov (GSE198250) | GSE297116 |
| --- | --- | --- |
| deposited reads | **raw** | **already preprocessed** |
| 5' UMI | in the sequence (10 nt) | pruned before upload |
| 3' adapter | present | trimmed (1 read in 400k retains it) |
| 3' UMI remnant | present | none |
| UMI recoverable? | yes | **no** — not in the read names, not a technical read in SRA |
| stage 2 | run it | **skip it** |
| stage 3 | `IS_UMI=TRUE` | `IS_UMI=FALSE` |

`check_read_state.sh` decides this from the reads themselves, on four independent signals
(uniform vs variable read length; adapter presence; the position of the 5'-U spike, since
mature miRNAs start with U and the miRNA starts the insert, so the T-spike marks where the
insert begins; and the offset at which reads start matching a mature miRNA, which *is* the
length of whatever still precedes it). Signals may abstain, and the majority wins, so one
weak signal cannot flip the call.

It is validated in both directions: run on Manakov (`SRR18281063`) it returns **RAW 4–0**,
run on GSE297116 (`SRR33558257`) **ALREADY TRIMMED 4–0**.

> Do **not** test for a UMI by asking whether the first 10 nt are compositionally *flat*.
> Real randomers are not flat — Manakov's is C-skewed (up to 37% C) — and that test calls
> "no UMI" on reads that carry one. Diversity is the discriminator: Manakov's first 10 nt
> yield 262,918 distinct 10-mers in 400,000 reads (176,722 singletons). An earlier version
> of this script made exactly this mistake.

### Datasets without UMIs

If the UMIs are unrecoverable, HybriDetector cannot remove PCR duplicates, and with
`is_umi=FALSE` it emits only `Ndups` — **no `Nunique` column at all** — which miRBench's
`filtering.py` indexes directly. `patches/mirbench-nunique-na.patch` makes `Nunique` fall
back to **NA**, in both `filtering.py` and `make_neg_sets.py` (which otherwise hardcodes
`Nunique = 0` for negatives, which would leave a column that is NA for one class and 0
for the other: populated-looking, and a perfect label separator). The patch is a no-op on
the canonical datasets, which all carry UMIs.

How much does the missing deduplication actually cost? Less than it sounds:

- `Nunique` is a **pure passthrough**. It is never used in any filter, threshold, cluster,
  split, or negative-sampling decision — it is copied through and zeroed for negatives.
  (In the Manakov training set *every* negative has `Nunique = 0` and every positive has
  `Nunique > 0`: it is a perfect label leak, not a feature — don't train on it.)
- HybriDetector's filters select on alignment length, alignment quality, RepeatMasker
  overlap and base composition — **no read-count thresholds**; "high confidence" means a
  mismatch-free miRNA alignment, not depth.
- PCR duplicates are byte-identical reads, so they collapse into the *same* chimera row
  and merely sum into its count.

So the set of (miRNA, target) pairs called is essentially unaffected; what is lost is the
read-support number. Don't compare `Nunique` across datasets or pool it with Manakov's.

The final files land in `data/myset/step6_add_conservation/` and carry the v7 schema:

```
gene  noncodingRNA  noncodingRNA_name  noncodingRNA_fam  feature  label
chr  start  end  strand  Nunique  dominant_region  regions_present
read_start_in_sel_tx_1based  read_end_in_sel_tx_1based  gene_cluster_ID
gene_phyloP  gene_phastCons
```

The `*.train.*.conservation.tsv` and `*.test.*.conservation.tsv` files are the finished
product — a labelled, train/test-split dataset ready to hand to whatever model you train.

## How much of a new series is actually new?

Before treating a new series' positives as fresh training signal, ask how many are just
re-observations of interactions the existing datasets already contain (Manakov alone is
~1.4M positives). `positive_overlap.py` answers that:

```bash
micromamba run -n mirbench_pp python chimeric_eclip/positive_overlap.py \
    --input data/myset/step6_add_conservation/*.conservation.tsv \
    --merged-out data/myset_positives_merged_v7.tsv \
    --novel-out  data/myset_positives_novel_v7.tsv
```

`--input` takes the step-6 v7 files (label==1 kept); `--ref` defaults to the six miRBench
files (Manakov train/test/leftout, Hejret train/test, Klimentova test) but takes any
baseline you pass. It reports three definitions of "already seen", per reference dataset
and combined, then writes the novel subset:

| definition | what it requires |
| --- | --- |
| exact coords + miRNA | identical chr/start/end/strand **and** same miRNA — almost always an undercount, since HybriDetector anchors its 50 nt window on the read, so two experiments rarely produce byte-identical intervals |
| **site overlap + same miRNA** | the target intervals overlap **and** the miRNA matches — **the one to quote**, robust to that window jitter |
| site overlap, any miRNA | the target intervals overlap, regardless of miRNA |

The gap between the last two is the interesting part: sites where a *different* miRNA hits
the *same* AGO2-bound locus — novelty that lives in the assignment, not the target. For
GSE304955 (SAEC), 30% of positives were same-miRNA re-hits of a Manakov site while 73.5%
were the same locus with a different miRNA, so `--novel-out` kept ~105k of 150k. `--novel-out`
is written under the middle definition. Both output files hold **positives only** — the
matching negatives are dropped, so it is a novelty audit, not a ready-to-train set.

## Things worth knowing before you run it

**Protocol assumptions live in stage 2.** The defaults are chim-eCLIP (Manakov):
a 10 nt 5′ UMI and a 10 nt 3′ UMI remnant. A CLASH or CLEAR-CLIP library has a
different barcode layout, and getting this wrong is silent — HybriDetector is
"very sensitive to every single base difference in the input reads", so a stray
barcode base does not error, it just quietly costs you chimeras. Override
`UMI_PATTERN` / `UMI3_LEN` / `MIN_LEN` per protocol.

**Memory — and why the cluster path exists.** This is the one place where *where you run*
changes *what you get*.

A dense hg38 STAR index needs ~30 GB resident, both to build and to load for every
alignment. On the 31 GB laptop it is not merely slow, it is **impossible**: STAR was
OOM-killed at 27.9 GB resident during the `packing SA` step (kernel memcg log), after
5 CPU-hours. The only way to run locally is a **sparse suffix array**
(`STAR_SA_SPARSE_D=2`), which halves index RAM — 16 GB instead of an OOM.

**Sparse is not free, and I measured it.** Dense and sparse produce *different*
alignments. On chr21, with 164,171 real chr21-mapping reads and HybriDetector's own
alignment parameters, counting the distinct **reads** whose placement changes (name, flag,
chrom, pos, CIGAR, NH):

| comparison | differing reads |
| --- | --- |
| dense vs dense (same index, twice) | 1,224 (0.75%) — the noise floor, from `--outMultimapperOrder Random` |
| dense vs **sparse** | **6,931 (4.22%) — 5.7× the noise floor** |

So the difference is real, not run-to-run jitter. (The count is per *read*, not per SAM
record: a multimapper whose alignment changes emits several differing records, so a naive
record diff over-counts by 3–4× here — `validate_sparse_index.sh` collapses to distinct
reads. The within-chr21 estimate is well-sampled, but chr21 is ~1.5% of the genome, so
the exact genome-wide magnitude is not established — enough to refute neutrality, not
to pin the number.)

Consequence: **miRBench built Manakov's index dense.** A sparse local run is therefore
*not* mapping-identical to Manakov. If that matters for your comparison, use the cluster.

**The dense/sparse choice is not a tracked Snakemake input.** `STAR_gen_index` reads the
sparsity from the environment, but the rule declares only the genome FASTA as input — so if
an `index/STAR/` already exists, Snakemake reuses it *whatever sparsity you now ask for*,
and the run's `DENSE`/`SPARSE` banner would then describe the request, not the bytes on disk
(the trap: a leftover sparse index reused under a `DENSE` claim, silently not comparable to
Manakov). This is guarded in two places, both of which read the `genomeSAsparseD` STAR
recorded in `index/STAR/genomeParameters.txt` and **abort** on a mismatch, telling you to
`rm -rf .../index/STAR` to rebuild or to set `STAR_SA_SPARSE_D` to match the existing index:
`03_run_hybridetector*.sh` check before launching, and — since Snakemake itself won't catch
it — the patched Snakefile carries an `onstart` handler (`hybridetector-fixes.patch`) that
raises the identical error, so calling `snakemake` on HybriDetector directly is guarded too.

The Snakefile and STAR wrapper now **default to upstream** — dense index, `mem` 200/50/34 —
so nothing about this laptop follows the pipeline anywhere else. `03_run_hybridetector.sh`
**sizes itself to the host it runs on**: it reads `MemTotal`, budgets 90% of it, and takes
the dense path at ≥40 GB and the sparse path below that. It prints which one it chose. So
the *same command* is right on the laptop and on a big server; only the machine changes.

## Running somewhere bigger

**On a plain server (no scheduler)** — nothing special. Same command as local; it detects
the RAM and goes dense on its own:

```bash
bash chimeric_eclip/setup.sh --prebuild-hd-envs
bash chimeric_eclip/01_download_geo.sh --samplesheet chimeric_eclip/samplesheets/<sheet>.tsv --out data/raw/myset
IS_UMI=FALSE bash chimeric_eclip/03_run_hybridetector.sh data/raw/myset
```

It will announce `-> budget NNN GB, DENSE STAR index`. Override `MEM_GB` if you are sharing
the box with someone and shouldn't take 90% of it; override `CORES` (default `nproc`) likewise.

**On a SLURM cluster** — use the scheduler runner, which submits one `sbatch` job per rule:

```bash
$EDITOR chimeric_eclip/slurm/config.yaml     # partition=compute -> your partition
bash chimeric_eclip/setup.sh --prebuild-hd-envs
bash chimeric_eclip/03_run_hybridetector_slurm.sh data/raw/myset
```

It exports **neither** `STAR_SA_SPARSE_D` nor `HD_MEM_*`, so the workflow runs stock upstream.
`slurm/config.yaml` maps the Snakefile's `resources: mem` (GB) and `threads` onto
`sbatch --mem` / `--cpus-per-task`, and asks for a 64 GB slot for the genome index and 48 GB
for the alignments (upstream's literal `mem = 200` is a scheduling token, not a real
requirement — STAR needs ~32 GB). `latency-wait` is 120 s because shared filesystems lag
behind job completion; do not lower it.

Either way, stages 1, 2 and 4 are cheap and unchanged — only stage 3 needs the big machine.

**First run of stage 3 is expensive.** It downloads the Ensembl release-90 primary
assembly (~900 MB) and builds ~10 STAR indices — hours, and ~40 GB of disk, once.
The `.snakemake/conda` envs add ~4 GB.

**Which HybriDetector output is the right one.** Stage 4 consumes
`*.unified_length_all_types_unique_high_confidence.tsv` — the one with raw column names
(`seq.g`, `noncodingRNA_seq`, `chr.g`, `Nunique`). There is a sibling
`*_high_confidence_finalout.tsv` holding the same rows under the prettier column names
documented in HybriDetector's README (`Genomic fragment sequence`, …); miRBench's
`filtering.py` does not understand those, so stage 4 always reads the raw-column file.

**How the labels are made** (stage 4, step 3) is the substance of the miRBench paper.
Negatives are not random miRNA–site shuffles: for each positive, a binding site is
sampled from a *different gene cluster*, and each miRNA family keeps the same share of
the negative class as it has of the positive class. That is what removes the miRNA
frequency-class bias — without it a model scores well by learning which miRNAs are
merely common. A consequence worth knowing: a positive and its negative
twin share a target site, so any feature that is a property of the site alone (mapped
coordinates, target conservation, accessibility) is constant within a twin pair and
cannot separate them.

**Train/test split is by chromosome**, not at random: step 0 sets `test = (chr == "1")`,
so chr1 sites are the test set. Step 2 additionally holds out miRNA families that occur
nowhere else, which is what produced the `*_leftout_v7.tsv` files.

## Deviations from upstream, and why

- **`merge_replicates` is bypassed** — multi-sample runs now take the per-sample path and are
  concatenated in stage 4. The rule is broken in *both* released branches and always has been:
  the Snakefile has declared `hyb_pairs/Merged.{type}.tsv` since commit `cc6939a`, but
  `merge_replicates.R` still writes `Merged.hybrids_deduplicated_filtered_collapsed_{type}.tsv`,
  so the job exits 0 and Snakemake fails it with `MissingOutputException`. **miRBench never used
  it**: `concat_HD_output.sh` globs the *per-sample* `*.unified_length_all_types_unique_high_confidence.tsv`
  and concatenates, with dedup happening in post-processing step 0. Taking the same path is both
  the working option and the one that reproduces Manakov v7. `HD_MERGE_REPLICATES=1` opts back in.
- **`merge_replicates.R` no-UMI crash fixed** (`Object 'Nunique' not found`). The `is_umi=FALSE`
  branch drops `Nunique`, then an unguarded `order(-Nunique)` sorts by it. Guarded to match the
  same function's other sort (`-Nunique` with UMIs, `-Ndups` without). Only reachable via
  `HD_MERGE_REPLICATES=1`, but fixed so that path isn't doubly broken.
- **cutadapt `-f fastq` dropped** (stage 2). The option was removed in cutadapt 2.0 and
  is a hard error on the 5.x in `eclip_pp`; the format is inferred anyway.
- **SLURM removed.** Upstream's preprocessing and post-processing scripts are `sbatch`
  job arrays with hardcoded `#SBATCH --account` directives. Stages 2 and 4 are local
  loops over the same commands with the same parameters.
- **`jq` added to three HybriDetector rule envs.** The post-link script of
  `bioconductor-bsgenome.hsapiens.ucsc.hg38` shells out to `jq` via `yq` but does not
  declare it as a dependency, so the `filter_and_collapse` / `unify_length` /
  `merge_replicates` envs fail to build without it.
- **`python <3.12` pinned in every HybriDetector rule env.** Snakemake 7.18's script
  runner imports `distutils`, which was removed from the stdlib in Python 3.12 (PEP 632).
  The rule envs pin no Python, so conda hands them the newest one and every `script:`
  rule dies with `ModuleNotFoundError: No module named 'distutils'`. This is latent in
  upstream — it only began to bite once conda-forge's default Python crossed 3.12 — and
  it does not show up at env-build time, only when the first rule actually runs.
- **Genome-length count fixed in `STAR_gen_index`.** It derived STAR's
  `--genomeSAindexNbases` from `grep -v '>' <genome> | wc -m`, but the whole-genome input
  is *gzipped* (the rule gunzips it later), so grep saw binary, said so on stderr, and
  `wc -m` returned 0 → `math.log(0)` → `ValueError: math domain error`. The ncRNA indices
  were unaffected because their inputs are plain `.fa`. Now reads through `zcat`, which
  measures 3.15 Gb and yields the correct value of 14 for hg38.
- **Ensembl download switched from FTP to HTTPS** (same file; FTP is blocked or slow on
  many networks).
- **`defaults` / `hcc` channels dropped** from the miRBench env (Anaconda ToS;
  everything resolves from conda-forge/bioconda). The upstream download env also pins
  ~120 exact builds and hardcodes a University of Malta `ftp_proxy`; `eclip_dl` is a
  clean rewrite.
