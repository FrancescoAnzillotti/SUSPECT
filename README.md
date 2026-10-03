# RAFT

RAFT is a fast **pre-screening** tool for association testing when samples are related to each other (e.g. GWAS).

Instead of fitting a Linear Mixed Model (LMM) for every candidate, RAFT computes a cheap weighted Z-score for each one and keeps only the promising candidates. An exact LMM ([GEMMA](https://github.com/genetics-statistics/GEMMA)) can then be run on the survivors only. RAFT is also built for the **subset problem**: when you want to analyze only some of the samples, the expensive work from the full run is reused instead of recomputed.

This repository contains the code used in the paper:

| File | What it is |
|---|---|
| `z_score_general.c` | The RAFT implementation (C). |
| `run_subsets_saige.py` | Python wrapper that runs the benchmark: full run, then many random sample subsets, optionally against SAIGE. |
| `Makefile` | Makefile to automatically compile the C code. |
---

## How it works

1. **Standardize** every candidate (column mean 0, sd 1). Candidates with zero variance are dropped.
2. **Kinship matrix** `K = (1/m) · W Wᵀ`, which says how related each pair of samples is.
3. **Propensity scores**: `p = mean(y) + K (K + δI)⁻¹ (y − mean(y))`. This is the part of the outcome explained by relatedness (the LMM "BLUP"). `δ = 1` and `p` is clipped to `[0.01, 0.99]`.
   - Approximate mode solves this with a single Cholesky solve (no eigendecomposition).
   - Exact mode computes the eigendecomposition, because GEMMA needs it.
4. **Surprise scores** (weights): `1/p` for cases and `1/(1−p)` for controls, normalized within each class.
5. **Weighted Z-score** for every candidate, from four weighted counts (with a +0.5 correction). Candidates with `|Z| > threshold` are kept.
6. *(Exact mode only)* GEMMA is run on the kept candidates, reusing the eigendecomposition.

The Z threshold is a Bonferroni-style value that depends on the number of candidates, and is scaled by `sqrt(n' / n)` when only `n'` of the `n` samples are analyzed.

### Subset runs

The first run on a dataset is the **reference (full) run**. It saves three files in its folder:

- `specs.txt`: number of samples/candidates and the threshold
- `C_COMPUTED_kinship_matrix.bin`: the full kinship matrix
- `C_COMPUTED_propensity.bin`: the full-run propensity scores

A **subset run** with these files present does not rebuild `K`. It only cuts out `K[S,S]` for the selected samples and then either:

- recomputes the weights with one Cholesky solve (default), or
- looks up the weights from the full run (`--reuse`, fastest).

A subset is defined by setting the outcome of the removed samples to `NA`. The genotype matrix is never changed.

---

## Requirements

- Linux (the code uses `mmap`, and the wrapper uses `taskset`)
- gcc with OpenMP
- OpenBLAS and LAPACKE (`cblas.h`, `lapacke.h`)
- Python 3.8+ (standard library only)
- A pre-compiled [GEMMA](https://github.com/genetics-statistics/GEMMA) binary (see below)
- *Only for the SAIGE comparison:* [pixi](https://pixi.sh) with SAIGE installed, plus `bgzip` and `tabix` (htslib)

On Ubuntu/Debian:

```bash
sudo apt install build-essential libopenblas-dev liblapacke-dev
```

## Build

```bash
make
```

## Folder layout

The C program calls the provided GEMMA binary as `../../gemma` from inside each run folder. So the `gemma` binary must sit **two levels above the run folder**, i.e. in the directory that contains `results/`:

```
.
├── gemma                  <- GEMMA binary 
├── z_score_general        <- compiled RAFT
├── run_subsets_saige.py
├── genomics_data/
│   ├── geno.subset.txt
│   └── pheno.subset.txt
└── results/               <- created by the wrapper
    ├── FULL/
    └── remove10_cases25_rep1/ ...
```

---

## Input files

**Genotype / candidate matrix** (`-m`): one row per candidate, one column per sample. Separators can be spaces, tabs or commas.

- With `-g` (BIMBAM, GEMMA style): `id  allele1  allele2  v_1 ... v_N`
- Without `-g`: `id  v_1 ... v_N`

Values are 0, 1, 2 (dosages in between are rounded) or `NA`.

**Outcome** (`-o`): one value per sample, in the same order as the matrix columns.

- Default: `0` / `1`
- With `-2`: `1` / `2` (GEMMA coding)
- Anything else (e.g. `NA`) means the sample is excluded.

---

## Quick start

### 1. Run RAFT on one dataset

```bash
mkdir -p results/FULL
cd results
../z_score_general FULL 0 -m /abs/path/geno.txt -o /abs/path/pheno.txt -g -2
```

Usage:

```
./z_score_general <folder> <0|1> -m <matrix> -o <outcome> [-g] [-2] [--reuse]
```

| Argument | Meaning |
|---|---|
| `<folder>` | Existing folder where outputs go. Relative `-m` / `-o` paths are looked up inside it. |
| `0` / `1` | `0` = approximate (RAFT alone), `1` = exact (RAFT, then GEMMA on the survivors). |
| `-g` | Matrix is in BIMBAM format (two allele columns). |
| `-2` | Outcome is coded 1/2 instead of 0/1. |
| `--reuse` | On a subset run, take the weights from the full run instead of recomputing them. |

### 2. Run RAFT on a subset by hand

From inside `results/`, link the three reference files into a new folder and pass an outcome file where removed samples are `NA`:

```bash
mkdir subset1
for f in specs.txt C_COMPUTED_kinship_matrix.bin C_COMPUTED_propensity.bin; do
  ln -s "$PWD/FULL/$f" subset1/$f
done
../z_score_general subset1 0 -m /abs/path/geno.txt -o /abs/path/pheno_with_NA.txt -g -2 --reuse
```

### 3. Run the whole benchmark with the wrapper

```bash
python3 run_subsets_saige.py \
  --matrix  ./genomics_data/geno.subset.txt \
  --pheno   ./genomics_data/pheno.subset.txt \
  --binary  ./z_score_general \
  --workdir ./workdir \
  --results ./results \
  -g -2 \
  --removal 10,30,50,70 \
  --removed-cases 25,50,75 \
  --replicates 3 \
  --seed 1 \
  --cores 10 \
  --reuse
```

The wrapper runs the full dataset first, then every combination of removal percentage × removed-case percentage × replicate. Each run is pinned to `--cores` CPU cores.

| Option | Meaning |
|---|---|
| `--removal` | % of samples removed, comma-separated (one run set per value). |
| `--removed-cases` | Of the removed samples, the % that are cases. Comma-separated. |
| `--replicates` | Random repetitions per combination. |
| `--seed` | Seed for choosing the removed samples. |
| `--reuse` | Pass `--reuse` to the C binary (reuse full-run weights). |
| `--saige` | Also run SAIGE (Step 1 + Step 2) on every subset. Needs `-g`. |
| `--saige-pixi-dir` | Folder with the `pixi.toml` that has SAIGE installed. |
| `--saige-alpha` | Significance level for the SAIGE Bonferroni threshold (default 0.05). |
| `--keep-saige-inputs` | Keep the generated VCF/PLINK files. |

To include SAIGE, add:

```bash
  --saige --saige-pixi-dir /path/to/your/pixi/project
```

The wrapper converts each subset to the formats SAIGE needs and runs the official scripts from scratch. Only the SAIGE analysis time is reported as `saige_time_ms`. Input conversion is reported separately as `saige_prepare_time_ms`.

---

## Outputs

Inside each run folder (e.g. `results/FULL`):

| File | Content |
|---|---|
| `z_scores.txt` | Weighted Z-score of every valid candidate. |
| `output/gemma_output.assoc.txt` | GEMMA-only results (the baseline). |
| `exact/` | Exact mode only: pruned data, eigen-decomposition and GEMMA results on the survivors. |
| `specs.txt`, `*.bin` | Reference-run caches (see "Subset runs"). |

Inside the directory where the binary is launched (`results/` when using the wrapper), three logs are appended to:

- `relevant_log.csv`: how many candidates RAFT and GEMMA flag
- `time_log.csv`: running times in milliseconds
- `fn_log.csv`: candidates GEMMA flags but RAFT discards (false negatives)

The wrapper merges everything into **`results/report.csv`** and **`results/report.tsv`**, one row per experiment. The main columns:

| Column | Meaning |
|---|---|
| `zscore_retained_percentage` | % of candidates RAFT keeps. |
| `recall`, `false_negative_rate` | How many of GEMMA's hits RAFT keeps / loses. |
| `gemma_time_ms`, `zscore_time_ms` | GEMMA alone vs. RAFT alone. |
| `zscore_plus_gemma_time_ms` | RAFT followed by GEMMA (exact mode only). |
| `zscore_speedup_vs_gemma` | `gemma_time_ms / zscore_time_ms`. |
| `saige_time_ms`, `saige_accuracy` | SAIGE time, and % of GEMMA's hits that SAIGE also flags (with `--saige`). |

---

## Notes

- Always start from a clean reference folder. If `specs.txt` already exists in a folder, the program treats it as a subset run. (The wrapper deletes stale reference files for you.)
- Significance thresholds (both RAFT's Z threshold and the GEMMA `-log10(p)` threshold used for comparison) are set automatically from the number of candidates.
- The full-run memory use is roughly a few copies of the `M × N` matrix as doubles, so large `M` needs a large machine.