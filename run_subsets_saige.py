#!/usr/bin/env python3

# Example:
# python3 run_subsets_saige.py \
#   --matrix ./genomics_data/geno.subset.txt \
#   --pheno ./genomics_data/pheno.subset.txt \
#   --binary ./z_score_general \
#   --workdir ./workdir \
#   --results ./results \
#   -g \
#   -2 \
#   --removal 10,20,30,40 \
#   --removed-cases 25,50,75 \
#   --replicates 3 \
#   --seed 1 \
#   --cores 10 \
#   --reuse \
#   --saige \
#   --saige-pixi-dir /path/to/your/pixi/project

import argparse
import csv
import math
import os
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

MISSING_VALUES = {"", "NA", "NaN", "NAN", "nan", ".", "?"}

# Synthetic coordinates for the benchmark VCF/BIM.
CHROM, START_POS = "1", 1

# Files the C binary leaves in the FULL folder; subset folders link to them.
REFERENCE_FILES = (
    "specs.txt",
    "C_COMPUTED_kinship_matrix.bin",
    "C_COMPUTED_propensity.bin",
)

LOG_KEYS = (
    "zscore_relevant_snps",
    "zscore_processed_snps",
    "zscore_percentage",
    "gemma_relevant_snps",
    "gemma_tested_snps",
    "gemma_percentage",
    "gemma_after_zscore_relevant_snps",
    "gemma_after_zscore_tested_snps",
    "gemma_time_ms",
    "zscore_time_ms",
    "gemma_after_zscore_time_ms",
    "zscore_plus_gemma_time_ms",
)

TIME_KEYS = {
    "GEMMA": "gemma_time_ms",
    "Z_SCORE": "zscore_time_ms",
    "GEMMA_EXACT_AFTER_ZSCORE": "gemma_after_zscore_time_ms",
    "ZSCORE+GEMMA": "zscore_plus_gemma_time_ms",
}

# PLINK BED 2-bit codes for BIMBAM dosage 0 (A2/A2), 1 (A1/A2), 2 (A1/A1), missing.
PLINK_CODES = {0: 0b11, 1: 0b10, 2: 0b00, None: 0b01}
VCF_GENOTYPES = ("0/0", "0/1", "1/1")


# ============================================================
# INPUT & PREP
# ============================================================

def read_pheno(path):
    with open(path) as f:
        return f.read().split()


def read_matrix(path, genotype):
    n_meta = 3 if genotype else 1
    rows = []

    with open(path) as f:
        for line_number, line in enumerate(f, start=1):
            if not line.strip():
                continue

            # BIMBAM files commonly have a trailing comma.
            x = [val.strip() for val in line.rstrip().rstrip(",").split(",")]

            if len(x) <= n_meta:
                raise ValueError(
                    f"{path}:{line_number}: expected {n_meta} ID/allele "
                    "column(s) followed by genotype values."
                )

            rows.append((x[:n_meta], x[n_meta:]))

    return rows


def labels(pheno, pheno12):
    low = "1" if pheno12 else "0"
    high = "2" if pheno12 else "1"
    return [0 if x == low else 1 if x == high else None for x in pheno]


def write_matrix(rows, matrix_out):
    with open(matrix_out, "w") as f:
        for meta, values in rows:
            f.write("\t".join(meta + values) + "\n")


def write_pheno(pheno, keep, pheno_out):
    """Full-length outcome file; samples not in `keep` are written as NA."""
    kept = set(keep)

    with open(pheno_out, "w") as f:
        for i, value in enumerate(pheno):
            f.write((value if i in kept else "NA") + "\n")


# ============================================================
# RUN C
# ============================================================

def run_c(binary, run_dir, matrix, pheno, genotype, pheno12, exact, cwd,
          n_cores=10, reuse=False):
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "taskset", "-c", f"0-{n_cores - 1}",
        str(binary), str(run_dir), str(exact),
        "-m", str(matrix),
        "-o", str(pheno),
    ]

    if genotype:
        cmd.append("-g")

    if pheno12:
        cmd.append("-2")

    if reuse:
        cmd.append("--reuse")

    env = {
        **os.environ,
        "OPENBLAS_NUM_THREADS": str(n_cores),
        "OMP_NUM_THREADS": str(n_cores),
    }

    start = time.perf_counter()
    result = subprocess.run(cmd, cwd=cwd, env=env)
    elapsed = time.perf_counter() - start

    return result.returncode, elapsed


# ============================================================
# SAIGE INPUT FORMAT
# ============================================================

def parse_dosage(value, snp_id, sample_index):
    value = value.strip()

    if value in MISSING_VALUES:
        return None

    try:
        dosage = float(value)
    except ValueError as exc:
        raise ValueError(
            f"{snp_id}: sample {sample_index}: invalid dosage '{value}'"
        ) from exc

    if not 0.0 <= dosage <= 2.0:
        raise ValueError(
            f"{snp_id}: sample {sample_index}: dosage {dosage} is outside [0,2]"
        )

    return dosage


def best_guess_from_dosage(dosage):
    # Nearest hard genotype: <0.5 -> 0, 0.5..<1.5 -> 1, >=1.5 -> 2.
    if dosage is None:
        return None
    return 0 if dosage < 0.5 else 1 if dosage < 1.5 else 2


def write_saige_pheno(pheno, keep, output_path):
    """Binary trait as 0/1; missing phenotypes as NA."""
    with open(output_path, "w") as f:
        f.write("IID\tPHENO\n")

        for sample_number, original_index in enumerate(keep, start=1):
            label = pheno[original_index]
            f.write(f"sample_{sample_number}\t{'NA' if label is None else label}\n")


def write_plink_bed_bim_fam(rows, keep, output_prefix):
    """
    Hard-call PLINK 1 files, used ONLY for SAIGE Step 1 (GRM construction).
    Step 2 uses the original dosages through the VCF DS field.
    """
    n_samples = len(keep)

    with open(f"{output_prefix}.bed", "wb") as bed, \
         open(f"{output_prefix}.bim", "w") as bim, \
         open(f"{output_prefix}.fam", "w") as fam:

        # PLINK BED magic bytes + SNP/variant-major mode.
        bed.write(bytes([0x6C, 0x1B, 0x01]))

        for sample_number in range(1, n_samples + 1):
            iid = f"sample_{sample_number}"
            fam.write(f"{iid}\t{iid}\t0\t0\t0\t-9\n")

        for variant_index, (meta, values) in enumerate(rows):
            snp_id, allele1, allele2 = meta
            bim.write(
                f"{CHROM}\t{snp_id}\t0\t{START_POS + variant_index}\t"
                f"{allele1}\t{allele2}\n"
            )

            codes = [
                PLINK_CODES[
                    best_guess_from_dosage(parse_dosage(values[i], snp_id, s))
                ]
                for s, i in enumerate(keep, start=1)
            ]

            # Four samples per byte; the last byte is zero-padded.
            for j in range(0, n_samples, 4):
                bed.write(bytes([
                    sum(code << (2 * k) for k, code in enumerate(codes[j:j + 4]))
                ]))


def write_saige_id_file(rows, output_path):
    """Variant IDs expected by SAIGE for VCF input: chrom:pos:REF:ALT."""
    with open(output_path, "w") as f:
        for variant_index, (meta, _) in enumerate(rows):
            f.write(f"{CHROM}:{START_POS + variant_index}:{meta[2]}:{meta[1]}\n")


def write_vcf(rows, keep, output_path):
    """
    VCF whose DS field counts BIMBAM allele1 (REF = allele2, ALT = allele1).
    """
    with open(output_path, "w") as out:
        out.write(
            "##fileformat=VCFv4.3\n"
            "##source=BIMBAM_to_SAIGE_benchmark\n"
            '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage of ALT allele">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Best-guess genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
            + "".join(f"\tsample_{s}" for s in range(1, len(keep) + 1))
            + "\n"
        )

        for variant_index, (meta, values) in enumerate(rows):
            snp_id, allele1, allele2 = meta
            out.write(
                f"{CHROM}\t{START_POS + variant_index}\t{snp_id}\t"
                f"{allele2}\t{allele1}\t.\tPASS\t.\tGT:DS"
            )

            for sample_number, i in enumerate(keep, start=1):
                dosage = parse_dosage(values[i], snp_id, sample_number)

                if dosage is None:
                    out.write("\t./.:.")
                else:
                    gt = VCF_GENOTYPES[best_guess_from_dosage(dosage)]
                    out.write(f"\t{gt}:{dosage:.8f}")

            out.write("\n")


def compress_and_index_vcf(vcf_path):
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} was not found in PATH. Install htslib first.")

    gz_path = Path(str(vcf_path) + ".gz")

    with gz_path.open("wb") as out:
        subprocess.run(
            ["bgzip", "-c", str(vcf_path)],
            stdout=out, stderr=subprocess.PIPE, check=True,
        )

    subprocess.run(
        ["tabix", "--csi", "-p", "vcf", str(gz_path)],
        capture_output=True, check=True,
    )

    vcf_path.unlink()
    csi_path = Path(str(gz_path) + ".csi")

    if not csi_path.exists():
        raise RuntimeError(
            f"tabix completed but did not create expected index: {csi_path}"
        )

    return gz_path, csi_path


def prepare_saige_inputs(rows, pheno, keep, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = output_dir / "genotypes.vcf"
    pheno_path = output_dir / "pheno.txt"
    plink_prefix = output_dir / "grm_genotypes"

    start = time.perf_counter()

    write_vcf(rows, keep, vcf_path)
    vcf_gz, vcf_csi = compress_and_index_vcf(vcf_path)
    write_plink_bed_bim_fam(rows, keep, plink_prefix)
    write_saige_pheno(pheno, keep, pheno_path)

    return {
        "vcf": vcf_gz,
        "vcf_index": vcf_csi,
        "plink_prefix": plink_prefix,
        "pheno": pheno_path,
        "prepare_time_ms": (time.perf_counter() - start) * 1000.0,
    }


# ============================================================
# SAIGE DISCOVERY / EXECUTION
# ============================================================

def run_pixi_r(pixi_dir, r_args, n_cores, log_path):
    """Launch `pixi run Rscript ...` pinned to cores 0..n_cores-1."""
    cmd = ["taskset", "-c", f"0-{n_cores - 1}", "pixi", "run", "Rscript"] + r_args

    env = {
        **os.environ,
        "OPENBLAS_NUM_THREADS": str(n_cores),
        "OMP_NUM_THREADS": str(n_cores),
    }

    with open(log_path, "w") as log:
        start = time.perf_counter()
        result = subprocess.run(
            cmd, cwd=pixi_dir, env=env, stdout=log, stderr=subprocess.STDOUT
        )
        elapsed = time.perf_counter() - start

    return result.returncode, elapsed


def discover_saige(pixi_dir):
    """
    Get the SAIGE version through the pixi environment and locate the
    Step 1 / Step 2 wrapper scripts (they live in <pixi_dir>/extdata/ in the
    standard checkout, not necessarily in the installed R package).
    """
    result = subprocess.run(
        [
            "pixi", "run", "R", "--vanilla", "--slave", "-e",
            'library(SAIGE); cat(as.character(packageVersion("SAIGE"))); '
            'cat("\\n"); cat(find.package("SAIGE")); cat("\\n")',
        ],
        cwd=pixi_dir, capture_output=True, text=True,
    )

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]

    if result.returncode != 0 or len(lines) < 2:
        raise RuntimeError(
            "Could not load SAIGE through 'pixi run R'.\n\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )

    version, package_dir = lines[0], Path(lines[1])

    def find_script(name):
        for d in (pixi_dir / "extdata", package_dir / "extdata", pixi_dir):
            if (d / name).exists():
                return (d / name).resolve()
        return next(pixi_dir.rglob(name), None)

    step1_script = find_script("step1_fitNULLGLMM.R")
    step2_script = find_script("step2_SPAtests.R")

    if step1_script is None or step2_script is None:
        raise RuntimeError(
            f"SAIGE {version} loaded, but step1_fitNULLGLMM.R / "
            f"step2_SPAtests.R were not found under {pixi_dir} "
            f"or {package_dir}."
        )

    return version, step1_script.resolve(), step2_script.resolve()


def run_saige(rows, pheno, keep, run_dir, tmp_saige_dir, saige_pixi_dir,
              step1_script, step2_script, n_cores, p_threshold):
    """
    Official SAIGE Step 1 (null model, GRM from hard-call PLINK) and
    Step 2 (single-variant tests on VCF DS), run through `pixi run Rscript`.
    """
    run_dir.mkdir(parents=True, exist_ok=True)
    tmp_saige_dir.mkdir(parents=True, exist_ok=True)

    inputs = prepare_saige_inputs(rows, pheno, keep, tmp_saige_dir)
    print(f"    SAIGE input preparation: {inputs['prepare_time_ms'] / 1000:.2f}s")

    step1_prefix = run_dir / "saige_step1"
    step2_output = run_dir / "saige_step2.txt"
    step2_id_file = tmp_saige_dir / "saige_step2_ids.txt"
    write_saige_id_file(rows, step2_id_file)

    model_file = Path(str(step1_prefix) + ".rda")
    variance_ratio_file = Path(str(step1_prefix) + ".varianceRatio.txt")

    def run_step(name, args, log):
        print(f"    SAIGE {name} starting...")
        rc, elapsed = run_pixi_r(saige_pixi_dir, args, n_cores, log)
        print(f"    SAIGE {name} completed in {elapsed:.2f}s")

        if rc != 0:
            raise RuntimeError(
                f"SAIGE {name} failed for {run_dir.name} "
                f"with exit code {rc}. See {log}"
            )

        return rc, elapsed

    step1_rc, step1_elapsed = run_step("Step 1", [
        str(step1_script),
        f"--plinkFile={inputs['plink_prefix']}",
        f"--phenoFile={inputs['pheno']}",
        "--phenoCol=PHENO",
        "--sampleIDColinphenoFile=IID",
        "--traitType=binary",
        f"--outputPrefix={step1_prefix}",
        f"--nThreads={n_cores}",
        "--LOCO=FALSE",
        "--minMAFforGRM=0.01",
        "--tol=0.02",
        "--maxiter=20",
        "--tolPCG=1e-5",
        "--maxiterPCG=500",
        "--numRandomMarkerforVarianceRatio=30",
        "--IsOverwriteVarianceRatioFile=TRUE",
        "--isCateVarianceRatio=FALSE",
    ], run_dir / "saige_step1.log")

    for path in (model_file, variance_ratio_file):
        if not path.exists():
            raise RuntimeError(f"SAIGE Step 1 completed but {path} is missing")

    step2_rc, step2_elapsed = run_step("Step 2", [
        str(step2_script),
        f"--vcfFile={inputs['vcf']}",
        f"--vcfFileIndex={inputs['vcf_index']}",
        "--vcfField=DS",
        f"--SAIGEOutputFile={step2_output}",
        f"--chrom={CHROM}",
        "--minMAF=0.0",
        "--minMAC=0.5",
        "--maxMissing=0.15",
        "--impute_method=best_guess",
        "--LOCO=FALSE",
        f"--GMMATmodelFile={model_file}",
        f"--varianceRatioFile={variance_ratio_file}",
        "--markers_per_chunk=10000",
        "--is_output_moreDetails=FALSE",
        "--is_overwrite_output=TRUE",
        "--SPAcutoff=2.0",
        f"--idstoIncludeFile={step2_id_file}",
        f"--nThreads={n_cores}",
    ], run_dir / "saige_step2.log")

    if not step2_output.exists():
        raise RuntimeError(
            f"SAIGE Step 2 completed but output file is missing: {step2_output}"
        )

    saige_elapsed_ms = (step1_elapsed + step2_elapsed) * 1000.0

    # Compare against GEMMA's association output generated by the C benchmark.
    # Input conversion is deliberately outside the timed region above.
    gemma_assoc = find_gemma_assoc_file(run_dir)

    if gemma_assoc is None:
        raise RuntimeError(
            f"Could not find GEMMA *.assoc.txt output in {run_dir}; "
            "cannot compute saige_accuracy."
        )

    saige_hits = parse_significant_snps(
    step2_output,
    ("marker", "rs", "snp", "snpid", "markerid"),
    "p.value",
    p_threshold,
    )

    gemma_hits = parse_significant_snps(
        gemma_assoc,
        ("rs",),
        "p_wald",
        p_threshold,
    )

    saige_relevant_snps = len(saige_hits)

    saige_percentage = (
        saige_relevant_snps / len(rows) * 100.0
        if rows
        else None
    )

    return {
    "saige_prepare_time_ms": inputs["prepare_time_ms"],
    "saige_time_ms": saige_elapsed_ms,
    "saige_total_time_ms": (
        inputs["prepare_time_ms"] + saige_elapsed_ms
    ),
    "saige_relevant_snps": saige_relevant_snps,
    "saige_percentage": saige_percentage,
    "saige_accuracy": (
        len(saige_hits & gemma_hits) / len(gemma_hits) * 100.0
        if gemma_hits else None
    ),
    "saige_step1_return_code": step1_rc,
    "saige_step2_return_code": step2_rc,
}


def find_gemma_assoc_file(run_dir):
    """Find the GEMMA association result produced by the C benchmark."""
    for path in sorted(run_dir.rglob("*.assoc.txt")):
        try:
            with open(path) as f:
                header = [x.strip().lower() for x in f.readline().split("\t")]
        except OSError:
            continue

        if "p_wald" in header and "rs" in header:
            return path

    return None


def parse_significant_snps(path, id_names, p_name, p_threshold):
    """SNP IDs whose p-value column (matched by substring) is <= p_threshold."""
    hits = set()

    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")

        try:
            header = [x.strip().lower() for x in next(reader)]
        except StopIteration:
            raise RuntimeError(f"Output is empty: {path}")

        id_index = next(
            (header.index(name) for name in id_names if name in header), None
        )
        p_index = next((i for i, col in enumerate(header) if p_name in col), None)

        if id_index is None or p_index is None:
            raise RuntimeError(
                f"Could not find SNP ID ({'/'.join(id_names)}) or {p_name} "
                f"column in {path}: {header}"
            )

        for fields in reader:
            if len(fields) <= max(id_index, p_index):
                continue

            snp_id = fields[id_index].strip()

            try:
                p = float(fields[p_index])
            except ValueError:
                continue

            if snp_id and math.isfinite(p) and p <= p_threshold:
                hits.add(snp_id)

    return hits


# ============================================================
# PARSE GLOBAL LOGS
# ============================================================

def read_tsv_lines(path):
    if path.exists():
        with open(path) as f:
            for line in f:
                yield line.strip().split("\t")


def parse_global_logs(results_dir):
    data = {}

    def get_row(folder):
        if folder not in data:
            data[folder] = {**dict.fromkeys(LOG_KEYS), "false_negatives_total": 0}
        return data[folder]

    for parts in read_tsv_lines(results_dir / "relevant_log.csv"):
        if len(parts) < 4:
            continue

        try:
            row = get_row(Path(parts[0]).resolve())
            method = parts[1]

            if method == "ZSCORE":
                row["zscore_relevant_snps"] = int(parts[2])
                row["zscore_processed_snps"] = int(parts[3])
                row["zscore_percentage"] = float(parts[4].rstrip("%"))

            elif method == "GEMMA":
                row["gemma_relevant_snps"] = int(parts[2])
                row["gemma_tested_snps"] = int(parts[3])
                row["gemma_percentage"] = float(parts[4].rstrip("%"))

            elif method == "GEMMA_AFTER_ZSCORE":
                row["gemma_after_zscore_relevant_snps"] = int(parts[2])
                row["gemma_after_zscore_tested_snps"] = int(parts[3])

        except (ValueError, IndexError):
            continue

    for parts in read_tsv_lines(results_dir / "time_log.csv"):
        if len(parts) < 3:
            continue

        try:
            value = float(parts[2])
            row = get_row(Path(parts[0]).resolve())

            if parts[1] in TIME_KEYS:
                row[TIME_KEYS[parts[1]]] = value

        except (ValueError, IndexError):
            continue

    for parts in read_tsv_lines(results_dir / "fn_log.csv"):
        if len(parts) >= 2:
            get_row(Path(parts[0]).resolve())["false_negatives_total"] += 1

    return data


# ============================================================
# COMPUTE EXPERIMENT METRICS
# ============================================================

def ratio(num, den, scale=1.0):
    """num / den * scale, or None when undefined."""
    return None if num is None or not den else num / den * scale


def compute_metrics(params, log_data, input_snps):
    row = {**params, **{k: log_data.get(k) for k in LOG_KEYS}}
    row["false_negatives_total"] = log_data.get("false_negatives_total", 0)

    gemma_rel = row["gemma_relevant_snps"]
    fn = row["false_negatives_total"]
    zscore_proc = row["zscore_processed_snps"]
    zscore_rel = row["zscore_relevant_snps"]
    gemma_ms = row["gemma_time_ms"]

    row["false_negative_rate"] = ratio(fn, gemma_rel, 100.0)
    row["recall"] = (
        ratio(gemma_rel - fn, gemma_rel, 100.0) if gemma_rel else None
    )

    row["zscore_retained_percentage"] = ratio(zscore_rel, zscore_proc, 100.0)

    row["variance_retention_percentage"] = ratio(zscore_proc, input_snps, 100.0)
    row["variance_removed_percentage"] = (
        None if row["variance_retention_percentage"] is None
        else 100.0 - row["variance_retention_percentage"]
    )

    row["zscore_speedup_vs_gemma"] = ratio(gemma_ms, row["zscore_time_ms"])
    row["total_pipeline_speedup_vs_gemma"] = ratio(
        gemma_ms, row["zscore_plus_gemma_time_ms"]
    )

    row["gemma_candidate_retention_percentage"] = ratio(
        zscore_rel, row["gemma_tested_snps"], 100.0
    )
    row["gemma_candidate_reduction_percentage"] = (
        None if row["gemma_candidate_retention_percentage"] is None
        else 100.0 - row["gemma_candidate_retention_percentage"]
    )

    return row


# ============================================================
# MAIN
# ============================================================

def main():
    p = argparse.ArgumentParser()

    p.add_argument("--matrix", required=True)
    p.add_argument("--pheno", required=True)
    p.add_argument("--binary", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--results", required=True)

    p.add_argument("--removal", default="")
    p.add_argument("--removed-cases", default="")
    p.add_argument("--replicates", type=int, default=1)
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("-g", "--genotype-file", action="store_true")
    p.add_argument("-2", "--pheno12", action="store_true")

    p.add_argument(
        "--cores", type=int, default=10,
        help="CPU cores for C and SAIGE (cores 0..N-1 are pinned with taskset).",
    )
    p.add_argument(
        "--reuse", action="store_true",
        help="Pass --reuse to the C binary: subset runs load the FULL run's "
             "propensities instead of recomputing them.",
    )

    p.add_argument(
        "--saige", action="store_true",
        help="Run SAIGE Step 1 + Step 2 for every experiment.",
    )
    p.add_argument(
        "--saige-pixi-dir", default=".",
        help="Directory containing the pixi.toml used for SAIGE (default: .).",
    )
    p.add_argument(
        "--saige-alpha", type=float, default=0.05,
        help="Significance level for Bonferroni correction (default: 0.05).",
    )
    p.add_argument(
        "--keep-saige-inputs", action="store_true",
        help="Keep generated SAIGE VCF/PLINK/phenotype inputs in each result directory.",
    )

    args = p.parse_args()

    matrix = Path(args.matrix).resolve()
    pheno_file = Path(args.pheno).resolve()
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    results = Path(args.results).resolve()
    saige_pixi_dir = Path(args.saige_pixi_dir).resolve()

    workdir.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)

    available = len(os.sched_getaffinity(0))

    if not 1 <= args.cores <= available:
        raise ValueError(f"--cores must be between 1 and {available}, got {args.cores}")

    if not 0 < args.saige_alpha < 1:
        raise ValueError("--saige-alpha must be between 0 and 1")

    if shutil.which("taskset") is None:
        raise RuntimeError("taskset was not found in PATH.")

    # Clean legacy global logs so the C binary does not append to old data.
    for log_name in ("relevant_log.csv", "time_log.csv", "fn_log.csv"):
        (results / log_name).unlink(missing_ok=True)

    pheno = read_pheno(pheno_file)
    rows = read_matrix(matrix, args.genotype_file)

    if any(len(values) != len(pheno) for _, values in rows):
        raise ValueError("Matrix and phenotype file have different sample counts.")

    y = labels(pheno, args.pheno12)

    cases = [i for i, value in enumerate(y) if value == 1]
    controls = [i for i, value in enumerate(y) if value == 0]

    n = len(pheno)
    input_snps = len(rows)
    saige_p_threshold = args.saige_alpha / input_snps

    saige_step1_script = saige_step2_script = None

    if args.saige:
        if not args.genotype_file:
            raise ValueError(
                "--saige requires --genotype-file/-g because SAIGE "
                "needs SNP alleles and genotype input."
            )

        if shutil.which("pixi") is None:
            raise RuntimeError("--saige was requested, but 'pixi' was not found in PATH.")

        if not saige_pixi_dir.exists():
            raise RuntimeError(f"SAIGE pixi directory does not exist: {saige_pixi_dir}")

        saige_version, saige_step1_script, saige_step2_script = (
            discover_saige(saige_pixi_dir)
        )

        print(f"SAIGE version: {saige_version}")
        print(f"SAIGE Step 1 script: {saige_step1_script}")
        print(f"SAIGE Step 2 script: {saige_step2_script}")

    print(
        f"Samples: {n} | Cases: {len(cases)} | "
        f"Controls: {len(controls)} | SNPs: {input_snps}"
    )

    all_experiments = []

    with tempfile.TemporaryDirectory(dir=workdir, prefix="subset_tmp_") as tmp:
        tmpdir = Path(tmp)

        def run_saige_experiment(name, run_dir, keep):
            print(f"  {name}: running SAIGE...")

            saige_input_dir = (
                run_dir / "saige_inputs"
                if args.keep_saige_inputs
                else tmpdir / f"saige_{name}"
            )

            metrics = run_saige(
                rows=rows,
                pheno=y,
                keep=keep,
                run_dir=run_dir,
                tmp_saige_dir=saige_input_dir,
                saige_pixi_dir=saige_pixi_dir,
                step1_script=saige_step1_script,
                step2_script=saige_step2_script,
                n_cores=args.cores,
                p_threshold=saige_p_threshold,
            )

            accuracy = metrics["saige_accuracy"]

            print(
                f"    SAIGE: "
                f"{metrics['saige_time_ms'] / 1000:.2f}s | "
                f"{metrics['saige_relevant_snps']} SNPs | "
                f"{metrics['saige_percentage']:.2f}% of {input_snps} total SNPs | "
                f"accuracy: "
                f"{'NA' if accuracy is None else f'{accuracy:.2f}%'}"
            )

            return metrics

        # The matrix is written once. Every run, FULL or subset, reads all of
        # it; subsets differ only in which outcome values are NA.
        full_matrix = tmpdir / "full_matrix.txt"
        write_matrix(rows, full_matrix)

        # =====================================================
        # FULL DATASET
        # =====================================================

        print("\n=== FULL DATASET ===")

        full_dir = results / "FULL"
        full_pheno = tmpdir / "full_pheno.txt"
        write_pheno(pheno, range(n), full_pheno)

        # Stale reference files would make the C binary treat FULL as a subset run.
        for file in REFERENCE_FILES:
            (full_dir / file).unlink(missing_ok=True)

        rc, elapsed = run_c(
            binary, full_dir, full_matrix, full_pheno,
            args.genotype_file, args.pheno12, 0,
            cwd=results, n_cores=args.cores,
        )

        if rc != 0:
            raise RuntimeError(f"Full run failed with exit code {rc}")

        if not (full_dir / "specs.txt").exists():
            raise RuntimeError("C program did not produce FULL/specs.txt")

        full_saige_metrics = (
            run_saige_experiment("FULL", full_dir, list(range(n)))
            if args.saige else {}
        )

        all_experiments.append({
            "experiment": "FULL",
            "folder": full_dir.resolve(),
            "removal_pct": 0.0,
            "removed_case_pct": 0.0,
            "replicate": 0,
            "total_samples": n,
            "removed_samples": 0,
            "removed_cases": 0,
            "removed_controls": 0,
            "remaining_samples": n,
            "input_snps": input_snps,
            "python_wall_time_seconds": elapsed,
            "return_code": rc,
            **full_saige_metrics,
        })

        print(f"Full C run completed in {elapsed:.2f}s")

        # =====================================================
        # REMOVAL EXPERIMENTS
        # =====================================================

        if args.removal and args.removed_cases:
            rng = random.Random(args.seed)

            removals = [float(x) for x in args.removal.split(",")]
            removed_case_pcts = [float(x) for x in args.removed_cases.split(",")]

            for removal_pct in removals:
                n_remove = round(n * removal_pct / 100)

                print(f"\n=== REMOVE {removal_pct:g}% ({n_remove} samples) ===")

                for case_pct in removed_case_pcts:
                    n_remove_cases = round(n_remove * case_pct / 100)
                    n_remove_controls = n_remove - n_remove_cases

                    if n_remove_cases > len(cases) or n_remove_controls > len(controls):
                        continue

                    for rep in range(1, args.replicates + 1):
                        name = f"remove{removal_pct:g}_cases{case_pct:g}_rep{rep}"
                        run_dir = results / name

                        removed = set(
                            rng.sample(cases, n_remove_cases)
                            + rng.sample(controls, n_remove_controls)
                        )
                        keep = [i for i in range(n) if i not in removed]

                        # Removed samples become NA in the outcome file.
                        pheno_out = tmpdir / f"{name}_pheno.txt"
                        write_pheno(pheno, keep, pheno_out)

                        # Link the FULL run's specs and caches so the C binary
                        # recognises this as a subset run (see its SUBSET REUSE note).
                        run_dir.mkdir(parents=True, exist_ok=True)

                        for file in REFERENCE_FILES:
                            (run_dir / file).unlink(missing_ok=True)
                            (run_dir / file).symlink_to(full_dir / file)

                        rc, elapsed = run_c(
                            binary, run_dir, full_matrix, pheno_out,
                            args.genotype_file, args.pheno12, 0,
                            cwd=results, n_cores=args.cores, reuse=args.reuse,
                        )

                        saige_metrics = (
                            run_saige_experiment(name, run_dir, keep)
                            if args.saige else {}
                        )

                        all_experiments.append({
                            "experiment": name,
                            "folder": run_dir.resolve(),
                            "removal_pct": removal_pct,
                            "removed_case_pct": case_pct,
                            "replicate": rep,
                            "total_samples": n,
                            "removed_samples": n_remove,
                            "removed_cases": n_remove_cases,
                            "removed_controls": n_remove_controls,
                            "remaining_samples": len(keep),
                            "input_snps": input_snps,
                            "python_wall_time_seconds": elapsed,
                            "return_code": rc,
                            **saige_metrics,
                        })

                        print(f"  {name} C completed in {elapsed:.2f}s (code {rc}) \n\n\n")

        else:
            print(
                "\nNo --removal / --removed-cases provided "
                "— running FULL dataset only."
            )

    # =========================================================
    # PARSE GLOBAL LOGS AND COMBINE METRICS
    # =========================================================

    print("\nMerging data and generating reports...")

    parsed_logs = parse_global_logs(results)

    final_results = []

    for exp in all_experiments:
        folder = exp.pop("folder")
        final_results.append(
            compute_metrics(exp, parsed_logs.get(folder, {}), exp["input_snps"])
        )

    # =========================================================
    # EXPORT RESULTS
    # =========================================================

    preferred_columns = [
        # Experiment design
        "experiment", "removal_pct", "removed_case_pct", "replicate",
        "total_samples", "removed_samples", "removed_cases", "removed_controls",
        "remaining_samples", "input_snps",

        # Z-score
        "zscore_processed_snps", "zscore_relevant_snps", "zscore_percentage",
        "zscore_retained_percentage", "variance_retention_percentage",
        "variance_removed_percentage",

        # GEMMA
        "gemma_tested_snps", "gemma_relevant_snps", "gemma_percentage",
        "gemma_after_zscore_tested_snps", "gemma_after_zscore_relevant_snps",
        "gemma_candidate_retention_percentage",
        "gemma_candidate_reduction_percentage",

        # C / GEMMA accuracy
        "false_negatives_total", "false_negative_rate", "recall",

        # Times
        "gemma_time_ms", "zscore_time_ms", "gemma_after_zscore_time_ms",
        "zscore_plus_gemma_time_ms", "zscore_speedup_vs_gemma",
        "total_pipeline_speedup_vs_gemma",

        ## SAIGE
        "saige_prepare_time_ms",
        "saige_time_ms",
        "saige_total_time_ms",
        "saige_relevant_snps",
        "saige_percentage",
        "saige_accuracy",
        "saige_step1_return_code",
        "saige_step2_return_code",

        # Python/C status
        "python_wall_time_seconds", "return_code",
    ]

    fieldnames = [
        col for col in preferred_columns
        if any(col in row for row in final_results)
    ]

    for report, delimiter in (("report.csv", ","), ("report.tsv", "\t")):
        with open(results / report, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(final_results)

    print("========================================")
    print("Done.")
    print(f"CSV report: {results / 'report.csv'}")
    print(f"TSV report: {results / 'report.tsv'}")
    print("========================================")


if __name__ == "__main__":
    main()