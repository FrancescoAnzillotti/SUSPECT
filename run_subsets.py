


#  python3 run_subsets.py     
# --matrix ./genomics_data/geno.subset.txt     
# --pheno ./genomics_data/pheno.subset.txt     
# --binary ./z_score_general     
# --workdir ./workdir     
# --results ./results     
# -g 
# -2    
# --removal 10,20,30,40     
# --removed-cases 25,50,75     
# --replicates 3     
# --seed 1



import argparse
import csv
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

# ============================================================
# INPUT & PREP
# ============================================================

def read_pheno(path):
    with open(path) as f:
        return f.read().split()

def read_matrix(path, genotype):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            x = [val.strip() for val in line.split(',')]
            meta = x[:3] if genotype else x[:1]
            values = x[3:] if genotype else x[1:]
            rows.append((meta, values))
    return rows

def labels(pheno, pheno12):
    low = "1" if pheno12 else "0"
    high = "2" if pheno12 else "1"
    return [0 if x == low else 1 if x == high else None for x in pheno]

def write_subset(rows, pheno, keep, matrix_out, pheno_out):
    matrix_out.parent.mkdir(parents=True, exist_ok=True)
    with open(matrix_out, "w") as f:
        for meta, values in rows:
            f.write("\t".join(meta + [values[i] for i in keep]) + "\n")
    with open(pheno_out, "w") as f:
        for i in keep:
            f.write(pheno[i] + "\n")

# ============================================================
# RUN C
# ============================================================

def run_c(binary, run_dir, matrix, pheno, genotype, pheno12, exact, cwd):
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), str(run_dir), str(exact), "-m", str(matrix), "-o", str(pheno)]
    if genotype:
        cmd.append("-g")
    if pheno12:
        cmd.append("-2")
    start = time.time()
    result = subprocess.run(cmd, cwd=cwd)
    return result.returncode, time.time() - start

# ============================================================
# PARSE GLOBAL LOGS
# ============================================================

def parse_global_logs(results_dir):
    data = {}
    
    def get_row(folder):
        if folder not in data:
            data[folder] = {
                "zscore_relevant_snps": None, "zscore_processed_snps": None, "zscore_percentage": None,
                "gemma_relevant_snps": None, "gemma_tested_snps": None, "gemma_percentage": None,
                "gemma_after_zscore_relevant_snps": None, "gemma_after_zscore_tested_snps": None,
                "gemma_time_ms": None, "zscore_time_ms": None,
                "gemma_after_zscore_time_ms": None, "zscore_plus_gemma_time_ms": None,
                "false_negatives_total": 0
            }
        return data[folder]

    # Parse relevant_log.csv
    rel_path = results_dir / "relevant_log.csv"
    if rel_path.exists():
        with open(rel_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4: continue
                try:
                    folder = Path(parts[0]).resolve()
                    method = parts[1]
                    row = get_row(folder)
                    
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

    # Parse time_log.csv
    time_path = results_dir / "time_log.csv"
    if time_path.exists():
        with open(time_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3: continue
                try:
                    folder = Path(parts[0]).resolve()
                    method = parts[1]
                    val = float(parts[2])
                    row = get_row(folder)
                    
                    if method == "GEMMA": row["gemma_time_ms"] = val
                    elif method == "Z_SCORE": row["zscore_time_ms"] = val
                    elif method == "GEMMA_EXACT_AFTER_ZSCORE": row["gemma_after_zscore_time_ms"] = val
                    elif method == "ZSCORE+GEMMA": row["zscore_plus_gemma_time_ms"] = val
                except (ValueError, IndexError):
                    continue

    # Parse fn_log.csv
    fn_path = results_dir / "fn_log.csv"
    if fn_path.exists():
        with open(fn_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2: continue
                try:
                    folder = Path(parts[0]).resolve()
                    get_row(folder)["false_negatives_total"] += 1
                except Exception:
                    continue

    return data

# ============================================================
# COMPUTE EXPERIMENT METRICS
# ============================================================

def compute_metrics(params, log_data, input_snps):
    row = dict(params)
    
    keys = [
        "zscore_relevant_snps", "zscore_processed_snps", "zscore_percentage",
        "gemma_relevant_snps", "gemma_tested_snps", "gemma_percentage",
        "gemma_after_zscore_relevant_snps", "gemma_after_zscore_tested_snps",
        "gemma_time_ms", "zscore_time_ms", "gemma_after_zscore_time_ms",
        "zscore_plus_gemma_time_ms"
    ]
    for k in keys:
        row[k] = log_data.get(k)
        
    row["false_negatives_total"] = log_data.get("false_negatives_total", 0)
    
    gemma_rel = row["gemma_relevant_snps"]
    fn = row["false_negatives_total"]
    if gemma_rel and gemma_rel > 0:
        row["false_negative_rate"] = fn / gemma_rel * 100.0
        row["recall"] = (gemma_rel - fn) / gemma_rel * 100.0
    else:
        row["false_negative_rate"] = row["recall"] = None

    zscore_proc = row["zscore_processed_snps"]
    zscore_rel = row["zscore_relevant_snps"]
    
    if zscore_proc and zscore_proc > 0 and zscore_rel is not None:
        row["zscore_retained_percentage"] = zscore_rel / zscore_proc * 100.0
    else:
        row["zscore_retained_percentage"] = None

    if zscore_proc is not None and input_snps > 0:
        row["variance_retention_percentage"] = zscore_proc / input_snps * 100.0
        row["variance_removed_percentage"] = 100.0 - row["variance_retention_percentage"]
    else:
        row["variance_retention_percentage"] = row["variance_removed_percentage"] = None

    # if row["gemma_time_milliseconds"] is not None:
    #     row["gemma_time_ms"] = row["gemma_time_milliseconds"]
    # else:
    #     row["gemma_time_ms"] = None

    gemma_ms = row["gemma_time_ms"]
    if gemma_ms is not None and row["zscore_time_ms"] and row["zscore_time_ms"] > 0:
        row["zscore_speedup_vs_gemma"] = gemma_ms / row["zscore_time_ms"]
    else:
        row["zscore_speedup_vs_gemma"] = None

    if gemma_ms is not None and row["zscore_plus_gemma_time_ms"] and row["zscore_plus_gemma_time_ms"] > 0:
        row["total_pipeline_speedup_vs_gemma"] = gemma_ms / row["zscore_plus_gemma_time_ms"]
    else:
        row["total_pipeline_speedup_vs_gemma"] = None

    gemma_tested = row["gemma_tested_snps"]
    if gemma_tested and gemma_tested > 0 and zscore_rel is not None:
        row["gemma_candidate_retention_percentage"] = zscore_rel / gemma_tested * 100.0
        row["gemma_candidate_reduction_percentage"] = 100.0 - row["gemma_candidate_retention_percentage"]
    else:
        row["gemma_candidate_retention_percentage"] = row["gemma_candidate_reduction_percentage"] = None

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
    p.add_argument("--removal", default="10,20,30,40")
    p.add_argument("--removed-cases", default="25,50,75")
    p.add_argument("--replicates", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("-g", "--genotype-file", action="store_true")
    p.add_argument("-2", "--pheno12", action="store_true")
    args = p.parse_args()
    
    matrix = Path(args.matrix).resolve()
    pheno_file = Path(args.pheno).resolve()
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    results = Path(args.results).resolve()
    
    workdir.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)

    # Clean legacy global logs if they exist to prevent appending to old data
    for log_name in ["relevant_log.csv", "time_log.csv", "fn_log.csv"]:
        log_file = results / log_name
        if log_file.exists():
            log_file.unlink()

    pheno = read_pheno(pheno_file)
    rows = read_matrix(matrix, args.genotype_file)
    if any(len(values) != len(pheno) for _, values in rows):
        raise ValueError("Matrix and phenotype file have different sample counts.")
        
    y = labels(pheno, args.pheno12)
    cases = [i for i, value in enumerate(y) if value == 1]
    controls = [i for i, value in enumerate(y) if value == 0]
    n = len(pheno)
    input_snps = len(rows)
    
    print(f"Samples: {n} | Cases: {len(cases)} | Controls: {len(controls)} | SNPs: {input_snps}")
    
    all_experiments = []

    with tempfile.TemporaryDirectory(dir=workdir, prefix="subset_tmp_") as tmp:
        tmpdir = Path(tmp)

        # =====================================================
        # FULL DATASET
        # =====================================================
        print("\n=== FULL DATASET ===")
        full_dir = results / "FULL"
        full_matrix = tmpdir / "full_matrix.txt"
        full_pheno = tmpdir / "full_pheno.txt"
        write_subset(rows, pheno, list(range(n)), full_matrix, full_pheno)

        rc, elapsed = run_c(binary, full_dir, full_matrix, full_pheno, args.genotype_file, args.pheno12, 1, cwd=results)
        if rc != 0:
            raise RuntimeError(f"Full run failed with exit code {rc}")
            
        specs_path = full_dir / "specs.txt"
        if not specs_path.exists():
            raise RuntimeError("C program did not produce FULL/specs.txt")

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
            "return_code": rc
        })
        print(f"Full run completed in {elapsed:.2f}s")

        # =====================================================
        # REMOVAL EXPERIMENTS
        # =====================================================
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
                    
                    removed = set(rng.sample(cases, n_remove_cases) + rng.sample(controls, n_remove_controls))
                    keep = [i for i in range(n) if i not in removed]

                    matrix_out = tmpdir / f"{name}_matrix.txt"
                    pheno_out = tmpdir / f"{name}_pheno.txt"
                    write_subset(rows, pheno, keep, matrix_out, pheno_out)

                    run_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(specs_path, run_dir / "specs.txt")

                    rc, elapsed = run_c(binary, run_dir, matrix_out, pheno_out, args.genotype_file, args.pheno12, 1, cwd=results)

                    # Clean up temporary dataset files immediately after C finishes reading them
                    matrix_out.unlink(missing_ok=True)
                    pheno_out.unlink(missing_ok=True)
                    
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
                        "return_code": rc
                    })
                    print(f"  {name} completed in {elapsed:.2f}s (code {rc})")

    # =========================================================
    # PARSE GLOBAL LOGS AND COMBINE METRICS
    # =========================================================
    print("\nMerging data and generating reports...")
    parsed_logs = parse_global_logs(results)
    final_results = []
    
    for exp in all_experiments:
        folder = exp.pop("folder")
        log_data = parsed_logs.get(folder, {})
        final_results.append(compute_metrics(exp, log_data, exp["input_snps"]))

    if not final_results:
        print("No results generated.")
        return

    # =========================================================
    # EXPORT RESULTS
    # =========================================================
    preferred_columns = [
        "experiment", "removal_pct", "removed_case_pct", "replicate", "total_samples",
        "removed_samples", "removed_cases", "removed_controls", "remaining_samples",
        "input_snps", "zscore_processed_snps", "zscore_relevant_snps", "zscore_percentage",
        "zscore_retained_percentage", "variance_retention_percentage", "variance_removed_percentage",
        "gemma_tested_snps", "gemma_relevant_snps", "gemma_percentage",
        "gemma_after_zscore_tested_snps", "gemma_after_zscore_relevant_snps",
        "gemma_candidate_retention_percentage", "gemma_candidate_reduction_percentage",
        "false_negatives_total", "false_negative_rate", "recall",
        "gemma_time_milliseconds", "gemma_time_ms", "zscore_time_ms", "gemma_after_zscore_time_ms",
        "zscore_plus_gemma_time_ms", "zscore_speedup_vs_gemma", "total_pipeline_speedup_vs_gemma",
        "python_wall_time_seconds", "return_code"
    ]
    
    fieldnames = [col for col in preferred_columns if any(col in row for row in final_results)]

    report_csv = results / "report.csv"
    with open(report_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_results)

    report_tsv = results / "report.tsv"
    with open(report_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(final_results)

    print("========================================")
    print("Done.")
    print(f"CSV report: {report_csv}")
    print(f"TSV report: {report_tsv}")
    print("========================================")

if __name__ == "__main__":
    main()