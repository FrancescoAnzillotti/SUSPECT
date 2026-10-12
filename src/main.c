#include "raft.h"

/* Program entry point and end-to-end workflow. */

int main(int argc, char *argv[]) {
        
    FILE *relevant_log_file = fopen("relevant_log.csv", "a");
    FILE *fn_log_file = fopen("fn_log.csv", "a");
    FILE *time_log_file = fopen("time_log.csv", "a");

    if (!relevant_log_file || !fn_log_file || !time_log_file) 
    {
        perror("opening log files");
        exit(1);
    }
    
    const char *target_dir = NULL;
    const char *matrix_path = NULL;
    const char *outcome_path = NULL;
    int is_genotype_file = 0;
    int exact = -1;
    int pheno_offset = 0; // default: outcome file already 0/1, no shift
    int reuse_weights = 0; // --reuse: subset run takes the reference-run propensities instead of recomputing them from K[S,S]
    int run_baseline = 1;   // --no-baseline skips the full GEMMA run


    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            matrix_path = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outcome_path = argv[++i];
        } else if (strcmp(argv[i], "-g") == 0) {
            is_genotype_file = 1;
        } 
        else if (strcmp(argv[i], "--reuse") == 0) {
            reuse_weights = 1;
        }
        else if (strcmp(argv[i], "--no-baseline") == 0) {
            run_baseline = 0;
        }
        else if (strcmp(argv[i], "-2") == 0) {
            pheno_offset = 1; // outcome file is 1/2-coded, subtract 1 to normalize
        } else if (target_dir == NULL) {
            target_dir = argv[i];
        } else if (exact == -1) {
            exact = atoi(argv[i]);
        }
    }

    if (target_dir == NULL || exact == -1 || matrix_path == NULL || outcome_path == NULL) {
        fprintf(stderr, "Usage: %s <folder_path> <0=approximate|1=exact> -m <matrix.txt> -o <outcome.txt> [-g] [-2] [--reuse]\n", argv[0]);
        return 1;
    }

    printf("Resolving paths for matrix and outcome files...\n");
    char matrix_abs[PATH_MAX];
    char outcome_abs[PATH_MAX];
    resolve_input_path(target_dir, matrix_path, matrix_abs, sizeof(matrix_abs));
    resolve_input_path(target_dir, outcome_path, outcome_abs, sizeof(outcome_abs));


    node_z_score *z_score_list = NULL;
    folder_geno_pheno *data = NULL;

    if (!is_directory(target_dir)) {
        fprintf(stderr, "Error: '%s' is not a valid directory.\n", target_dir);
        return 1;
    }

    printf("Reading files...\n");
    data = read_files((char *)target_dir, matrix_abs, outcome_abs, is_genotype_file, pheno_offset, &z_score_list);
    if (data == NULL) {
        fprintf(stderr, "Error: Failed to read files in %s\n", target_dir);
        return 1;
    }

    double gemma_thr = gemma_threshold_for_n(data->snp_size);
    printf("SNPs: %d -> Z thr %.2f, GEMMA thr %.2f\n", data->snp_size, z_threshold_for_n(data->snp_size), gemma_thr);

    double t0 = now_seconds();


    // ------------------------------------------------ FULL AND EXACT GEMMA ---------------------------------------------

    char **gemma_relevant = NULL;
    if(run_baseline)
    {
        // Construct gemma commands
        char gemma_cmd1[16384];
        char gemma_cmd2[16384];


        snprintf(gemma_cmd1, sizeof(gemma_cmd1), 
            " cd %s && ../../gemma -g %s -p %s -gk 2 -maf 0 -miss 1 -o kinship > /dev/null 2>&1",
            target_dir, matrix_abs, outcome_abs);
        
        snprintf(gemma_cmd2, sizeof(gemma_cmd2), 
            "cd %s && ../../gemma -g %s -p %s -k output/kinship.sXX.txt -maf 0 -miss 1 -lmm 4 -n 1 -o gemma_output > /dev/null 2>&1",
            target_dir, matrix_abs, outcome_abs);

        printf("Executing GEMMA kinship...\n");
        int result = system(gemma_cmd1);
        if(result == -1) {
            perror("Failed to execute gemma kinship command");
            exit(EXIT_FAILURE);
        }
        
        printf("Executing GEMMA association testing...\n");
        result = system(gemma_cmd2);
        if(result == -1) {
            perror("Failed to execute gemma command");
            exit(EXIT_FAILURE);
        }

        double t_gemma_end = now_seconds();
        char **gemma_relevant = read_gemma_relevant_snps(target_dir, gemma_thr, relevant_log_file, 0);
        printf("Time taken for full and exact GEMMA execution: %.1f milliseconds\n\n", (t_gemma_end - t0));
        fprintf(time_log_file, "%s\tGEMMA\t%.1f\tms\n", target_dir, (t_gemma_end - t0));
    }
    

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// ------------------------------------------------ EIGEN-SHRUNK Z-SCORE APPROACH ---------------------------------------------------------------------------------------------

    double t_zscore_start = now_seconds();



    const double delta = 1.0;  

    // specs.txt + a matching kinship cache -> this is a subset run on the reference dataset
    int use_saved_K = kinship_cache_is_valid(target_dir, (size_t)data->original_sample_size, data->snp_size);

    double *spectral_pheno = NULL;
    const char *weights_source = "Cholesky solve on K[S,S]";
    if (reuse_weights && use_saved_K)
    {
        spectral_pheno = load_propensity_subset(target_dir, (size_t)data->original_sample_size, data->skipped_samples, (size_t)data->pruned_sample_size);
        if (spectral_pheno != NULL) weights_source = "reference-run propensities (--reuse)";
    }

    // allocate keep mask for candidate pruning 
    data->snp_keep_mask = (int *)calloc((size_t)data->snp_size, sizeof(int));
    if (data->snp_keep_mask == NULL) 
    {
        perror("calloc snp_keep_mask");
        exit(EXIT_FAILURE);
    }

    int n_valid_snps = 0;
    standardize_snp_matrix(data, NULL, &n_valid_snps);   // stats + validity mask only

    data->n_valid_snps = n_valid_snps;
    data->valid_snp_indices = collect_valid_snp_indices_for_GRM(data->snp_keep_mask, data->snp_size, &data->n_valid_snps);


    double *GRM = NULL;
    if (!use_saved_K)
    {
        GRM = compute_GRM_blocked(data->full_geno_matrix, data->original_sample_size, data->snp_means_including_skipped_samples, data->snp_stds_including_skipped_samples, 1.0, data->valid_snp_indices, data->n_valid_snps);
        if (GRM == NULL)
        {
            fprintf(stderr, "compute_GRM failed\n");
            exit(1);
        }
        write_kinship_binary(target_dir, GRM, (size_t)data->original_sample_size);
    }

    // K[S,S]: needed for the eigendecomposition handed to GEMMA (exact mode) and whenever the propensities
    // are not taken from the cache
    double *analyzed_GRM = NULL;
    if (exact || spectral_pheno == NULL)
    {
        if (use_saved_K)
            analyzed_GRM = read_kinship_submatrix(target_dir, (size_t)data->original_sample_size, data->skipped_samples, (size_t)data->pruned_sample_size);
        else if (data->pruned_sample_size == data->original_sample_size)
        {
            analyzed_GRM = GRM;   
            GRM = NULL;
        }
        else
            analyzed_GRM = extract_analyzed_submatrix(GRM, data->original_sample_size, data->skipped_samples, data->pruned_sample_size);

        if (analyzed_GRM == NULL) 
        {
            fprintf(stderr, "could not build the analyzed kinship submatrix\n");
            exit(1);
        }
    }

    free(GRM);
    GRM = NULL;

    double *eigvals = NULL;
    if (exact)
    {
        // GEMMA needs U and Lambda of K[S,S]. (U[S,:] is NOT the eigenbasis of K[S,S]: a new decomposition is required)
        eigvals = (double *)malloc((size_t)data->pruned_sample_size * sizeof(double));
        if (eigvals == NULL) 
        {
            perror("malloc eigvals");
            exit(1);
        }

        int info = compute_eigendecomposition(analyzed_GRM, data->pruned_sample_size, eigvals, exact, data);
        if (info != 0) 
        {
            fprintf(stderr, "LAPACKE_dsyevd failed, info=%d\n", info);
            exit(1);
        }

        if (spectral_pheno == NULL) 
        {
            weights_source = "eigendecomposition of K[S,S]";
            spectral_pheno = spectral_transform_phenotype(data->pheno_vector, analyzed_GRM, eigvals, data->pruned_sample_size, delta);
        }
    }
    else if (spectral_pheno == NULL)
    {
        // approximate mode: same operator, no eigendecomposition
        spectral_pheno = spectral_transform_phenotype_cholesky(data->pheno_vector, analyzed_GRM, data->pruned_sample_size, delta);
    }

    if (spectral_pheno == NULL) 
    {
        fprintf(stderr, "propensity computation failed\n");
        exit(1);
    }

    if (!use_saved_K)
        save_propensity_cache(target_dir, spectral_pheno, data->skipped_samples, (size_t)data->original_sample_size);

    printf("Kinship: %s | weights: %s | analyzed samples: %d of %d\n",
           use_saved_K ? "loaded from cache and subset" : "computed (reference run) and cached",
           weights_source, data->pruned_sample_size, data->original_sample_size);

    double *surprise_scores = compute_clipped_IPW_surprise_score(spectral_pheno, data->pheno_vector, data->pruned_sample_size);

    compute_z_scores(data, surprise_scores, &z_score_list, relevant_log_file, exact);
    write_to_file(target_dir, z_score_list, data->snp_keep_mask, data->snp_size);


    double t_zscore_end = now_seconds();

    fprintf(time_log_file, "%s\tZ_SCORE\t%.1f\tms\n", target_dir, (t_zscore_end - t_zscore_start));
    printf("Time taken for EIGEN-SHRUNK Z-SCORE APPROACH on folder %s: %.1f milliseconds\n", target_dir, (t_zscore_end - t_zscore_start));

    if(gemma_relevant) 
        match_z_scores_with_gemma(z_score_list, gemma_relevant, fn_log_file, data);

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




// ------------------------------------------------ EXACT GEMMA AFTER Z-SCORE APPROACH ---------------------------------------------------------------------------------------------

    char **after_z_gemma_relevant = NULL;

    if(exact)
    {
        char exact_after_zscore_cmd[16384];
        snprintf(exact_after_zscore_cmd, sizeof(exact_after_zscore_cmd), "cd %s/exact && ../../../gemma -g ./C_pruned_geno_matrix.txt -p ./C_pruned_pheno.txt -u eigenvectors.txt -d eigenvalues.txt -maf 0 -miss 1 -lmm 4 -n 1 -o gemma_output > /dev/null 2>&1", target_dir);
        int result = system(exact_after_zscore_cmd);
        if(result == -1) {
            perror("Failed to execute GEMMA exact command after z-score computation");
            exit(EXIT_FAILURE); 
        }

        double t_exact_end = now_seconds();
        fprintf(time_log_file, "%s\tGEMMA_EXACT_AFTER_ZSCORE\t%.1f\tms\n", target_dir, (t_exact_end - t_zscore_end));
        printf("Time taken for exact GEMMA execution after EIGEN-SHRUNK Z-SCORE APPROACH on folder %s: %.1f milliseconds\n", target_dir, (t_exact_end - t_zscore_end));
        printf("Total time taken for EIGEN-SHRUNK Z-SCORE APPROACH and exact GEMMA execution on remaining positions %s: %.1f milliseconds \n", target_dir, (t_exact_end - t_zscore_start));
        fprintf(time_log_file, "%s\tZSCORE+GEMMA\t%.1f\tms\n\n\n\n", target_dir, (t_exact_end - t_zscore_start));


        after_z_gemma_relevant = read_gemma_relevant_snps(target_dir, gemma_thr, relevant_log_file, 1);
        if(gemma_relevant) 
            match_gemmas(after_z_gemma_relevant, gemma_relevant, fn_log_file, data);

        free(after_z_gemma_relevant);
    }
    
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    if (data != NULL) 
    {
        free(data->pheno_vector);
        free(data->folder_name);
        free(data->full_geno_matrix); 
        free(data->skipped_samples);
        free(data->snp_keep_mask);
        free(data->valid_snp_indices);
        free(GRM);
        free(analyzed_GRM);   
        free(eigvals);
        free(spectral_pheno);
        free(surprise_scores);
        free(data);
    }

    if(z_score_list != NULL)
    {
        free_z_score_list(z_score_list);
    }    


    fflush(time_log_file);

    free(gemma_relevant);


    fclose(relevant_log_file);
    fclose(fn_log_file);
    fclose(time_log_file);

    return 0;
}