#include "raft.h"

/* phenotype transforms, IPW surprise weights, and propensity cache. */

double mean_vector( int *y, int pruned_sample_size)
{
    double s = 0.0;
    for (int i = 0; i < pruned_sample_size; i++) {
        s += y[i];
    }
    return s / (double)pruned_sample_size;
}

int count_cases(int *y, int pruned_sample_size)
{
    int count = 0;
    for (int i = 0; i < pruned_sample_size; i++) {
        if (y[i] == 1) {
            count++;
        }
    }
    return count;
}



double *spectral_transform_phenotype(int *pheno_vector, const double *U, const double *eigvals, int pruned_sample_size, double delta)
{
    double mean_y = mean_vector(pheno_vector, pruned_sample_size);

    double *y_centered = malloc(pruned_sample_size * sizeof(double));
    double *tmp       = malloc(pruned_sample_size * sizeof(double));
    double *tmp2      = malloc(pruned_sample_size * sizeof(double));
    double *p       = malloc(pruned_sample_size * sizeof(double));

    if (!y_centered || !tmp || !tmp2 || !p) {
        free(y_centered);
        free(tmp);
        free(tmp2);
        free(p);
        return NULL;
    }

    /* y_centered = y - mean_y */
    for (int i = 0; i < pruned_sample_size; i++) {
        y_centered[i] = pheno_vector[i] - mean_y;
    }

    /* tmp = U^T * y_centered   (row-major U: contiguous GEMV instead of a column-strided loop) */
    cblas_dgemv(CblasRowMajor, CblasTrans, pruned_sample_size, pruned_sample_size,
                1.0, U, pruned_sample_size, y_centered, 1, 0.0, tmp, 1);

    /* tmp2[i] = D[i] * tmp[i] */
    for (int i = 0; i < pruned_sample_size; i++) {
        double d = eigvals[i] / (eigvals[i] + delta);
        tmp2[i] = d * tmp[i];
    }

    /* out = U * tmp2 */
    cblas_dgemv(CblasRowMajor, CblasNoTrans, pruned_sample_size, pruned_sample_size,
                1.0, U, pruned_sample_size, tmp2, 1, 0.0, p, 1);

    for (int i = 0; i < pruned_sample_size; i++) {
        p[i] = p[i] + mean_y;
        if(p[i] < 0.01) p[i] = 0.01;
        if(p[i] > 0.99) p[i] = 0.99;
    }

    free(y_centered);
    free(tmp);
    free(tmp2);

    return p;
}

double *clip_surprise_scores(double *surprise_scores, int pruned_sample_size, int *pheno_vector, double sum_weights_cases, double sum_weights_ctrls)
{
    double n_cases = (double)count_cases(pheno_vector, pruned_sample_size);
    double n_ctrls = (double)pruned_sample_size - n_cases;

    double correction_cases = n_cases / sum_weights_cases;
    double correction_ctrls = n_ctrls / sum_weights_ctrls;

    for(int i= 0; i<pruned_sample_size; i++)
    {
        surprise_scores[i] = (pheno_vector[i] == 1) ? (surprise_scores[i] * correction_cases) : (surprise_scores[i] * correction_ctrls);
    }

    return surprise_scores;
}

double *compute_clipped_IPW_surprise_score(double *spectral_pheno, int *pheno_vector, int pruned_sample_size)
{
    double *surprise_scores = malloc(pruned_sample_size * sizeof(double));
    if (!surprise_scores) {
        perror("Memory allocation for surprise_scores failed");
        return NULL;
    }

    double sum_weights_cases = 0.0;
    double sum_weights_ctrls = 0.0;

    for(int i = 0; i< pruned_sample_size; i++)
    {
        double p_i = spectral_pheno[i];
        int y_i = pheno_vector[i];

        if(y_i != 0 && y_i != 1)
        {
            fprintf(stderr, "Error: phenotype value at index %d is not binary (0 or 1): %d\n", i, y_i);
            free(surprise_scores);
            return NULL;
        }

        if(y_i == 1)
        {
            surprise_scores[i] = 1.0 / p_i;
            sum_weights_cases += surprise_scores[i];
        }
        if(y_i == 0)
        {
            surprise_scores[i] = 1.0 / (1.0 - p_i);
            sum_weights_ctrls += surprise_scores[i];
        }
    }

    double *normalized_surprises = clip_surprise_scores(surprise_scores, pruned_sample_size, pheno_vector, sum_weights_cases, sum_weights_ctrls);

    return normalized_surprises; 
}




void save_propensity_cache(const char *folder_path, const double *p, const int *keep_mask, size_t original_sample_size)
{
    double *p_full = (double *)malloc(original_sample_size * sizeof(double));
    if (p_full == NULL) { perror("malloc propensity cache"); return; }

    size_t k = 0;
    for (size_t i = 0; i < original_sample_size; i++)
        p_full[i] = keep_mask[i] ? p[k++] : NAN;

    char path[8192];
    propensity_cache_path(folder_path, path, sizeof(path));
    write_binary_atomic(path, p_full, original_sample_size);
    free(p_full);
}

/* Reference-run propensities restricted to the kept samples. NULL if the cache is unusable. */
double *load_propensity_subset(const char *folder_path, size_t original_sample_size, const int *keep_mask, size_t pruned_n)
{
    char path[8192];
    propensity_cache_path(folder_path, path, sizeof(path));

    struct stat st;
    if (stat(path, &st) != 0 || (size_t)st.st_size != original_sample_size * sizeof(double))
        return NULL;

    FILE *f = fopen(path, "rb");
    if (f == NULL) return NULL;

    double *p_full = (double *)malloc(original_sample_size * sizeof(double));
    double *p = (double *)malloc(pruned_n * sizeof(double));
    if (p_full == NULL || p == NULL || fread(p_full, sizeof(double), original_sample_size, f) != original_sample_size) {
        fclose(f); free(p_full); free(p);
        return NULL;
    }
    fclose(f);

    size_t k = 0;
    for (size_t i = 0; i < original_sample_size; i++) {
        if (!keep_mask[i]) continue;
        if (isnan(p_full[i])) {                     // sample had no phenotype in the reference run
            fprintf(stderr, "--reuse: kept sample %zu has no cached propensity; computing it from the kinship instead.\n", i);
            free(p_full); free(p);
            return NULL;
        }
        p[k++] = p_full[i];
    }

    free(p_full);
    return p;
}


double *spectral_transform_phenotype_cholesky(int *pheno_vector, double *analyzed_GRM, int pruned_sample_size, double delta)
{
    size_t n = (size_t)pruned_sample_size;
    double mean_y = mean_vector(pheno_vector, pruned_sample_size);

    double *p = (double *)malloc(n * sizeof(double));
    if (p == NULL) return NULL;

    for (size_t i = 0; i < n; i++) {
        analyzed_GRM[i * n + i] += delta;
        p[i] = pheno_vector[i] - mean_y;            // right-hand side, overwritten by the solution
    }

    lapack_int info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', (lapack_int)n, 1, analyzed_GRM, (lapack_int)n, p, (lapack_int)n);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dposv failed, info=%d\n", (int)info);
        free(p);
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        p[i] = (double)pheno_vector[i] - delta * p[i];
        if(p[i] < 0.01) p[i] = 0.01;
        if(p[i] > 0.99) p[i] = 0.99;
    }
    return p;
}
