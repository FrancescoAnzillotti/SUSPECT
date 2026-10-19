#include "suspect.h"
#define GRM_BLOCK 1024

/* GRM construction, analyzed submatrices, eigendecomposition, and kinship cache. */


double *compute_GRM_blocked(const double *full_geno, int original_sample_size, const double *means, const double *stds, double sample_fraction, const int *valid_snp_indices, int n_valid_snps)
{
    size_t n = (size_t)original_sample_size;

    /* SNPs that enter the kinship: all valid ones, or a random sample of them */
    const int *used_indices = valid_snp_indices;
    int n_used = n_valid_snps;
    int *sampled = NULL;

    if (sample_fraction < 1.0)
    {
        sampled = sample_snp_row_indices_from_valid(valid_snp_indices, n_valid_snps, sample_fraction, &n_used, 42);
        if (sampled == NULL || n_used <= 0) { free(sampled); return NULL; }
        used_indices = sampled;
    }

    double *GRM = (double *)calloc(n * n, sizeof(double));
    if (GRM == NULL) { free(sampled); return NULL; }

    size_t block_rows = GRM_BLOCK;
    if (block_rows > (size_t)n_used) block_rows = (size_t)n_used;

    double *block = (double *)malloc(block_rows * n * sizeof(double));
    if (block == NULL) { free(GRM); free(sampled); return NULL; }

    for (size_t start = 0; start < (size_t)n_used; start += block_rows)
    {
        size_t rows = (size_t)n_used - start;
        if (rows > block_rows) rows = block_rows;

        #pragma omp parallel for schedule(static)
        for (size_t r = 0; r < rows; r++)
        {
            int j = used_indices[start + r];
            const double *src = full_geno + (size_t)j * n;
            double *dst = block + r * n;
            double mu = means[j], sd = stds[j];
            for (size_t i = 0; i < n; i++)
                dst[i] = (src[i] >= 0.0) ? (src[i] - mu) / sd : 0.0;   // missing -> 0 
        }

        // beta = 1.0: accumulate this block into the upper triangle 
        cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans,
                    (int)n, (int)rows, 1.0, block, (int)n, 1.0, GRM, (int)n);
    }
    free(block);
    free(sampled);

    double inv = 1.0 / (double)n_used;
    #pragma omp parallel for schedule(dynamic, 16)
    for (size_t i = 0; i < n; i++)
        for (size_t j = i; j < n; j++) {
            double v = GRM[i * n + j] * inv;
            GRM[i * n + j] = v;
            GRM[j * n + i] = v;
        }
    return GRM;
}













/*----------------------------------------------------------
  Extract the pruned_n x pruned_n submatrix of a full_n x full_n
  matrix, keeping only rows/cols where keep_mask[i] == 1.
  This matches GEMMA's internal behavior of restricting the
  kinship matrix to phenotype-complete individuals before
  running -lmm.
----------------------------------------------------------*/
double *extract_analyzed_submatrix(const double *full_matrix, int full_n, const int *keep_mask, int pruned_n)
{
    size_t n_full = (size_t)full_n;
    size_t n_pruned = (size_t)pruned_n;

    double *sub = malloc(n_pruned * n_pruned * sizeof(double));
    if (sub == NULL) {
        perror("malloc analyzed submatrix");
        exit(1);
    }

    if (n_pruned == n_full) {                       // nothing dropped
        memcpy(sub, full_matrix, n_full * n_full * sizeof(double));
        return sub;
    }

    size_t *kept_idx = malloc(n_pruned * sizeof(size_t));
    if (kept_idx == NULL) {
        perror("malloc kept_idx");
        exit(1);
    }

    size_t k = 0;
    for (size_t i = 0; i < n_full; i++) {
        if (keep_mask[i]) kept_idx[k++] = i;
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_pruned; i++) {
        const double *src_row = full_matrix + kept_idx[i] * n_full;
        double *dst_row = sub + i * n_pruned;
        for (size_t j = 0; j < n_pruned; j++) {
            dst_row[j] = src_row[kept_idx[j]];
        }
    }

    free(kept_idx);
    return sub;
}




void write_eigen_decomposition_to_file(double *U, double *eigvals, int n, const char *folder_path)
{
    char eigenvectors_output_path[8192];
    char eigenvalues_output_path[8192];

    snprintf(eigenvectors_output_path, sizeof(eigenvectors_output_path),
             "%s/exact/eigenvectors.txt", folder_path);
    snprintf(eigenvalues_output_path, sizeof(eigenvalues_output_path),
             "%s/exact/eigenvalues.txt", folder_path);

    FILE *eigenvectors_output_file = fopen(eigenvectors_output_path, "w");
    if (!eigenvectors_output_file) {
        perror(eigenvectors_output_path);
        return;
    }

    FILE *eigenvalues_output_file = fopen(eigenvalues_output_path, "w");
    if (!eigenvalues_output_file) {
        perror(eigenvalues_output_path);
        fclose(eigenvectors_output_file);
        return;
    }

    for (int i = 0; i < n; i++) {
        fprintf(eigenvalues_output_file, "%.17g\n", eigvals[i]);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(eigenvectors_output_file, j == 0 ? "%.17g" : "\t%.17g",
                    U[(size_t)i * n + j]);
        }
        fputc('\n', eigenvectors_output_file);
    }

    fclose(eigenvectors_output_file);
    fclose(eigenvalues_output_file);
}


void reverse_eigen_decomposition(double *U, double *eigvals, int n)
{
    for(int i = 0; i < n / 2; i++) {
        double temp = eigvals[i];
        eigvals[i] = eigvals[n - 1 - i];
        eigvals[n - 1 - i] = temp;
    }
    for(size_t row = 0; row < (size_t)n; row++) {
        for(size_t col = 0; col < (size_t)n / 2; col++) {
            double temp = U[row * n + col];
            U[row * n + col] = U[row * n + (n - 1 - col)];
            U[row * n + (n - 1 - col)] = temp;
        }
    }
}





/*----------------------------------------------------------
  Compute Eigendecomposition of the GRM

      GRM = U * diag(eigvals) * U^T

  Input:
      GRM
      Shape: pruned_sample_size × pruned_sample_size
      eigvals (empty array of size pruned_sample_size)

  Output:
      grm overwritten with eigenvectors
      eigvals filled with eigenvalues
----------------------------------------------------------*/

int compute_eigendecomposition(double *grm, int pruned_sample_size, double *eigvals, int exact, folder_geno_pheno *data)
{
    int info;

    /* grm will be overwritten with eigenvectors */
    info = LAPACKE_dsyevd(
        LAPACK_ROW_MAJOR,
        'V',   /* eigenvalues + eigenvectors */
        'U',   /* use upper triangle of symmetric matrix */
        pruned_sample_size,
        grm,
        pruned_sample_size,
        eigvals
    );
    if(info == 0) 
    {
        reverse_eigen_decomposition(grm, eigvals, pruned_sample_size);
    }

    if(exact)
    {
        write_eigen_decomposition_to_file(grm, eigvals, pruned_sample_size, data->folder_name);
    }

    return info;
}





void write_kinship(const char *folder_path, double *GRM, int original_sample_size) 
{
    char kinship_output_path[8192];
    snprintf(kinship_output_path, sizeof(kinship_output_path), "%s/C_COMPUTED_kinship_matrix.txt", folder_path);

    FILE *kinship_file = fopen(kinship_output_path, "w");
    if (kinship_file == NULL) {
        perror(kinship_output_path);
        return;
    }

    for (int i = 0; i < original_sample_size; i++) {
        for (int j = 0; j < original_sample_size; j++) {
            fprintf(kinship_file, "%.8g", GRM[(size_t)i * original_sample_size + j]);
            if (j < original_sample_size - 1) {
                fprintf(kinship_file, "\t");
            }
        }
        fprintf(kinship_file, "\n");
    }

    fclose(kinship_file);
}



/*
  SUBSET REUSE
  --------------------------------------------------------------------------------------
  Reference run  (no specs.txt, or no valid cache): standardize + GRM (O(M N^2)), and the
                 full original_sample_size x original_sample_size kinship is written, raw binary doubles,
                 to C_COMPUTED_kinship_matrix.bin together with the propensity scores
                 (C_COMPUTED_propensity.bin, used only by --reuse).
  Subset run     (specs.txt + valid cache): no standardized matrix, no GRM. The kept samples are those
                 with a valid phenotype (skipped_samples, i.e. the NA pattern of the outcome file,
                 which is also what GEMMA reads). The kinship is memory-mapped and only the kept
                 rows/columns are copied out: K[S,S] is exactly the kinship of the subset.
*/

void kinship_cache_path(const char *folder_path, char *out, size_t out_size)
{
    snprintf(out, out_size, "%s/C_COMPUTED_kinship_matrix.bin", folder_path);
}

void propensity_cache_path(const char *folder_path, char *out, size_t out_size)
{
    snprintf(out, out_size, "%s/C_COMPUTED_propensity.bin", folder_path);
}

/* 1 if specs.txt exists AND the kinship cache matches this dataset (same number of SNPs and samples) */
int kinship_cache_is_valid(const char *folder_path, size_t original_sample_size, int snp_size)
{
    int ref_n_samples = 0, ref_n_candidates = 0;
    double ref_z_threshold = 0.0;

    if (!read_specs(folder_path, &ref_n_samples, &ref_n_candidates, &ref_z_threshold))
        return 0;                                   // no specs.txt: this is the reference run

    char path[8192];
    kinship_cache_path(folder_path, path, sizeof(path));

    struct stat st;
    if (stat(path, &st) != 0 ||
        (size_t)st.st_size != original_sample_size * original_sample_size * sizeof(double) ||
        ref_n_candidates != snp_size)
    {
        fprintf(stderr, "specs.txt found but kinship cache %s is missing or does not match this dataset: recomputing it.\n", path);
        return 0;
    }
    return 1;
}

int write_binary_atomic(const char *path, const double *values, size_t n_values)
{
    char tmp_path[8192 + 8];
    snprintf(tmp_path, sizeof(tmp_path), "%s.tmp", path);

    FILE *f = fopen(tmp_path, "wb");
    if (f == NULL) { perror(tmp_path); return -1; }

    if (fwrite(values, sizeof(double), n_values, f) != n_values) {
        perror(tmp_path);
        fclose(f);
        remove(tmp_path);
        return -1;
    }
    fclose(f);

    if (rename(tmp_path, path) != 0) { perror(path); return -1; }   // atomic: readers never see a half-written file
    return 0;
}

int write_kinship_binary(const char *folder_path, const double *GRM, size_t original_sample_size)
{
    char path[8192];
    kinship_cache_path(folder_path, path, sizeof(path));
    return write_binary_atomic(path, GRM, original_sample_size * original_sample_size);
}

/* Memory-maps the cached kinship and copies out K[kept, kept]. Dropped rows are never touched. */
double *read_kinship_submatrix(const char *folder_path, size_t full_n, const int *keep_mask, size_t pruned_n)
{
    char path[8192];
    kinship_cache_path(folder_path, path, sizeof(path));

    int fd = open(path, O_RDONLY);
    if (fd < 0) { perror(path); return NULL; }

    size_t bytes = full_n * full_n * sizeof(double);
    const double *K = (const double *)mmap(NULL, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (K == MAP_FAILED) { perror("mmap kinship"); return NULL; }

    double *sub = (double *)malloc(pruned_n * pruned_n * sizeof(double));
    size_t *kept_idx = (size_t *)malloc(pruned_n * sizeof(size_t));
    if (sub == NULL || kept_idx == NULL) {
        perror("malloc kinship submatrix");
        free(sub); free(kept_idx);
        munmap((void *)K, bytes);
        return NULL;
    }

    size_t k = 0;
    for (size_t i = 0; i < full_n; i++) {
        if (keep_mask[i]) kept_idx[k++] = i;
    }

    if (pruned_n == full_n) {                       // nothing dropped
        memcpy(sub, K, bytes);
    } else {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < pruned_n; i++) {
            const double *src_row = K + kept_idx[i] * full_n;
            double *dst_row = sub + i * pruned_n;
            for (size_t j = 0; j < pruned_n; j++) {
                dst_row[j] = src_row[kept_idx[j]];
            }
        }
    }

    free(kept_idx);
    munmap((void *)K, bytes);
    return sub;
}

/* Propensity scores of the reference run, scattered to original sample positions (NA = no phenotype). */
