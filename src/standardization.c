#include "suspect.h"

/* SNP validation, standardization, and SNP-row sampling/extraction. */

int *collect_valid_snp_indices_for_GRM(const int *mask, int snp_size, int *n_valid_snps)
{
    int count = 0;
    for (int j = 0; j < snp_size; ++j) {
        if (mask[j]) count++;
    }

    if (count <= 0) {
        fprintf(stderr, "No valid SNPs remain after variance filtering.\n");
        exit(EXIT_FAILURE);
    }

    int *indices = (int *)malloc((size_t)count * sizeof(int));
    if (indices == NULL) {
        perror("malloc valid_snp_indices");
        exit(EXIT_FAILURE);
    }

    int k = 0;
    for (int j = 0; j < snp_size; ++j) {
        if (mask[j]) indices[k++] = j;
    }

    *n_valid_snps = count;
    return indices;
}

int *sample_snp_row_indices_from_valid(const int *valid_snp_indices, int n_valid_snps, double sample_fraction, int *sampled_snp_size, unsigned int seed)
{
    *sampled_snp_size = (int)floor(sample_fraction * (double)n_valid_snps);
    if (*sampled_snp_size < 1) *sampled_snp_size = 1;
    if (*sampled_snp_size > n_valid_snps) *sampled_snp_size = n_valid_snps;

    int *indices = (int *)malloc((size_t)n_valid_snps * sizeof(int));
    if (indices == NULL) 
    {
        perror("malloc sample indices");
        exit(EXIT_FAILURE);
    }

    memcpy(indices, valid_snp_indices, (size_t)n_valid_snps * sizeof(int));

    srand(seed);
    for (int i = n_valid_snps - 1; i > 0; --i) 
    {
        int j = rand() % (i + 1);
        int tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }

    int *sampled = (int *)malloc((size_t)(*sampled_snp_size) * sizeof(int));
    if (sampled == NULL) 
    {
        perror("malloc sampled indices");
        free(indices);
        exit(EXIT_FAILURE);
    }

    memcpy(sampled, indices, (size_t)(*sampled_snp_size) * sizeof(int));
    free(indices);
    return sampled;
}

void process_snp_row(double *row, int original_sample_size, int *skipped_samples, double *mean_incl_skipped, double *std_incl_skipped, double *mean_no_skipped, double *std_no_skipped)
{
    double sum_incl = 0.0;
    double sum_sq_incl = 0.0;
    int count_incl = 0;
    int n_miss_incl = 0;

    double sum_no = 0.0;
    double sum_sq_no = 0.0;
    int count_no = 0;
    int n_miss_no = 0;

    for (int i = 0; i < original_sample_size; ++i) {
        double val = row[i];

        if (val < 0.0) {
            n_miss_incl++;
            if (skipped_samples[i]) {
                n_miss_no++;
            }
            continue;
        }

        sum_incl += val;
        sum_sq_incl += val * val;
        count_incl++;

        if (skipped_samples[i]) {
            sum_no += val;
            sum_sq_no += val * val;
            count_no++;
        }
    }

    if (count_incl > 0) {
        *mean_incl_skipped = sum_incl / count_incl;

        double centered_ss_incl = sum_sq_incl - (sum_incl * sum_incl) / count_incl;
        int total_incl = count_incl + n_miss_incl;
        *std_incl_skipped = (total_incl > 0 && centered_ss_incl > 0.0)
            ? sqrt(centered_ss_incl / (double)total_incl)
            : 0.0;
    } else {
        *mean_incl_skipped = 0.0;
        *std_incl_skipped = 0.0;
    }

    if (count_no > 0) {
        *mean_no_skipped = sum_no / count_no;

        double centered_ss_no = sum_sq_no - (sum_no * sum_no) / count_no;
        int total_no = count_no + n_miss_no;
        *std_no_skipped = (total_no > 0 && centered_ss_no > 0.0)
            ? sqrt(centered_ss_no / (double)total_no)
            : 0.0;
    } else {
        *mean_no_skipped = 0.0;
        *std_no_skipped = 0.0;
    }
}





void standardize_snp_matrix_row(double *row, int original_sample_size, double mean, double std, double std_no_skipped, int *keep_mask, double *std_row)
{
    /* std_row == NULL: only the validity mask is computed (subset-reuse path, kinship comes from the cache) */
    if (std_row != NULL)
        for (int i = 0; i < original_sample_size; ++i)
            std_row[i] = 0.0;

    if (std_no_skipped <= 0.0 || std <= 0.0) {
        *keep_mask = 0;
        return;
    }

    *keep_mask = 1;
    if (std_row == NULL) return;
    for (int i = 0; i < original_sample_size; ++i) {
        if (row[i] >= 0.0) std_row[i] = (row[i] - mean) / std;
    }
}


void standardize_snp_matrix(folder_geno_pheno *data, double *std_matrix, int *n_valid_snps)
{
    int valid_snp_count = 0;

    data->snp_means_including_skipped_samples = (double *)calloc(data->snp_size, sizeof(double));
    data->snp_stds_including_skipped_samples = (double *)calloc(data->snp_size, sizeof(double));
    data->snp_means_no_skipped_samples = (double *)calloc(data->snp_size, sizeof(double));
    data->snp_stds_no_skipped_samples = (double *)calloc(data->snp_size, sizeof(double));

    if(!data->snp_means_including_skipped_samples || !data->snp_stds_including_skipped_samples ||
       !data->snp_means_no_skipped_samples || !data->snp_stds_no_skipped_samples) 
    {
        perror("Memory allocation for SNP means and stds failed");
        exit(EXIT_FAILURE);
    }
 
    #pragma omp parallel for reduction(+:valid_snp_count)
    for (int j = 0; j < data->snp_size; ++j) 
    {
        double *row = data->full_geno_matrix + (size_t)j * (size_t)data->original_sample_size;
 
        
        process_snp_row(row, data->original_sample_size, 
                        data->skipped_samples,
                        &data->snp_means_including_skipped_samples[j], 
                        &data->snp_stds_including_skipped_samples[j], 
                        &data->snp_means_no_skipped_samples[j], 
                        &data->snp_stds_no_skipped_samples[j]);

        standardize_snp_matrix_row(row, data->original_sample_size, 
                        data->snp_means_including_skipped_samples[j], 
                        data->snp_stds_including_skipped_samples[j], 
                        data->snp_stds_no_skipped_samples[j],   
                        &data->snp_keep_mask[j],
                        (std_matrix != NULL) ? std_matrix + (size_t)j * (size_t)data->original_sample_size : NULL);
        
        if(data->snp_keep_mask[j])
        {
            valid_snp_count++;
        }


    }    
 
    *n_valid_snps = valid_snp_count;
}



int *sample_snp_row_indices(int snp_size, double sample_fraction, int *sampled_snp_size, unsigned int seed)
{
    *sampled_snp_size = (int)(sample_fraction * snp_size);

    int *indices = malloc(snp_size * sizeof(int));
    for (int i = 0; i < snp_size; i++)
        indices[i] = i;

    srand(seed);

    /* Fisher-Yates shuffle */
    for (int i = snp_size - 1; i > 0; i--) 
    {
        int j = rand() % (i + 1);
        int tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }

    int *sampled = malloc((*sampled_snp_size) * sizeof(int));
    memcpy(sampled, indices, (*sampled_snp_size) * sizeof(int));

    free(indices);

    return sampled;
}



double *extract_subset_std_matrix(const double *std_matrix, int original_sample_size, const int *sampled_snp_indices, int sampled_snp_size)
{
    double *subset_std_matrix = (double *)malloc((size_t)sampled_snp_size * (size_t)original_sample_size * sizeof(double));
    if (subset_std_matrix == NULL) 
    {
        perror("malloc subset_std_matrix");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sampled_snp_size; i++) 
    {
        memcpy(subset_std_matrix + (size_t)i * (size_t)original_sample_size, std_matrix + (size_t)sampled_snp_indices[i] * (size_t)original_sample_size, (size_t)original_sample_size * sizeof(double));
    }

    return subset_std_matrix;
}


