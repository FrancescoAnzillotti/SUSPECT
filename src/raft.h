#ifndef RAFT_H
#define RAFT_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <cblas.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <lapacke.h>
#include <time.h>
#include <omp.h>
#include <limits.h>
#include <sys/mman.h>
#include <fcntl.h>

#define max_rows 2000
#define buffer_size 57344
#define pheno_buffer 128

typedef struct
{
    char *folder_name;
    int *pheno_vector;
    int pruned_sample_size;
    //int *geno_matrix;            // pruned, integer
    double *full_geno_matrix;    // full dosages
    int snp_size;

    int *snp_keep_mask;          // 1 = keep SNP, 0 = drop SNP
    int *valid_snp_indices;      // raw SNP row indices that survived filtering
    int n_valid_snps;            // number of surviving SNPs

    int *skipped_samples;        // keep mask for samples
    int original_sample_size;

    double *snp_means_including_skipped_samples;
    double *snp_stds_including_skipped_samples;
    double *snp_means_no_skipped_samples;
    double *snp_stds_no_skipped_samples;

    double adjusted_z_score_threshold; // adjusted threshold for z-score filtering

} folder_geno_pheno;

typedef struct node_z_score
{
    char *snp_name;
    float real_z_score;
    struct node_z_score *next;
} node_z_score;

/* OpenBLAS helpers retained from the original implementation. */
extern void openblas_set_num_threads(int);
extern int openblas_get_num_threads(void);

/* time.c */
double now_seconds(void);

/* utils.c */
int is_directory(const char *path);
int is_integer_and_positive(int value);
int parse_double_token(const char *token, double *out);
int *allocate_pheno_vector(void);
//int *allocate_geno_matrix(int pruned_sample_size, int n_snps);
void print_array(int *array, int length);
void free_z_score_list(node_z_score *head);
void free_mem_close_files(FILE *geno, FILE *pheno, folder_geno_pheno *data, node_z_score *head);
void create_node_and_head_insert(node_z_score **head, char *snp_name);
void print_list(node_z_score *head);
void resolve_input_path(const char *target_dir, const char *given_path, char *out_abs, size_t out_size);

/* data_io.c */
int get_n_snps_in_file(FILE *geno_file);
int *read_pheno_vector(FILE *pheno_file, int *pruned_sample_size, int **keep_mask, int *original_sample_size, int pheno_offset);
void *read_geno_matrix(FILE *geno_file, int *pruned_sample_size, int *keep_mask, int *original_sample_size, int *snp_size, node_z_score **head, double **out_full_geno, int is_genotype_file);
folder_geno_pheno *read_files(char *folder_path, const char *matrix_path, const char *outcome_path, int is_genotype_file, int pheno_offset, node_z_score **head);
void write_geno_col_to_file(FILE *output_file, folder_geno_pheno *data, int snp_index, char *snp_name);
void write_pruned_pheno_to_file(folder_geno_pheno *data);
void write_specs_if_absent(const char *folder_path, int n_samples, int n_candidates, double z_threshold);
int read_specs(const char *folder_path, int *ref_n_samples, int *ref_n_candidates, double *ref_z_threshold);

/* standardization.c */
int *collect_valid_snp_indices_for_GRM(const int *mask, int snp_size, int *n_valid_snps);
int *sample_snp_row_indices_from_valid(const int *valid_snp_indices, int n_valid_snps, double sample_fraction, int *sampled_snp_size, unsigned int seed);
void process_snp_row(double *row, int original_sample_size, int *skipped_samples, double *mean_incl_skipped, double *std_incl_skipped, double *mean_no_skipped, double *std_no_skipped);
void standardize_snp_matrix_row(double *row, int original_sample_size, double mean, double std, double std_no_skipped, int *keep_mask, double *std_row);
void standardize_snp_matrix(folder_geno_pheno *data, double *std_matrix, int *n_valid_snps);
int *sample_snp_row_indices(int snp_size, double sample_fraction, int *sampled_snp_size, unsigned int seed);
double *extract_subset_std_matrix(const double *std_matrix, int original_sample_size, const int *sampled_snp_indices, int sampled_snp_size);

/* kinship.c */
//void compute_grm_CALL(double *subset_std_matrix, int sampled_snp_size, int original_sample_size, double *grm, int m_divisor);
//double *compute_GRM(const double *std_matrix, int snp_size, int original_sample_size, double sample_fraction, const int *valid_snp_indices, int n_valid_snps);
double *compute_GRM_blocked(const double *full_geno, int original_sample_size, const double *means, const double *stds, double sample_fraction, const int *valid_snp_indices, int n_valid_snps);
double *extract_analyzed_submatrix(const double *full_matrix, int full_n, const int *keep_mask, int pruned_n);
void write_eigen_decomposition_to_file(double *U, double *eigvals, int n, const char *folder_path);
void reverse_eigen_decomposition(double *U, double *eigvals, int n);
int compute_eigendecomposition(double *grm, int pruned_sample_size, double *eigvals, int exact, folder_geno_pheno *data);
void write_kinship(const char *folder_path, double *GRM, int original_sample_size);
void kinship_cache_path(const char *folder_path, char *out, size_t out_size);
int kinship_cache_is_valid(const char *folder_path, size_t original_sample_size, int snp_size);
int write_binary_atomic(const char *path, const double *values, size_t n_values);
int write_kinship_binary(const char *folder_path, const double *GRM, size_t original_sample_size);
double *read_kinship_submatrix(const char *folder_path, size_t full_n, const int *keep_mask, size_t pruned_n);

/* propensity.c */
double mean_vector(int *y, int pruned_sample_size);
int count_cases(int *y, int pruned_sample_size);
double *spectral_transform_phenotype(int *pheno_vector, const double *U, const double *eigvals, int pruned_sample_size, double delta);
double *clip_surprise_scores(double *surprise_scores, int pruned_sample_size, int *pheno_vector, double sum_weights_cases, double sum_weights_ctrls);
double *compute_clipped_IPW_surprise_score(double *spectral_pheno, int *pheno_vector, int pruned_sample_size);
void propensity_cache_path(const char *folder_path, char *out, size_t out_size);
void save_propensity_cache(const char *folder_path, const double *p, const int *keep_mask, size_t original_sample_size);
double *load_propensity_subset(const char *folder_path, size_t original_sample_size, const int *keep_mask, size_t pruned_n);
double *spectral_transform_phenotype_cholesky(int *pheno_vector, double *analyzed_GRM, int pruned_sample_size, double delta);

/* zscore.c */
//void compute_four_counts(folder_geno_pheno *data, int snp_index, double *counts, char *snp_col, double *surprise_scores);
void compute_four_counts(folder_geno_pheno *data, int snp_index, double *counts, double *surprise_scores);
double z_threshold_for_n(int n_snps);
double gemma_threshold_for_n(int n_snps);
void compute_z_scores(folder_geno_pheno *data, double *surprise_scores, node_z_score **head, FILE *relevant_log_file, int exact);
void write_to_file(const char *folder_path, node_z_score *head, int *snp_keep_mask, int snp_size);

/* gemma.c */
char **read_gemma_relevant_snps(const char *folder_path, double thr, FILE *relevant_log_file, int exact);
void match_z_scores_with_gemma(node_z_score *z_score_list, char **gemma_relevant, FILE *fn_log_file, folder_geno_pheno *data);
void match_gemmas(char **after_z_gemma_relevant, char **gemma_relevant, FILE *fn_log_file, folder_geno_pheno *data);

#endif /* RAFT_H */
