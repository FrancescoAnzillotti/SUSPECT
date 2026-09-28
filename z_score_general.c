#define _GNU_SOURCE

#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <lapacke.h>
#include <time.h>
#include <omp.h>
#include <limits.h> 

#define max_rows 2000
//#define max_snps 100000
#define buffer_size 57344
#define pheno_buffer 128


typedef struct
{
    char *folder_name;
    int *pheno_vector;
    int pruned_sample_size;
    int *geno_matrix;            // pruned, integer
    double *full_geno_matrix;    // full dosages
    int snp_size;

    int *snp_keep_mask; // 1 = keep SNP, 0 = drop SNP
    int *valid_snp_indices;       // raw SNP row indices that survived filtering
    int n_valid_snps;             // number of surviving SNPs

    int *skipped_samples;         // keep mask for samples
    int original_sample_size;

    double *snp_means_including_skipped_samples;
    double *snp_stds_including_skipped_samples;
    double *snp_means_no_skipped_samples;
    double *snp_stds_no_skipped_samples;

} folder_geno_pheno;

typedef struct node_z_score
{
    char *snp_name;
    float real_z_score;
    float linear_z_score;
    float corrected_linear_z_score;
    struct node_z_score *next;
} node_z_score;



extern void openblas_set_num_threads(int);
extern int openblas_get_num_threads(void);

double now_seconds()  //millisecons
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    // Convert seconds to ms, and nanoseconds to ms, then add them
    return (ts.tv_sec * 1000.0) + (ts.tv_nsec / 1e6);
}

int is_directory(const char *path) 
{
    struct stat st;
    return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

int is_integer_and_positive(int value)
{
    return value>=0;
}

int parse_double_token(const char *token, double *out)
{
    char *end = NULL;
    errno = 0;

    double x = strtod(token, &end);
    if (end == token) return 0;              // no conversion
    while (*end && isspace((unsigned char)*end)) end++;
    if (*end != '\0') return 0;              // extra junk after number
    if (errno == ERANGE) return 0;

    *out = x;
    return 1;
}


int * allocate_pheno_vector(void)
{
    int *pheno = (int *)malloc(max_rows*sizeof(int));

    if (pheno == NULL) 
    {
        perror("Initial phenotype memory allocation failed");
        exit(1);
    }

    return pheno;
}

int *allocate_geno_matrix(int pruned_sample_size, int n_snps)
{
    int *geno = (int *)malloc(pruned_sample_size * n_snps * sizeof(int));

    if (geno == NULL) 
    {
        perror("Initial genotype memory allocation failed");
        exit(1);
    }

    return geno;
}

void print_array(int *array, int length)
{
    for(int i = 0; i < length; i++)
    {
        printf("%i ", array[i]);
    }
}

void free_z_score_list(node_z_score *head)
{
    node_z_score *current = head;
    while (current != NULL) 
    {
        node_z_score *next_node = current->next;
        free(current->snp_name);  // Free the SNP name string
        free(current);            // Free the node itself
        current = next_node;
    }
}

void free_mem_close_files(FILE *geno, FILE *pheno, folder_geno_pheno *data, node_z_score *head)
{
    if (geno != NULL) fclose(geno);
    if (pheno != NULL) fclose(pheno);

    if (data != NULL) {
        free(data->pheno_vector);
        free(data->geno_matrix);
        free(data->folder_name);
        free(data);
    }

    if (head != NULL)
        free_z_score_list(head);

}

void create_node_and_head_insert(node_z_score **head, char *snp_name )
{
    node_z_score *new_node = (node_z_score *)malloc(sizeof(node_z_score));
    
    new_node->snp_name = snp_name;
    new_node->real_z_score = 0.00;
    new_node->linear_z_score = 0.00;
    new_node->corrected_linear_z_score = 0.00;
    new_node->next = *head;

    *head = new_node;
}

void print_list(node_z_score *head)
{
    node_z_score *current = head;
    while(current != NULL)
    {
        printf("SNP name: %s \t real_z_score: %f \t linear_z_score: %f \t corrected_linear_z_score: %f \n", current->snp_name, current->real_z_score, current->linear_z_score, current->corrected_linear_z_score);
        current = current->next;
    }

}

void resolve_input_path(const char *target_dir, const char *given_path, char *out_abs, size_t out_size)
{
    char joined[PATH_MAX];

    if (given_path[0] == '/') {
        snprintf(joined, sizeof(joined), "%s", given_path);   // already absolute, use as-is
    } else {
        snprintf(joined, sizeof(joined), "%s/%s", target_dir, given_path); // relative -> inside target_dir, matching every other file in this pipeline
    }

    if (realpath(joined, out_abs) == NULL) {
        perror(joined);
        exit(1);
    }
}


int get_n_snps_in_file(FILE *geno_file)
{
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    int snp_count = 0;

    while ((len = getline(&line, &cap, geno_file)) != -1) {
        if (len > 0) snp_count++;
    }

    free(line);
    rewind(geno_file);
    return snp_count;
}

int *read_pheno_vector(FILE *pheno_file, int *pruned_sample_size, int **keep_mask, int *original_sample_size, int pheno_offset)
{
    char buffer[pheno_buffer];
    int n_indiv = 0;
    int mask_index = 0;
    int n_doublings_original = 1;
    int n_doublings_pruned = 1;

    int *pheno_vector = allocate_pheno_vector();

    const char *low_token  = pheno_offset ? "1" : "0"; // which two literal tokens count as valid, non-NA values
    const char *high_token = pheno_offset ? "2" : "1";


    while (fscanf(pheno_file, "%127s", buffer) == 1) {
        
        
        if(mask_index == n_doublings_original*max_rows)
        {
            int *temp = realloc(*keep_mask, sizeof(int)*2*n_doublings_original*max_rows);
            if(temp == NULL)
            {
                perror("Memory reallocation for keep mask failed");
                free(*keep_mask);
                exit(1);
            }
            *keep_mask = temp;
            n_doublings_original++;

        }
        
        
        if (strcmp(buffer, low_token) != 0 && strcmp(buffer, high_token) != 0)
        {
            (*keep_mask)[mask_index] = 0;
            mask_index++;
            continue;  // it's NA
        }  
        
        
        int value = atoi(buffer) - pheno_offset;     
        
        if(is_integer_and_positive(value))
        {
            
            if(n_indiv == n_doublings_pruned*max_rows)
            {
                int * temp = realloc(pheno_vector, sizeof(int)*2*n_doublings_pruned*max_rows);
                if (temp == NULL) 
                {
                    perror("Memory reallocation failed");
                    free(pheno_vector);
                    exit(1);
                }
                pheno_vector = temp;

                n_doublings_pruned++;
            }
            
            pheno_vector[n_indiv] = value;
            n_indiv++;
            (*keep_mask)[mask_index] = 1;  // Mark this sample as kept
        }
        else
            (*keep_mask)[mask_index] = 0;  // Mark this sample as skipped

        mask_index++;
        
    }

    *pruned_sample_size = n_indiv;
    *original_sample_size = mask_index;
    return pheno_vector;
}



int *read_geno_matrix(FILE *geno_file, int *pruned_sample_size, int *keep_mask, int *original_sample_size, int *snp_size, node_z_score **head, double **out_full_geno, int is_genotype_file)
{
    size_t n_snps = (size_t)(*snp_size);
    size_t n_orig = (size_t)(*original_sample_size);
    size_t n_pruned = (size_t)(*pruned_sample_size);

    int *geno_matrix = allocate_geno_matrix(*pruned_sample_size, *snp_size);

    double *full_geno = malloc(n_snps * n_orig * sizeof(double));
    if (!full_geno) 
    {
        perror("Full genotype memory allocation failed");
        exit(1);
    }
    *out_full_geno = full_geno;

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    int snp_index = 0;

    while ((len = getline(&line, &cap, geno_file)) != -1) 
    {
        char *saveptr = NULL;
        char *token = strtok_r(line, " ,\t\r\n", &saveptr);
        if (!token) continue;

        char *snp_name = strdup(token);
        create_node_and_head_insert(head, snp_name);

        if (is_genotype_file) 
        {                                       
            token = strtok_r(NULL, " ,\t\r\n", &saveptr); // ref allele
            token = strtok_r(NULL, " ,\t\r\n", &saveptr); // alt allele
        }
        token = strtok_r(NULL, " ,\t\r\n", &saveptr); // first value

        int original_sample_index = 0;
        int pruned_sample_index = 0;

        while (token != NULL && original_sample_index < *original_sample_size) {
            double g_val = -1.0;
            if (strcmp(token, "NA") != 0 && strcmp(token, "NaN") != 0) 
            {
                g_val = atof(token);
            }

            full_geno[(size_t)snp_index * n_orig + (size_t)original_sample_index] = g_val;

            if (keep_mask[original_sample_index]) 
            {
                int g_int = (g_val < 0.0) ? -1 : (int)(g_val + 0.5);
                geno_matrix[(size_t)snp_index * n_pruned + (size_t)pruned_sample_index] = g_int;
                pruned_sample_index++;
            }

            original_sample_index++;
            token = strtok_r(NULL, " ,\t\r\n", &saveptr);
        }

        snp_index++;
    }

    free(line);
    *snp_size = snp_index;
    return geno_matrix;
}



folder_geno_pheno *read_files(char *folder_path, const char *matrix_path, const char *outcome_path, int is_genotype_file, int pheno_offset, node_z_score **head){
    folder_geno_pheno *data = (folder_geno_pheno *)malloc(sizeof(folder_geno_pheno));
    if (data == NULL) 
    {
        perror("Memory allocation for folder_geno_pheno struct failed");
        exit(1);
    }

    int number_of_doublings = 1;
    data->folder_name = strdup(folder_path);
    data->skipped_samples = allocate_pheno_vector();
    data->pruned_sample_size = 0;
    data->original_sample_size = 0;
    data->snp_size = 0;

    int n_indiv = 0;



    FILE *pheno_file = fopen(outcome_path, "r");
    FILE *geno_file = fopen(matrix_path, "r");

    if (geno_file == NULL) 
    {
        perror(matrix_path);
        free_mem_close_files(NULL, pheno_file, data, *head);
        return NULL;
    }

    if (pheno_file == NULL) 
    {
        perror(outcome_path);
        free_mem_close_files(geno_file, NULL, data, *head);
        return NULL;
    }

    data->pheno_vector = read_pheno_vector(pheno_file, &data->pruned_sample_size, &data->skipped_samples, &data->original_sample_size, pheno_offset);
    
    data->snp_size = get_n_snps_in_file(geno_file);
    data->geno_matrix = read_geno_matrix(geno_file, &data->pruned_sample_size, data->skipped_samples, &data->original_sample_size, &data->snp_size, head, &data->full_geno_matrix, is_genotype_file);

    fclose(geno_file);
    fclose(pheno_file);
    return data;
}



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

void process_snp_row(double *row, int original_sample_size, int *skipped_samples,
                     double *mean_incl_skipped, double *std_incl_skipped,
                     double *mean_no_skipped, double *std_no_skipped)
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
    for (int i = 0; i < original_sample_size; ++i)
        std_row[i] = 0.0;

    if (std_no_skipped <= 0.0 || std <= 0.0) {
        *keep_mask = 0;
        return;
    }

    *keep_mask = 1;
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
        double *row = data->full_geno_matrix + j * data->original_sample_size;
 
        
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
                        std_matrix + j * data->original_sample_size);
        
        if(data->snp_keep_mask[j])
        {
            valid_snp_count++;
        }


    }    
 
    *n_valid_snps = valid_snp_count;
}



/*----------------------------------------------------------
  Randomly sample SNP row indices

  Input:
      snp_size          Total number of SNPs (rows)
      sample_fraction   e.g. 0.20 for 20%
      seed              Random seed

  Output:
      sampled_snp_size  Number of sampled SNPs
      return            Array of sampled row indices
----------------------------------------------------------*/
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



/*----------------------------------------------------------
  Extract only the sampled SNPs.

  Input:
      std_matrix          Full standardized matrix
                          Shape: snp_size × original_sample_size

      sampled_snp_indices Row indices returned above

  Output:
      subset_std_matrix
      Shape: sampled_snp_size × original_sample_size
----------------------------------------------------------*/

double *extract_subset_std_matrix(const double *std_matrix, int original_sample_size, const int *sampled_snp_indices, int sampled_snp_size)
{
    double *subset_std_matrix = (double *)malloc(sampled_snp_size * original_sample_size * sizeof(double));
    if (subset_std_matrix == NULL) 
    {
        perror("malloc subset_std_matrix");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sampled_snp_size; i++) 
    {
        memcpy(subset_std_matrix + i * original_sample_size, std_matrix + sampled_snp_indices[i] * original_sample_size, original_sample_size * sizeof(double));
    }

    return subset_std_matrix;
}


/*----------------------------------------------------------
  Compute the Genomic Relationship Matrix

      GRM = (1/m) X^TX

  Input:
      subset_std_matrix
      Shape: sampled_snp_size × original_sample_size

  Output:
      grm
      Shape: original_sample_size × original_sample_size
----------------------------------------------------------*/

void compute_grm_CALL(double *subset_std_matrix, int sampled_snp_size, int original_sample_size, double *grm, int m_divisor)
{
    cblas_dgemm(
        CblasRowMajor,
        CblasTrans,
        CblasNoTrans,
        original_sample_size,
        original_sample_size,
        sampled_snp_size,
        1.0,
        subset_std_matrix,
        original_sample_size,
        subset_std_matrix,
        original_sample_size,
        0.0,
        grm,
        original_sample_size
    );

    double inv = 1.0 / (double)m_divisor;
    for (int i = 0; i < original_sample_size * original_sample_size; i++) 
    {
        grm[i] *= inv;
    }
}



double *compute_GRM(const double *std_matrix, int snp_size, int original_sample_size, double sample_fraction, const int *valid_snp_indices, int n_valid_snps)
{
    double *GRM = (double *)calloc((size_t)original_sample_size * original_sample_size, sizeof(double));
    if (GRM == NULL) return NULL;


    if (sample_fraction == 1.0) {

        double *valid_std_matrix = extract_subset_std_matrix(std_matrix, original_sample_size, valid_snp_indices, n_valid_snps);

        cblas_dgemm(CblasRowMajor, 
                    CblasTrans, 
                    CblasNoTrans,
                    original_sample_size, 
                    original_sample_size, 
                    n_valid_snps,
                    1.0,
                    valid_std_matrix, 
                    original_sample_size,
                    valid_std_matrix, 
                    original_sample_size,
                    0.0,
                    GRM, 
                    original_sample_size);

        double inv = 1.0 / (double)n_valid_snps;
        for (int i = 0; i < original_sample_size * original_sample_size; i++) 
        {
            GRM[i] *= inv;
        }
        free(valid_std_matrix);
        return GRM;
    }

    int sampled_snp_size = 0;
    int *sampled_snp_indices = sample_snp_row_indices_from_valid(valid_snp_indices, n_valid_snps, sample_fraction, &sampled_snp_size, 42);

    if (sampled_snp_indices == NULL || sampled_snp_size <= 0) 
    {
        free(GRM);
        return NULL;
    }

    double *subset_std_matrix = extract_subset_std_matrix(std_matrix, original_sample_size, sampled_snp_indices, sampled_snp_size);

    compute_grm_CALL(subset_std_matrix, sampled_snp_size, original_sample_size, GRM, sampled_snp_size);

    free(sampled_snp_indices);
    free(subset_std_matrix);
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
    double *sub = malloc((size_t)pruned_n * pruned_n * sizeof(double));
    if (sub == NULL) {
        perror("malloc analyzed submatrix");
        exit(1);
    }

    int *kept_idx = malloc(pruned_n * sizeof(int));
    if (kept_idx == NULL) {
        perror("malloc kept_idx");
        exit(1);
    }

    int k = 0;
    for (int i = 0; i < full_n; i++) {
        if (keep_mask[i]) kept_idx[k++] = i;
    }

    for (int i = 0; i < pruned_n; i++) {
        for (int j = 0; j < pruned_n; j++) {
            sub[i * pruned_n + j] = full_matrix[kept_idx[i] * full_n + kept_idx[j]];
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
                    U[i * n + j]);
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
    for(int row = 0; row < n; row++) {
        for(int col = 0; col < n / 2; col++) {
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
        char exact_folder_path[8192];
        snprintf(exact_folder_path, sizeof(exact_folder_path), "%s/exact", data->folder_name);
        if(!is_directory(exact_folder_path))
        {
            if(mkdir(exact_folder_path, 0755) != 0)
            {
                perror("Failed to create folder for exact eigen decomposition");
                return -1;
            }
        }
        write_eigen_decomposition_to_file(grm, eigvals, pruned_sample_size, data->folder_name);
    }

    return info;
}





/* Mean of phenotype vector */
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


/*
  Computes:
      out = mean_y + U * D * U^T * (y - mean_y)

  where D[i] = eigvals[i] / (eigvals[i] + delta)

  U is row-major, size n x n
*/
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

    /* tmp = U^T * y_centered */
    for (int i = 0; i < pruned_sample_size; i++) {
        double s = 0.0;
        for (int j = 0; j < pruned_sample_size; j++) {
            s += U[j * pruned_sample_size + i] * y_centered[j];
        }
        tmp[i] = s;
    }

    /* tmp2[i] = D[i] * tmp[i] */
    for (int i = 0; i < pruned_sample_size; i++) {
        double d = eigvals[i] / (eigvals[i] + delta);
        tmp2[i] = d * tmp[i];
    }

    /* out = U * tmp2 */
    for (int i = 0; i < pruned_sample_size; i++) {
        double s = 0.0;
        for (int j = 0; j < pruned_sample_size; j++) {
            s += U[i * pruned_sample_size + j] * tmp2[j];
        }
        p[i] = s + mean_y;
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




double *compute_four_counts(folder_geno_pheno *data, int snp_index, double *counts, char *snp_col, double *surprise_scores)
{
    int snp_base = snp_index * data->pruned_sample_size;

    for (int i = 0; i < data->pruned_sample_size; i++) {
        int snp_value = data->geno_matrix[snp_base + i];

        if (snp_value < 0) 
        {
            snp_value = (int)(data->snp_means_no_skipped_samples[snp_index] + 0.5);
        }
        snp_col[i] = snp_value + '0';

        double w = surprise_scores[i];
        int pheno_value = data->pheno_vector[i];
        int not_pheno_value = !pheno_value;

        counts[0] += (double)snp_value * (double)pheno_value * w;
        counts[2] += (double)(2 - snp_value) * (double)(snp_value < 2) * (double)pheno_value * w;

        counts[1] += (double)snp_value * (double)not_pheno_value * w;
        counts[3] += (double)(2 - snp_value) * (double)(snp_value < 2) * (double)not_pheno_value * w;
    }
}



void write_geno_col_to_file(FILE *output_file, folder_geno_pheno *data, int snp_index, char *snp_name)
{
    // Write SNP header columns
    fprintf(output_file, "%s\tA\tT", snp_name);
    
    int base = snp_index * data->original_sample_size;
    
    for (int i = 0; i < data->original_sample_size; i++) {
        // Only output samples that pass the inclusion/keep mask
        if (data->skipped_samples[i]) { 
            double val = data->full_geno_matrix[base + i];
            
            if (val < 0.0) {
                fprintf(output_file, "\tNA"); // Restore original missing strings for GEMMA
            } else {
                fprintf(output_file, "\t%g", val); // Retain original float/dosage precision
            }
        }
    }
    fprintf(output_file, "\n");
}

void write_pruned_pheno_to_file(folder_geno_pheno *data) 
{
    char pheno_out_path[8192];
    snprintf(pheno_out_path, sizeof(pheno_out_path), "%s/exact/C_pruned_pheno.txt", data->folder_name);
    
    FILE *pf = fopen(pheno_out_path, "w");
    if (pf == NULL) {
        perror("Failed to create pruned phenotype file");
        exit(1);
    }
    
    for (int i = 0; i < data->pruned_sample_size; i++) {
        // Adding 1 back to keep GEMMA's expected 1/2 format
        fprintf(pf, "%d\n", data->pheno_vector[i] + 1); 
    }
    fclose(pf);
}


void write_specs_if_absent(const char *folder_path, int n_samples, int n_candidates, double z_threshold)
{
    char specs_path[8192];
    snprintf(specs_path, sizeof(specs_path), "%s/specs.txt", folder_path);

    if (access(specs_path, F_OK) == 0) return; // a reference run already exists, don't clobber it

    FILE *specs_file = fopen(specs_path, "w");
    if (specs_file == NULL) { perror(specs_path); return; }

    fprintf(specs_file, "n_samples\t%d\n", n_samples);
    fprintf(specs_file, "n_candidates\t%d\n", n_candidates);
    fprintf(specs_file, "z_score_threshold\t%.6f\n", z_threshold);

    fclose(specs_file);
}

int read_specs(const char *folder_path, int *ref_n_samples, int *ref_n_candidates, double *ref_z_threshold)
{
    char specs_path[8192];
    snprintf(specs_path, sizeof(specs_path), "%s/specs.txt", folder_path);

    FILE *specs_file = fopen(specs_path, "r");
    if (specs_file == NULL) return 0;

    char key[64];
    int ok = 1;
    ok &= (fscanf(specs_file, "%63s %d", key, ref_n_samples) == 2);
    ok &= (fscanf(specs_file, "%63s %d", key, ref_n_candidates) == 2);
    ok &= (fscanf(specs_file, "%63s %lf", key, ref_z_threshold) == 2);

    fclose(specs_file);
    return ok;
}


/* Bonferroni-equivalent thresholds (alpha = 0.05), chunked in half-decades of n.
   Value is taken at the lower edge of each chunk, so it is never anti-conservative. */
static double z_threshold_for_n(int n_snps)
{
    int k = (int)floor(2.0 * log10((double)n_snps));
    switch (k) {
        case 4:  return 3.29;   /*        100 -        316 */
        case 5:  return 3.60;   /*        316 -      1,000 */
        case 6:  return 3.89;   /*      1,000 -      3,162 */
        case 7:  return 4.16;   /*      3,162 -     10,000 */
        case 8:  return 4.41;   /*     10,000 -     31,623 */
        case 9:  return 4.65;   /*     31,623 -    100,000 */
        case 10: return 4.89;   /*    100,000 -    316,228 */
        case 11: return 5.11;   /*    316,228 -  1,000,000 */
        case 12: return 5.32;   /*  1,000,000 -  3,162,278 */
        case 13: return 5.53;   /*  3,162,278 - 10,000,000 */
        case 14: return 5.73;   /* 10,000,000 - 31,622,777 */
        case 15: return 5.92;   /* 31.6M - 100M */
        default:
            return (k < 4) ? 3.29   /* n < 100: not chunked, threshold is a floor here */
                           : 6.10;  /* n >= 100M: conservative lower edge */
    }
}

static double gemma_threshold_for_n(int n_snps)
{
    int k = (int)floor(2.0 * log10((double)n_snps));
    if (k < 4) k = 4;                 /* same floor as above */
    return 1.30 + 0.5 * k;            /* = floor(-log10(0.05 / 10^(k/2)), 2 dp) */
}



void compute_z_scores(folder_geno_pheno *data, double *surprise_scores, node_z_score **head, FILE *relevant_log_file, int exact)
{
    node_z_score *current = *head;
    int count_relevant = 0;

    FILE *output_file = NULL;
    if (exact) 
    {
        write_pruned_pheno_to_file(data);

        char output_path[8192];
        snprintf(output_path, sizeof(output_path), "%s/exact/C_pruned_geno_matrix.txt", data->folder_name);

        output_file = fopen(output_path, "w");
        if (output_file == NULL) 
        {
            perror("Failed to open exact pruned genotype output file");
            exit(EXIT_FAILURE);
        }
    }

    int ref_n_samples = 0, ref_n_candidates = 0;
    double ref_z_threshold = 0.0;
    int have_specs = read_specs(data->folder_name, &ref_n_samples, &ref_n_candidates, &ref_z_threshold);

    double reference_sample_size = (have_specs && ref_n_samples > 0) ? (double)ref_n_samples : data->original_sample_size; // fallback preserves original behavior when no baseline exists yet

    //double adj_thr = data->pruned_sample_size < 500?   4.45*sqrt(data->original_sample_size / reference_sample_size) : 4.45;
    //double adj_thr = data->pruned_sample_size < 500?   4.30*sqrt(data->original_sample_size / reference_sample_size) : 4.30;


    double base_thr = z_threshold_for_n(data->snp_size);
    double adj_thr  = data->pruned_sample_size < 500 ? base_thr * sqrt(data->original_sample_size / reference_sample_size): base_thr;


    write_specs_if_absent(data->folder_name, data->pruned_sample_size, data->snp_size, adj_thr);    
    
    int processed_snps = 0;
    for (int i = data->snp_size - 1; i >= 0; i--) 
    {
        
        node_z_score *node = current;
        current = current->next;

        if(data->snp_stds_no_skipped_samples[i] <= 1e-6 || data->snp_stds_including_skipped_samples[i] < 1e-6) 
        {
            continue;   // prune monomorphic / zero-variance SNPs

        }

        processed_snps++;
        double counts[4] = {0.0, 0.0, 0.0, 0.0};
        char snp_col[buffer_size];
        compute_four_counts(data, i, counts, snp_col, surprise_scores);

        double a = counts[0] + 0.5;
        double b = counts[1] + 0.5;
        double c = counts[2] + 0.5;
        double d = counts[3] + 0.5;

        node->real_z_score = fabs(log((a / b) / (c / d)) / sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d));

        if (node->real_z_score > adj_thr) 
        {
            count_relevant++;
            if (exact) 
            {
                write_geno_col_to_file(output_file, data, i, node->snp_name);
            }
        }
    }

    fprintf(relevant_log_file, "%s\tZSCORE\t%d\t%d\t%.2f%%\n", data->folder_name, count_relevant, processed_snps, (count_relevant / (double)processed_snps) * 100.0);

    fflush(relevant_log_file);
    printf("Number of relevant SNPs retrieved by EIGEN-SHRUNK-ZSCORE: %d out of %d kept SNPs. Percentage: %.2f%%\n", count_relevant, processed_snps, (count_relevant / (double)processed_snps) * 100.0);

    if (output_file) 
    {
        fflush(output_file);
        fclose(output_file);
    }
}

void write_to_file(const char *folder_path, node_z_score *head, int *snp_keep_mask, int snp_size)
{
    char output_path[8192];
    snprintf(output_path, sizeof(output_path), "%s/z_scores.txt", folder_path);

    FILE *output_file = fopen(output_path, "w");
    if (output_file == NULL) 
    {
        perror(output_path);
        return;
    }

    fprintf(output_file, "SNP_name\tReal_Z_Score\n");

    node_z_score *current = head;
    for (int i = snp_size - 1; i >= 0; i--) 
    {
        node_z_score *node = current;
        current = current->next;

        if (snp_keep_mask && !snp_keep_mask[i]) 
        {
            continue;
        }

        fprintf(output_file, "%s\t%f\n", node->snp_name, node->real_z_score);
    }

    fclose(output_file);
}


char **read_gemma_relevant_snps(const char *folder_path, double thr, FILE *relevant_log_file, int exact)
{
    int count_relevant = 0;
    int n_snps = 0;
    char gemma_output_path[8192];

    if(!exact)
        snprintf(gemma_output_path, sizeof(gemma_output_path), "%s/output/gemma_output.assoc.txt", folder_path);
    else
    {
        snprintf(gemma_output_path, sizeof(gemma_output_path), "%s/exact/output/gemma_output.assoc.txt", folder_path);
    }

    FILE *gemma_file = fopen(gemma_output_path, "r");
    if (gemma_file == NULL) {
        perror(gemma_output_path);
        return NULL;
    }

    char line[16384];

    /* Read header */
    if (fgets(line, sizeof(line), gemma_file) == NULL) 
    {
        fclose(gemma_file);
        return NULL;
    }

    int rs_col = -1;
    int p_col  = -1;

    int col = 0;
    char *token = strtok(line, "\t\r\n");

    while (token != NULL) 
    {

        if (strcmp(token, "rs") == 0)
            rs_col = col;

        if (strcmp(token, "p_lrt") == 0)
            p_col = col;

        token = strtok(NULL, "\t\r\n");
        col++;
    }

    if (rs_col == -1 || p_col == -1) 
    {
        fprintf(stderr, "Could not locate rs or p_lrt column in GEMMA output.\n");
        fclose(gemma_file);
        return NULL;
    }

    int capacity = 128;
    int n = 0;

    char **snps = malloc(capacity * sizeof(char *));
    if (snps == NULL) {
        fclose(gemma_file);
        return NULL;
    }

    /* Read data */
    while (fgets(line, sizeof(line), gemma_file) != NULL) {

        n_snps++;
        col = 0;
        char *rs = NULL;
        double p_lrt = NAN;

        token = strtok(line, "\t\r\n");

        while (token != NULL) {

            if (col == rs_col)
                rs = token;

            if (col == p_col)
                p_lrt = atof(token);

            token = strtok(NULL, "\t\r\n");
            col++;
        }

        if (rs == NULL || !isfinite(p_lrt))
            continue;

        if (-log10(p_lrt) > thr) {

            count_relevant++;

            if (n == capacity) {
                capacity *= 2;

                char **tmp = realloc(snps, capacity * sizeof(char *));
                if (tmp == NULL)
                    break;

                snps = tmp;
            }

            snps[n++] = strdup(rs);
        }
    }

    fclose(gemma_file);

    snps = realloc(snps, (n + 1) * sizeof(char *));
    snps[n] = NULL;

    fflush(relevant_log_file);

    if(!exact)
    {
        printf("Number of relevant SNPs  retrieved by GEMMA alone: %d \t out of %d in file. \t Percentage: %.2f%%\n", count_relevant, n_snps, (count_relevant / (double)n_snps) * 100);
        fprintf(relevant_log_file, "%s\tGEMMA\t%d\t%d\t%.2f%%\n\n", folder_path, count_relevant, n_snps, (count_relevant / (double)n_snps) * 100);
    }
    else
    {
        printf("Number of relevant SNPs  retrieved by GEMMA after Z-score filtering: %d \t out of %d in file.\n\n", count_relevant, n_snps);
        fprintf(relevant_log_file, "%s\tGEMMA_AFTER_ZSCORE\t%d\t%d\n\n\n\n", folder_path, count_relevant, n_snps);
    }

    return snps;
}


void match_z_scores_with_gemma(node_z_score *z_score_list, char **gemma_relevant, FILE *fn_log_file, folder_geno_pheno *data) 
{
    node_z_score *current = z_score_list;

    for (int i = 0; gemma_relevant[i] != NULL; i++) 
    {
        int found = 0;
        current = z_score_list;  // Reset to the head of the list for each GEMMA SNP
        while (current != NULL) 
        {
            if (strcmp(current->snp_name, gemma_relevant[i]) == 0) 
            {
                found = 1;
            }
            current = current->next;
        }

        if (!found) 
        {
            printf("GEMMA SNP %s not found in z_score_list.\n", gemma_relevant[i]);
            fprintf(fn_log_file, "%s\t%s\n", data->folder_name, gemma_relevant[i]);
            fflush(fn_log_file);
        }
    }
    fprintf(fn_log_file, "\n\n");
}


void match_gemmas( char ** after_z_gemma_relevant, char ** gemma_relevant, FILE *fn_log_file, folder_geno_pheno *data)
{
    int total_not_found = 0;
    for (int i = 0; gemma_relevant[i] != NULL; i++) 
    {
        int found = 0;
        for (int j = 0; after_z_gemma_relevant[j] != NULL; j++) 
        {
            if (strcmp(gemma_relevant[i], after_z_gemma_relevant[j]) == 0) 
            {
                found = 1;
                break;
            }
        }

        if (!found) 
        {
            printf("Original GEMMA SNP %s not found in after Z-score GEMMA relevant SNPs.\n", gemma_relevant[i]);
            fprintf(fn_log_file, "%s\t%s\n", data->folder_name, gemma_relevant[i]);
            fflush(fn_log_file);
            total_not_found++;
        }
    }

    if(!total_not_found)
    {
        printf("All Original GEMMA relevant SNPs have been found against after Z-score GEMMA relevant SNPs.\n\n\n\n");
    }
    else{
        printf("\n\n\n\n");
    }
    fprintf(fn_log_file, "\n\n");
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
            fprintf(kinship_file, "%.8g", GRM[i * original_sample_size + j]);
            if (j < original_sample_size - 1) {
                fprintf(kinship_file, "\t");
            }
        }
        fprintf(kinship_file, "\n");
    }

    fclose(kinship_file);
}






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
    int pheno_offset = 0; // default: phenotype file already 0/1, no shift


    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            matrix_path = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outcome_path = argv[++i];
        } else if (strcmp(argv[i], "-g") == 0) {
            is_genotype_file = 1;
        } 
        else if (strcmp(argv[i], "-2") == 0) {
            pheno_offset = 1; // phenotype file is 1/2-coded, subtract 1 to normalize
        } else if (target_dir == NULL) {
            target_dir = argv[i];
        } else if (exact == -1) {
            exact = atoi(argv[i]);
        }
    }

    if (target_dir == NULL || exact == -1 || matrix_path == NULL || outcome_path == NULL) {
        fprintf(stderr, "Usage: %s <folder_path> <0=approximate|1=exact> -m <matrix.txt> -o <outcome.txt> [-g] [-2]\n", argv[0]);
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
    // Construct gemma commands
    char gemma_cmd1[16384];
    char gemma_cmd2[16384];


    snprintf(gemma_cmd1, sizeof(gemma_cmd1), 
        " cd %s && ../../gemma -g %s -p %s -gk 2 -maf 0 -miss 1 -o kinship > /dev/null 2>&1",
        target_dir, matrix_abs, outcome_abs);
    //snprintf(gemma_cmd2, sizeof(gemma_cmd2), 
    //    "cd %s && /home/francesco/anaconda3/envs/gwas_test/bin/gemma -g ./matrix.txt -p ./outcome_vector.txt -a ./anno.subset.txt -k output/kinship.sXX.txt -lmm 4 -n 1 -o gemma_output > /dev/null 2>&1", target_dir); // ");

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
    printf("Time taken for full and exact GEMMA execution: %f seconds\n\n", (t_gemma_end - t0));
    fprintf(time_log_file, "%s\tGEMMA\t%f\tms\n", target_dir, (t_gemma_end - t0));
    

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// ------------------------------------------------ EIGEN-SHRUNK Z-SCORE APPROACH ---------------------------------------------------------------------------------------------

    double t_zscore_start = now_seconds();



    // allocate keep mask for SNP pruning 
    data->snp_keep_mask = (int *)calloc((size_t)data->snp_size, sizeof(int));
    if (data->snp_keep_mask == NULL) 
    {
        perror("calloc snp_keep_mask");
        exit(EXIT_FAILURE);
    }

    // standardization + monomorphic SNP pruning
    double *std_matrix = (double *)malloc((size_t)data->snp_size * data->original_sample_size * sizeof(double));
    if (std_matrix == NULL) {
        perror("malloc std_matrix");
        exit(1);
    }

    int n_valid_snps = 0;
    standardize_snp_matrix(data, std_matrix, &n_valid_snps);

    data->n_valid_snps = n_valid_snps;
    data->valid_snp_indices = collect_valid_snp_indices_for_GRM(data->snp_keep_mask, data->snp_size, &data->n_valid_snps);

    //printf("%d %d\n", data->snp_size, n_valid_snps);

    // GRM uses only valid SNPs 
    double *GRM = compute_GRM(std_matrix, data->snp_size, data->original_sample_size, 1.0, data->valid_snp_indices, data->n_valid_snps);

    if (GRM == NULL) 
    {
        fprintf(stderr, "compute_GRM failed\n");
        free(std_matrix);
        exit(1);
    }

    //write_kinship(target_dir, GRM, data->original_sample_size);

    double *analyzed_GRM = extract_analyzed_submatrix(GRM, data->original_sample_size, data->skipped_samples, data->pruned_sample_size);

    double *eigvals = (double *)malloc(data->pruned_sample_size * sizeof(double));
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

    double *spectral_pheno = spectral_transform_phenotype(data->pheno_vector, analyzed_GRM, eigvals, data->pruned_sample_size, 1);

    double *surprise_scores = compute_clipped_IPW_surprise_score(spectral_pheno, data->pheno_vector, data->pruned_sample_size);

    compute_z_scores(data, surprise_scores, &z_score_list, relevant_log_file, exact);
    write_to_file(target_dir, z_score_list, data->snp_keep_mask, data->snp_size);


    double t_zscore_end = now_seconds();

    fprintf(time_log_file, "%s\tZ_SCORE\t%f\tms\n", target_dir, (t_zscore_end - t_zscore_start));
    printf("Time taken for EIGEN-SHRUNK Z-SCORE APPROACH on folder %s: %f milliseconds\n", target_dir, (t_zscore_end - t_zscore_start));

    match_z_scores_with_gemma(z_score_list, gemma_relevant, fn_log_file, data);

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




// ------------------------------------------------ EXACT GEMMA AFTER Z-SCORE APPROACH ---------------------------------------------------------------------------------------------

    char **after_z_gemma_relevant = NULL;

    if(exact)
    {
        char exact_after_zscore_cmd[16384];
        //snprintf(exact_after_zscore_cmd, sizeof(exact_after_zscore_cmd), "cd %s/exact && /home/francesco/anaconda3/envs/gwas_test/bin/gemma -g ./C_pruned_geno_matrix.txt -p ./C_pruned_pheno.txt -u eigenvectors.txt -d eigenvalues.txt -lmm 4 -n 1 -o gemma_output > /dev/null 2>&1", target_dir);
        snprintf(exact_after_zscore_cmd, sizeof(exact_after_zscore_cmd), "cd %s/exact && ../../../gemma -g ./C_pruned_geno_matrix.txt -p ./C_pruned_pheno.txt -u eigenvectors.txt -d eigenvalues.txt -maf 0 -miss 1 -lmm 4 -n 1 -o gemma_output > /dev/null 2>&1", target_dir);
        int result = system(exact_after_zscore_cmd);
        if(result == -1) {
            perror("Failed to execute GEMMA exact command after z-score computation");
            exit(EXIT_FAILURE); 
        }

        double t_exact_end = now_seconds();
        fprintf(time_log_file, "%s\tGEMMA_EXACT_AFTER_ZSCORE\t%f\tms\n", target_dir, (t_exact_end - t_zscore_end));
        printf("Time taken for exact GEMMA execution after EIGEN-SHRUNK Z-SCORE APPROACH on folder %s: %f milliseconds\n", target_dir, (t_exact_end - t_zscore_end));
        printf("Total time taken for EIGEN-SHRUNK Z-SCORE APPROACH and exact GEMMA execution on remaining positions %s: %f milliseconds \n", target_dir, (t_exact_end - t_zscore_start));
        fprintf(time_log_file, "%s\tZSCORE+GEMMA\t%f\tms\n\n\n\n", target_dir, (t_exact_end - t_zscore_start));


        after_z_gemma_relevant = read_gemma_relevant_snps(target_dir, gemma_thr, relevant_log_file, 1);
        match_gemmas(after_z_gemma_relevant, gemma_relevant, fn_log_file, data);

        free(after_z_gemma_relevant);
    }
    
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    if (data != NULL) 
    {
        free(data->pheno_vector);
        free(data->folder_name);
        free(data->geno_matrix);
        free(data->full_geno_matrix); 
        free(data->skipped_samples);
        free(data->snp_keep_mask);
        free(data->valid_snp_indices);
        free(std_matrix);
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


