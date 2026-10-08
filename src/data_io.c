#include "raft.h"

/* input parsing and filesystem output helpers. */

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



void *read_geno_matrix(FILE *geno_file, int *pruned_sample_size, int *keep_mask, int *original_sample_size, int *snp_size, node_z_score **head, double **out_full_geno, int is_genotype_file)
{
    size_t n_snps = (size_t)(*snp_size);
    size_t n_orig = (size_t)(*original_sample_size);
    size_t n_pruned = (size_t)(*pruned_sample_size);

    //int *geno_matrix = allocate_geno_matrix(*pruned_sample_size, *snp_size);

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

            // if (keep_mask[original_sample_index]) 
            // {
            //     int g_int = (g_val < 0.0) ? -1 : (int)(g_val + 0.5);
            //     geno_matrix[(size_t)snp_index * n_pruned + (size_t)pruned_sample_index] = g_int;
            //     pruned_sample_index++;
            // }

            original_sample_index++;
            token = strtok_r(NULL, " ,\t\r\n", &saveptr);
        }

        snp_index++;
    }

    free(line);
    *snp_size = snp_index;
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
    read_geno_matrix(geno_file, &data->pruned_sample_size, data->skipped_samples, &data->original_sample_size, &data->snp_size, head, &data->full_geno_matrix, is_genotype_file);

    fclose(geno_file);
    fclose(pheno_file);
    return data;
}



void write_geno_col_to_file(FILE *output_file, folder_geno_pheno *data, int snp_index, char *snp_name)
{
    // Write SNP header columns
    fprintf(output_file, "%s\tA\tT", snp_name);
    
    size_t base = (size_t)snp_index * (size_t)data->original_sample_size;
    
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
