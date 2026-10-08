#include "raft.h"

/* GEMMA result parsing and RAFT/GEMMA matching. */

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
                if (current->real_z_score > data->adjusted_z_score_threshold) {
                    found = 1;
                    break;
                }
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


