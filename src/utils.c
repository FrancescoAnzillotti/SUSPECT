#include "raft.h"

/* generic utilities, allocation, linked-list helpers, and path handling. */

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
        //free(data->geno_matrix);
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
    new_node->next = *head;

    *head = new_node;
}

void print_list(node_z_score *head)
{
    node_z_score *current = head;
    while(current != NULL)
    {
        printf("SNP name: %s \t z_score: %f \n", current->snp_name, current->real_z_score);
        current = current->next;
    }

}

void resolve_input_path(const char *target_dir, const char *given_path, char *out_abs, size_t out_size)
{
    char joined[PATH_MAX];

    if (given_path[0] == '/') {
        snprintf(joined, sizeof(joined), "%s", given_path);   // already absolute
    } else {
        snprintf(joined, sizeof(joined), "%s/%s", target_dir, given_path); // relative -> inside target_dir, matching every other file in this pipeline
    }

    if (realpath(joined, out_abs) == NULL) {
        perror(joined);
        exit(1);
    }
}


