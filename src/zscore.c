#include "raft.h"

/* weighted association score, thresholds, and Z-score output. */

void compute_four_counts(folder_geno_pheno *data, int snp_index, double *counts, double *surprise_scores)
{
    const double *row = data->full_geno_matrix + (size_t)snp_index * (size_t)data->original_sample_size;
    int imputed = (int)(data->snp_means_no_skipped_samples[snp_index] + 0.5);
    int k = 0;   // index into the pruned vectors

    for (int i = 0; i < data->original_sample_size; i++) {
        if (!data->skipped_samples[i]) continue;

        double g = row[i];
        int snp_value = (g < 0.0) ? imputed : (int)(g + 0.5);

        double w = surprise_scores[k];
        int pheno_value = data->pheno_vector[k];
        int not_pheno_value = !pheno_value;
        k++;

        counts[0] += (double)snp_value * (double)pheno_value * w;
        counts[2] += (double)(2 - snp_value) * (double)(snp_value < 2) * (double)pheno_value * w;
        counts[1] += (double)snp_value * (double)not_pheno_value * w;
        counts[3] += (double)(2 - snp_value) * (double)(snp_value < 2) * (double)not_pheno_value * w;
    }
}



// Bonferroni-equivalent thresholds (alpha = 0.05), chunked in half-decades of n.

double z_threshold_for_n(int n_snps)
{
    int k = (int)floor(2.0 * log10((double)n_snps));
    switch (k) {
        case 4:  return 3.29;   //        100 -        316 
        case 5:  return 3.60;   //        316 -      1,000 
        case 6:  return 3.89;   //      1,000 -      3,162 
        case 7:  return 4.16;   //      3,162 -     10,000 
        case 8:  return 4.41;   //     10,000 -     31,623 
        case 9:  return 4.65;   //     31,623 -    100,000 
        case 10: return 4.89;   //    100,000 -    316,228 
        case 11: return 5.11;   //    316,228 -  1,000,000 
        case 12: return 5.32;   //  1,000,000 -  3,162,278 
        case 13: return 5.53;   //  3,162,278 - 10,000,000 
        case 14: return 5.73;   // 10,000,000 - 31,622,777 
        case 15: return 5.92;   // 31.6M - 100M 
        default:
            return (k < 4) ? 3.29   // n < 100: not chunked, threshold is a floor here
                           : 6.10;  // n >= 100M: conservative lower edge 
    }
}

double gemma_threshold_for_n(int n_snps)
{
    int k = (int)floor(2.0 * log10((double)n_snps));
    if (k < 4) k = 4;                 // same floor as above 
    return 1.30 + 0.5 * k;            // = floor(-log10(0.05 / 10^(k/2)), 2 dp) 
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

    double base_thr = z_threshold_for_n(data->snp_size);
    data->adjusted_z_score_threshold = sqrt(data->pruned_sample_size / reference_sample_size)* base_thr;

    write_specs_if_absent(data->folder_name, data->pruned_sample_size, data->snp_size, data->adjusted_z_score_threshold);    
    
    int processed_snps = 0;
    for (int i = data->snp_size - 1; i >= 0; i--) 
    {
        
        node_z_score *node = current;
        current = current->next;

        if(data->snp_stds_no_skipped_samples[i] <= 1e-6 || data->snp_stds_including_skipped_samples[i] < 1e-6) 
        {
            continue;   // prune monomorphic / zero-variance 

        }

        processed_snps++;
        double counts[4] = {0.0, 0.0, 0.0, 0.0};
        compute_four_counts(data, i, counts, surprise_scores);

        double a = counts[0] + 0.5;
        double b = counts[1] + 0.5;
        double c = counts[2] + 0.5;
        double d = counts[3] + 0.5;

        node->real_z_score = fabs(log((a / b) / (c / d)) / sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d));

        if (node->real_z_score > data->adjusted_z_score_threshold) 
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

