# #!/bin/bash

# # Define your list of target folders
# folders=(
#   "simulated_100010_snps_OR1.5_1000cases_1000controls"
#   "simulated_100010_snps_OR1.5_2000cases_2000controls"
#   "simulated_100010_snps_OR1.5_3000cases_3000controls"
#   "simulated_100010_snps_OR1.5_4000cases_4000controls"
#   "simulated_100010_snps_OR1.5_5000cases_5000controls"
#   "simulated_10010_snps_OR1.5_1000cases_1000controls"
#   "simulated_10010_snps_OR1.5_2000cases_2000controls"
#   "simulated_10010_snps_OR1.5_3000cases_3000controls"
#   "simulated_10010_snps_OR1.5_4000cases_4000controls"
#   "simulated_10010_snps_OR1.5_5000cases_5000controls"
#   "simulated_20010_snps_OR1.5_1000cases_1000controls"
#   "simulated_20010_snps_OR1.5_2000cases_2000controls"
#   "simulated_20010_snps_OR1.5_3000cases_3000controls"
#   "simulated_20010_snps_OR1.5_4000cases_4000controls"
#   "simulated_20010_snps_OR1.5_5000cases_5000controls"
#   "simulated_30010_snps_OR1.5_1000cases_1000controls"
#   "simulated_30010_snps_OR1.5_2000cases_2000controls"
#   "simulated_30010_snps_OR1.5_3000cases_3000controls"
#   "simulated_30010_snps_OR1.5_4000cases_4000controls"
#   "simulated_30010_snps_OR1.5_5000cases_5000controls"
#   "simulated_50010_snps_OR1.5_1000cases_1000controls"
#   "simulated_50010_snps_OR1.5_2000cases_2000controls"
#   "simulated_50010_snps_OR1.5_3000cases_3000controls"
#   "simulated_50010_snps_OR1.5_4000cases_4000controls"
#   "simulated_50010_snps_OR1.5_5000cases_5000controls"
# )

# # Python script location
# SCRIPT="run_subsets_saige.py"

# for dir in "${folders[@]}"; do
#     dir="${dir%/}"

#     echo "----------------------------------------"
#     echo "Processing: $dir"
#     echo "----------------------------------------"

#     if python3 "$SCRIPT" \
#         --matrix "$dir/sim_cc_matrix.txt" \
#         --pheno "$dir/sim_cc_outcome.txt" \
#         --binary ./z_score_general \
#         --workdir ./workdir \
#         --results ./results_approx \
#         --removal 10,30,50,70 \
#         --removed-cases 25,50,75 \
#         --replicates 3 \
#         --seed 42 \
#         --cores 40 \
#         -g \
#         --reuse \
#        --exact 0 

#     then
#         # Move results into the analyzed folder
#         if [ -d "results_approx" ]; then
#             if [ -e "$dir/results_approx" ]; then
#                 echo "ERROR: $dir/results_approx already exists."
#                 echo "Refusing to overwrite it."
#                 exit 1
#             fi

#             if mv results_approx "$dir/"; then
#                 echo "Moved results_approx -> $dir/results_approx"
#             else
#                 echo "ERROR: Failed to move results_approx -> $dir/results_approx"
#                 exit 1
#             fi
#         else
#             echo "ERROR: No 'results_approx' directory found."
#             exit 1
#         fi

#         # Remove workdir
#         if [ -d "workdir" ]; then
#             rm -rf workdir
#             echo "Cleaned up workdir."
#         fi
#     else
#         echo "Error: Execution failed for $dir. Skipping file move and cleanup."
#     fi
# done


#!/bin/bash

# Define your list of target folders
folders=(
  "simulated_10020_snps_OR1.5_1000cases_1000controls"
  "simulated_10020_snps_OR1.5_800cases_800controls"
)

# Python script location
SCRIPT="run_subsets_saige.py"

for dir in "${folders[@]}"; do
    dir="${dir%/}"

    echo "----------------------------------------"
    echo "Processing: $dir"
    echo "----------------------------------------"

    if python3 "$SCRIPT" \
        --matrix "$dir/sim_cc_matrix.txt" \
        --pheno "$dir/sim_cc_outcome.txt" \
        --binary ./z_score_general \
        --workdir ./workdir \
        --results ./results_approx \
        --removal 10,30,50,70 \
        --removed-cases 25,50,75 \
        --replicates 1 \
        --seed 42 \
        --cores 10 \
        -g \
        --reuse \
       --exact 0 \
       --no-baseline

    then
        # Move results into the analyzed folder
        if [ -d "results_approx" ]; then
            if [ -e "$dir/results_approx" ]; then
                echo "ERROR: $dir/results_approx already exists."
                echo "Refusing to overwrite it."
                exit 1
            fi

            if mv results_approx "$dir/"; then
                echo "Moved results_approx -> $dir/results_approx"
            else
                echo "ERROR: Failed to move results_approx -> $dir/results_approx"
                exit 1
            fi
        else
            echo "ERROR: No 'results_approx' directory found."
            exit 1
        fi

        # Remove workdir
        if [ -d "workdir" ]; then
            rm -rf workdir
            echo "Cleaned up workdir."
        fi
    else
        echo "Error: Execution failed for $dir. Skipping file move and cleanup."
    fi
done