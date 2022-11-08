#!/bin/bash

# OPenPedCan 2022
# Eric Wafula
set -e
set -o pipefail

printf "Start QC and Summary checks...\n\n"


# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


# Set up paths to output directory for the analysis
results_dir="results"

######################### Checking gene wise table ##########################
printf "Checking TPM Tables...\n\n"

printf "Checking gene wise TPM Tables...\n\n"

Rscript -e "rmarkdown::render(
            '02-check.Rmd',
            params = list(
              current_table_gene = 'current/long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz',
              previous_table_gene = 'previous/long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz'),
              clean = TRUE)"  
              
              
# remove intermediate files
rm ${results_dir}/*.tsv

printf "Checking group wise TPM Tables...\n\n"

Rscript -e "rmarkdown::render(
            '02-check.Rmd',
            params = list(
              current_table_gene = 'current/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz',
              previous_table_gene = 'previous/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz'),
              clean = TRUE)"  
              
              
# remove intermediate files
rm ${results_dir}/*.tsv

              
printf "Analysis Done...\n"  
