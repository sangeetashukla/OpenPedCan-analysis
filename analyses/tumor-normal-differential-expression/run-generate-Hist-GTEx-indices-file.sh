#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Module author: Sangeeta Shukla, Alvin Farrell
# Shell script author: Sangeeta Shukla
# 2021

#This script creates a subset of the histologies.tsv, to use for testing the deseq module.

#module load R/4.1.0
Rscript --vanilla run-generate-Hist-GTEx-indices-file.R \
        --hist_file ../../data/v10/histologies.tsv \
        --counts_file ../../data/v10/gene-counts-rsem-expected_count-collapsed.rds \
        --outdir Input_Data \
	--ind_allcohorts ../../data/v10/independent-specimens.rnaseq.primary.tsv \
	--ind_eachcohort ../../data/v10/independent-specimens.rnaseq.primary.eachcohort.tsv
