#!/bin/bash
# Module author: Komal S. Rathi
# 2022

# This script runs the steps for immune deconvolution using xCell and CIBERSORT (absolute mode). 

set -e
set -o pipefail

# Run original files - will not by default
RUN_CIBERSORT=${RUN_CIBERSORT:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# xCell
# generate deconvolution output
echo "Deconvolution..."
Rscript --vanilla 01-immune-deconv.R \
--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
--clin_file '../../data/histologies.tsv' \
--deconv_method 'xcell' \
--output_dir 'results'

# generate heatmaps of average normalized immune scores per cancer or gtex group
echo "Create summary plots"
Rscript --vanilla 02-summary-plots.R \
--deconv_output 'results/xcell_output.rds' \
--output_dir 'plots'

# CIBERSORT (absolute mode)
if [ "$RUN_CIBERSORT" -gt "0" ]; then

	# generate deconvolution output
	echo "Deconvolution..."
	Rscript --vanilla 01-immune-deconv.R \
	--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
	--clin_file '../../data/histologies.tsv' \
	--deconv_method 'cibersort_abs' \
	--cibersort_mat 'util/LM22.txt' \
	--cibersort_binary 'util/CIBERSORT.R' \
	--output_dir 'results'

	# generate heatmaps of average normalized immune scores per cancer or gtex group
	echo "Create summary plots"
	Rscript --vanilla 02-summary-plots.R \
	--deconv_output 'results/cibersort_abs_output.rds' \
	--output_dir 'plots'

fi