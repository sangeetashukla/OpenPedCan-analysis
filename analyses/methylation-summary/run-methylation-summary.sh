#!/bin/bash
# PediatricOpenTargets 2021
# Eric Wafula
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up path for results directory
results_path="results"

####################### Create methylation array probe annotations ########################
Rscript --vanilla 01-create-methylation-annotation-table.R

######################## Create Beta-value  and M-values matrices #########################
Rscript --vanilla 02-create-methylation-matrices.R

###################### Calculate probe-level Beta-value quantiles #########################
Rscript --vanilla 03-calculate-beta-quantiles.R

############ Calculate probe-level correlations between Beta and TPM values ###############
Rscript --vanilla 04-target-beta-tpm-correlation.R

############################ Create methylation summary table #############################
Rscript --vanilla 05-create-methylation-summary-table.R

########################### Convert JSON to JSON Lines (JSONL) ############################
printf '\nConverting JSON file to JSONL files...\n'

jq --compact-output '.[]' \
  ${results_path}/methylation-beta-values-summary.json \
  > ${results_path}/methylation-beta-values-summary.jsonl

################################# Remove JSON file ########################################
printf '\nRemoving JSON files...\n'

rm ${results_path}/*.json

################################## Compress JSONL files ###################################
printf '\nCompressing JSONL files...\n'

if ls ${results_path}/*.jsonl.gz &>/dev/null
then 
  rm ${results_path}/*.jsonl.gz 
fi

gzip -v --no-name ${results_path}/*.jsonl

printf '\nAnalysis Done...\n'