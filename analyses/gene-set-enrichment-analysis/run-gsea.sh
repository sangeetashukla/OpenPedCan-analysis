#########################################################################
# Stephanie J. Spielman for ALSF CCDL 2020
#
# Run the GSEA pipeline:
## 1. `01-conduct-gsea-analysis.R` to calculate scores
## 2. `02-exploratory-gsea.Rmd` to explore the scores, lightly for now
#
# Usage: bash run-gsea.sh
#
# Takes one environment variable, `OPENPBTA_TESTING`, which is 1 for running
# samples in CI for testing, 0 for running the full dataset (Default)
#
# Takes one environment variable, `BASE_SUBTYPING`, if value is 1 then
# only runs modules required for subtyping if value is 0 runs all modules (Default)
#########################################################################


set -e
set -o pipefail


IS_CI=${OPENPBTA_TESTING:-0}
RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

if [ "$RUN_FOR_SUBTYPING" == 0 ]; then
DATA_DIR="../../data"
RESULTS_DIR="results"

######## Calculate scores from expression data ############
INPUT_FILE="${DATA_DIR}/gene-expression-rsem-tpm-collapsed.rds"
OUTPUT_FILE="${RESULTS_DIR}/gsva_scores.tsv"
HIST_FILE="${DATA_DIR}/gsva_scores.tsv"

Rscript --vanilla 01-conduct-gsea-analysis.R --input ${INPUT_FILE} --output ${OUTPUT_FILE} --histology ${HIST_FILE}

######## Model GSVA scores ############
# Only run when pbta-histologies.tsv is generated which has harmonized_diagnosis
Rscript -e "rmarkdown::render('02-model-gsea.Rmd', clean = TRUE, params=list(is_ci = ${IS_CI}))"
else
DATA_DIR="../collapse-rnaseq/results"
RESULTS_DIR="results"


fi
