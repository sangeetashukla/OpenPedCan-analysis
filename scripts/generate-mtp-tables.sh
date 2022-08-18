#!/bin/bash
# Eric Wafula, OPenPedCan 2022
# Run OpenPedCan modules to generate MTP tables

set -e
set -o pipefail

printf "Start generating mtp tables...\n\n"

# This script should always run as if it were being called from
# the directory it lives in.
SCRIPT_DIR="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$SCRIPT_DIR" || exit

# Set directories for modules that generate mtp tables
FUSION_MODULE_DIR=$SCRIPT_DIR/../analyses/fusion-frequencies
SNV_MODULE_DIR=$SCRIPT_DIR/../analyses/snv-frequencies
CNV_MODULE_DIR=$SCRIPT_DIR/../analyses/cnv-frequencies

# Run fusion-frequencies module
printf "\n\nGenerating mtp fusion frequencies tables...\n\n"
cd $FUSION_MODULE_DIR
bash run-frequencies.sh

# Run snv-frequencies module
printf "\n\nGenerating mtp snv frequencies tables...\n\n"
cd $SNV_MODULE_DIR
bash run-snv-frequencies.sh

# Run cnv-frequencies module
printf "\n\nGenrating mtp cnv frequencies tables...\n\n"
cd $CNV_MODULE_DIR
bash run-cnv-frequencies-analysis.sh

cd $SCRIPT_DIR

printf "\nDone generating mtp tables...\n\n"
