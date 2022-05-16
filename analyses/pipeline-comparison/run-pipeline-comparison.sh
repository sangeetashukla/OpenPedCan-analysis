#!/bin/bash
# PediatricOpenTargets 2022
# Sangeeta Shukla
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


# Set up paths to data files consumed by analysis
data_path='input'


# RNA-Seq counts files from varying workflows needed 
name_wf1='KF'
name_wf2='JAX'
counts_wf1="${data_path}/PPTC-KF-gene-expression-rsem-tpm-collapsed-matrix.tsv"
counts_wf2="${data_path}/PPTC-JAX-gene-expression-rsem-tpm-collapsed-matrix.tsv"

# homologs file
homologs="${data_path}/Gene_ID_Matches.tsv"

#Create functions for utility
#Rscript utils/compare_RNASeq_pipelines.R

out_name="results"+${name_wf1}+"_"+${name_wf2}+"_comparison"

#Run QC and comparison on the input RNA-Seq counts files
Rscript 01-pipeline-comparison-all-genes.R --wf1_name $name_wf1 \
--wf1_file $counts_wf1 \
--wf2_name $name_wf2 \
--wf2_file $counts_wf2 \
--input_homologs $homologs \
--output_filename $out_name
