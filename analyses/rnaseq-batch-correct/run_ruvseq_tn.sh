# T-N differences
#!/bin/bash
# Module authors: Adam Kraya
# 2022

# This script runs the steps for DESeq2 tumor-only analysis with and without RUVg batch correction.

set -e
set -o pipefail

Rscript code/02-ruvseq-deseq-tn.R \
--dataset 'HGG-DMG_TN' \
--cohort_values 'PBTA,GTEx' \
--cancer_group_values 'Diffuse midline glioma,High-grade glioma/astrocytoma' \
--gtex_subgroups 'Brain - Cortex,Brain - Cerebellum' \
--k_value 5 \
--drop 2

Rscript code/03-ruvseq-summarization.R \
--dataset 'HGG-DMG_TN' \
--cohort_values 'PBTA,GTEx' \
--cancer_group_values 'Diffuse midline glioma,High-grade glioma/astrocytoma' \
--ruvg_dir 'deseq2_analysis' \
--pos_c 'reactome_cell_cycle_msigdb_v7.5.1.rds' \
--neg_c 'hk_genes_normals.rds' \
--analysis_type 'tumor-normal' \
--gtex_subgroups 'Brain - Cortex,Brain - Cerebellum'
