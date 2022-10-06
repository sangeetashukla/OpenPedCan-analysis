# molecular subtype differences
#!/bin/bash
# Module authors: Komal S. Rathi, Adam Kraya
# 2022

# This script runs the steps for DESeq2 tumor-only analysis with and without RUVg batch correction.

set -e
set -o pipefail

Rscript code/01-ruvseq-deseq.R \
--dataset 'target_nbl' \
--cohort_values 'TARGET' \
--cancer_group_values 'Neuroblastoma' \
--pos_c 'MYCN_targets_M2919_M18532.rds' \
--neg_c 'hk_genes_normals.rds' \
--k_value 5

Rscript code/01-ruvseq-deseq.R \
--dataset 'hgg_dmg' \
--cohort_values 'PBTA' \
--cancer_group_values 'High-grade glioma/astrocytoma,Diffuse midline glioma' \
--pos_c '12915_2022_1324_MOESM4_ESM.rds' \
--neg_c 'hk_genes_normals.rds' \
--k_value 5
