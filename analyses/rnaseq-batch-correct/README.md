## Batch Correction of RNA-seq expression matrices

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi)), Adam Kraya (@aadamk)

### Description

This module uses HRT Gene Atlas housekeeping genes to correct for batch effects for tumor-only and tumor-normal comparisons. 
Batch correction is performed as follows: 
1. Input genes from the HRT Atlas v1.0 (PMID: 32663312) are assumed to be non-differentially expressed between tumor subtypes and across tumor vs normal.
Differences in housekeeping genes are assumed to be the result of technical rather than true biological variation.  
(e.g. tumors are not expected to use housekeeping genes for a growth or proliferative advantage).
2. DESeq2 analysis is run without batch correction as a baseline. 
3. [RUVg](https://bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf), which uses negative control genes for correcting batch effects, is run using multiple possible factors (`k`)
4. Diagnostic plots showing p-value distributions and clustering before and after batch correction are generated. 
5. The module then selects the `k` with the maximum sensitivity/specificity balance by considering changes in p-values
in positive and negative control gene sets.
6. DESeq2 results of the optimal method from 5 are provided and normalized counts are generated
as standard output for use in downstream visualizations. 

### Analysis scripts

#### Run tumor-only RUVg analysis 

The script below runs a tumor-only analysis for a specific cancer_group, using molecular subtype
as the default design variable. 


Example run:

```sh
Rscript code/01-ruvseq-deseq.R \
--dataset 'target_nbl' \
--cohort_value 'TARGET' \
--cancer_group_value 'Neuroblastoma' \
--k_value 5
```

Output files containing normalized counts, RUVg objects, and dge results are in the `output/[dataset]/deseq2-analysis/` folder

```sh
output/[dataset]/deseq2-analysis
├──*dge.rds # RUVg output object
└── *dge.tsv # DESeq2 DGE results
└── *pvalues.tsv # p-value results
└── *output.rds # normalized counts at all k's
```

#### Initial QC plots using PCA and UMAP

Here the goal is to create QC clustering plots using PCA and UMAP with the most variable genes (> 90% quantile) as well as housekeeping genes.

p-value histograms, PCA, and UMAP plots (.pdf):

```sh
# PCA, UMAP
plots/[dataset]/
├── *histogram.pdf # p-value histogram at all k's and for optimal k
├── *clustering.pdf # PCA and UMAP clustering results for all k and for optimal k
├── *controls.pdf # p-value histogram of negative control genes
└── *transcriptome.pdf # p-value histogram of full transcriptome
```
#### Run tumor-normal RUVg analysis 

The script below runs a tumor-normal analysis for cancer groups (comma-separated) versus gtex subgroups (comma-separated)
creating a tumor-normal comparison in the design. A key parameter is `drop`, which is a user-specifiable
parameter which drops the initial n factors of variation to guard against removing true biological variation
across tumor and normal.  


Example run:

```sh
Rscript code/02-ruvseq-deseq-tn.R \
--dataset 'HGG-DMG_TN' \
--cohort_values 'PBTA,GTEx' \
--cancer_group_values 'Diffuse midline glioma,High-grade glioma/astrocytoma' \
--gtex_subgroups 'Brain - Cortex,Brain - Cerebellum' \
--k_value 5 \
--drop 2
```

Output files containing normalized counts, RUVg objects, and dge results are in the `output/[dataset]/deseq2-analysis/` folder

```sh
output/[dataset]/deseq2-analysis
├──*dge.rds # RUVg output object
└── *dge.tsv # DESeq2 DGE results
└── *pvalues.tsv # p-value results
└── *output.rds # normalized counts at all k's
```

#### Initial QC plots using PCA and UMAP

Here the goal is to create QC clustering plots using PCA and UMAP with the most variable genes (> 90% quantile) as well as housekeeping genes.

p-value histograms, PCA, and UMAP plots (.pdf):

```sh
# PCA, UMAP
plots/[dataset]/
├── *histogram.pdf # p-value histogram at all k's and for optimal k
├── *clustering.pdf # PCA and UMAP clustering results for all k and for optimal k
├── *controls.pdf # p-value histogram of negative control genes
└── *transcriptome.pdf # p-value histogram of full transcriptome
```


#### Select optimal results 

This script evaluates DESeq2-RUVg results by pulling in all analysis results across all k's, 
evaluating changes in p-value for negative control genes (HRT Atlas genes) and positive control genes
that are analysis-specific. The script looks for k's that maximize sensitivity and specificity by
selecting the k that resulted in the most p-values increasing for negative control genes and decreasing
for positive control genes when correcting for batch effects. 


Example run tumor-normal:

```sh
Rscript code/03-ruvseq-summarization.R \
--dataset 'HGG-DMG_TN' \
--cohort_values 'PBTA,GTEx' \
--cancer_group_values 'Diffuse midline glioma,High-grade glioma/astrocytoma' \
--ruvg_dir 'deseq2_analysis' \
--pos_c 'reactome_cell_cycle_msigdb_v7.5.1.rds' \
--neg_c 'hk_genes_normals.rds' \
--analysis_type 'tumor-normal' \
--gtex_subgroups 'Brain - Cortex,Brain - Cerebellum'
```

Example run tumor-only:

```sh
Rscript code/03-ruvseq-summarization.R \
--dataset 'hgg_dmg' \
--cohort_values 'PBTA' \
--cancer_group_values 'High-grade glioma/astrocytoma,Diffuse midline glioma' \
--ruvg_dir 'deseq2_analysis' \
--pos_c '12915_2022_1324_MOESM4_ESM.rds' \
--neg_c 'hk_genes_normals.rds'
--analysis_type 'tumor-only'
```

Output files containing normalized counts, RUVg objects, and dge results are in the `output/[dataset]/deseq2-analysis/` folder

```sh
output/[dataset]/deseq2-analysis
├── /normalized_counts/*normalized_counts.rds # RUVg output object
```

Plots corresponding to the optimal `k` are as follows: 

```sh
# PCA, UMAP
plots/[dataset]/
├── *optimal_histogram.pdf # p-value histogram  and for optimal k
├── final_clustering*.pdf # PCA and UMAP clustering results for optimal k
```

###Archive
#### RUVSeq analysis 

The goal of this section is to identify and remove factors of unwanted variation. We evaluate the following methods: 

1. `RUVg`: Estimating the factors of unwanted variation using control genes. RUVg was tested with both DESeq2 as well as edgeR. Three sets of control genes were evaluated:
a. Housekeeping genes identified in Normal tissues by the HRT protocol
b. Housekeeping genes identified collectively in Tumors + Normal tissues using the HRT protocol
c. Empirical control genes from the first-pass DGE analysis. These are “in-silico empirical” negative controls, e.g., least significantly DE genes (p-adjust > 0.05) based on a first-pass DE analysis performed prior to RUVg normalization.
d. Housekeeping genes identified in Normal tissues by the HRT protocol that are differentially expressed with p-adjust < 0.05 from the first-pass DGE analysis.

2. `RUVr`: use the residuals from a DGE analysis (without batch correction) as a source of unwanted variation to correct for in a second pass DGE. RUVr was tested only using edgeR.

For each of these methods, multiple Ks are evaluated from 1 to 5. A p-value histogram is created for the full transcriptome. A chisq test as well as a KS test is performed on the p-values of negative control genes obtained after second pass DEG analysis to identify if the distribution of the p-values is significantly different than uniform distribution.   

Example run:

```sh
# matched samples in PBTA cohort
# 4 CNS Embryonal tumor, 4 Diffuse midline glioma and 10 High-grade glioma/astrocytoma
# sample size is okay in HGG so we will use that as a test
Rscript code/03-ruvseq-matched-samples.R \
--dataset match_pbta_hgg \
--k_value 5

# matched samples in TARGET cohort
# 24 Acute Lymphoblastic Leukemia samples have been sequenced with > 1 technology
Rscript code/03-ruvseq-matched-samples.R \
--dataset match_target_all \
--k_value 5

# full run script
bash run_ruvseq_matched_samples.sh
```

General description of output files:

* RUVg analysis:

1. `{dge_empirical_genes, hk_genes_normals, hk_genes_tumor_normals, dge_empirical_hk_genes}_dge_ruvg_{edger, deseq2}_{chisq, ks}_pvalues.tsv`: Chisq/KS test of uniformity on p-values obtained using second pass DEG analysis for all Ks (only negative control gene set)

2. `{dge_empirical_genes, hk_genes_normals, hk_genes_tumor_normals, dge_empirical_hk_genes}_dge_ruvg_{edger, deseq2}_clustering.pdf`: PCA and UMAP after RUVg for all Ks. (full transcriptome)

3. `{dge_empirical_genes, hk_genes_normals, hk_genes_tumor_normals, dge_empirical_hk_genes}_dge_ruvg_{edger, deseq2}_histogram_{controls, full_transcriptome}.pdf`: Histogram of p-values for second-pass DESeq2 analysis for all Ks using negative control genes and full transcriptome.

* RUVr analysis:

1. `dge_ruvr_edger_{chisq, ks}_pvalues.tsv`: Chisq/KS test of uniformity on p-values obtained using second pass DEG analysis for all Ks (only negative control gene set)

2. `dge_ruvr_edger_clustering.pdf`: PCA and UMAP after RUVr for all Ks. (full transcriptome)

3. `dge_ruvr_edger_histogram_{controls, full_transcriptome}.pdf`: Histogram of p-values for second-pass DESeq2 analysis for all Ks using negative control genes and full transcriptome.

Output files:

1. Matched PBTA HGG samples:

```sh
## PCA and UMAP clustering as well as boxplots before running RUVg or RUVr
output/match_pbta_hgg/
├── boxplots_with_and_without_norm.pdf 
└── clustering_with_and_without_norm.pdf

## DESeq2 + RUVg analysis
# histogram of p-values in first pass DESeq2-based DEG analysis
output/match_pbta_hgg/deseq2_analysis
└── dge_deseq2_histogram.pdf 

# using empirical genes identified in first pass DESeq2-based DEG analysis (p-adjust > 0.05)
output/match_pbta_hgg/deseq2_analysis
├── dge_empirical_genes_dge_ruvg_deseq2_chisq_pvalues.tsv 
├── dge_empirical_genes_dge_ruvg_deseq2_ks_pvalues.tsv 
├── dge_empirical_genes_dge_ruvg_deseq2_clustering.pdf 
├── dge_empirical_genes_dge_ruvg_deseq2_histogram_controls.pdf 
└── dge_empirical_genes_dge_ruvg_deseq2_histogram_full_transcriptome.pdf 

# using housekeeping genes in normals (HRT atlas)
output/match_pbta_hgg/deseq2_analysis
├── hk_genes_normals_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_normals_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_normals_dge_ruvg_deseq2_clustering.pdf
├── hk_genes_normals_dge_ruvg_deseq2_histogram_controls.pdf
└── hk_genes_normals_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_pbta_hgg/deseq2_analysis
├── hk_genes_tumor_normals_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_deseq2_clustering.pdf
├── hk_genes_tumor_normals_dge_ruvg_deseq2_histogram_controls.pdf
└── hk_genes_tumor_normals_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass DESeq2-based DGE analysis (p-adjust < 0.05)
output/match_pbta_hgg/deseq2_analysis
├── dge_empirical_hk_genes_dge_ruvg_deseq2_chisq_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_deseq2_ks_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_deseq2_clustering.pdf
├── dge_empirical_hk_genes_dge_ruvg_deseq2_histogram_controls.pdf
└── dge_empirical_hk_genes_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

## edgeR + RUVg analysis
# histogram of p-values in first-pass edgeR-based DEG analysis
output/match_pbta_hgg/edger_analysis
└── dge_edger_histogram.pdf

# using empirical genes identified in first pass edgeR-based DEG analysis (p-adjust > 0.05)
output/match_pbta_hgg/edger_analysis
├── dge_empirical_genes_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_genes_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_genes_dge_ruvg_edger_clustering.pdf
├── dge_empirical_genes_dge_ruvg_edger_histogram_controls.pdf
└── dge_empirical_genes_dge_ruvg_edger_histogram_full_transcriptome.pdf

# using housekeeping genes in normals (HRT atlas)
output/match_pbta_hgg/edger_analysis
├── hk_genes_normals_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_normals_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_normals_dge_ruvg_edger_clustering.pdf
├── hk_genes_normals_dge_ruvg_edger_histogram_controls.pdf
└── hk_genes_normals_dge_ruvg_edger_histogram_full_transcriptome.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_pbta_hgg/edger_analysis
├── hk_genes_tumor_normals_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_edger_clustering.pdf
├── hk_genes_tumor_normals_dge_ruvg_edger_histogram_controls.pdf
└── hk_genes_tumor_normals_dge_ruvg_edger_histogram_full_transcriptome.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass edgeR-based DGE analysis (p-adjust < 0.05)
output/match_pbta_hgg/edger_analysis
├── dge_empirical_hk_genes_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_edger_clustering.pdf
├── dge_empirical_hk_genes_dge_ruvg_edger_histogram_controls.pdf
└── dge_empirical_hk_genes_dge_ruvg_edger_histogram_full_transcriptome.pdf

## edgeR + RUVr analysis
# using residuals identified in first pass edgeR-based DEG analysis
output/match_pbta_hgg/edger_analysis
├── dge_ruvr_edger_chisq_pvalues.tsv
├── dge_ruvr_edger_ks_pvalues.tsv
├── dge_ruvr_edger_clustering.pdf
├── dge_ruvr_edger_histogram_controls.pdf
└── dge_ruvr_edger_histogram_full_transcriptome.pdf
```

2. Matched TARGET ALL samples:

```sh
## PCA and UMAP clustering as well as boxplots before running RUVg or RUVr
output/match_target_all
├── boxplots_with_and_without_norm.pdf 
└── clustering_with_and_without_norm.pdf

## DESeq2 + RUVg analysis
# histogram of p-values in first pass DESeq2-based DEG analysis
output/match_target_all/deseq2_analysis
└── dge_deseq2_histogram.pdf

# using empirical genes identified in first pass DESeq2-based DEG analysis (p-adjust > 0.05)
output/match_target_all/deseq2_analysis
├── dge_empirical_genes_dge_ruvg_deseq2_chisq_pvalues.tsv
├── dge_empirical_genes_dge_ruvg_deseq2_ks_pvalues.tsv
├── dge_empirical_genes_dge_ruvg_deseq2_clustering.pdf
├── dge_empirical_genes_dge_ruvg_deseq2_histogram_controls.pdf
└── dge_empirical_genes_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

# using housekeeping genes in normals (HRT atlas)
output/match_target_all/deseq2_analysis
├── hk_genes_normals_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_normals_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_normals_dge_ruvg_deseq2_clustering.pdf
├── hk_genes_normals_dge_ruvg_deseq2_histogram_controls.pdf
└── hk_genes_normals_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_target_all/deseq2_analysis
├── hk_genes_tumor_normals_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_deseq2_clustering.pdf
├── hk_genes_tumor_normals_dge_ruvg_deseq2_histogram_controls.pdf
└── hk_genes_tumor_normals_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass DESeq2-based DGE analysis (p-adjust < 0.05)
output/match_target_all/deseq2_analysis
├── dge_empirical_hk_genes_dge_ruvg_deseq2_chisq_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_deseq2_ks_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_deseq2_clustering.pdf
├── dge_empirical_hk_genes_dge_ruvg_deseq2_histogram_controls.pdf
└── dge_empirical_hk_genes_dge_ruvg_deseq2_histogram_full_transcriptome.pdf

## edgeR + RUVg analysis
# histogram of p-values in first-pass edgeR-based DEG analysis
output/match_target_all/edger_analysis
└── dge_edger_histogram.pdf

# using empirical genes identified in first pass edgeR-based DEG analysis (p-adjust > 0.05)
output/match_target_all/edger_analysis
├── dge_empirical_genes_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_genes_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_genes_dge_ruvg_edger_clustering.pdf
├── dge_empirical_genes_dge_ruvg_edger_histogram_controls.pdf
└── dge_empirical_genes_dge_ruvg_edger_histogram_full_transcriptome.pdf

# using housekeeping genes in normals (HRT atlas)
output/match_target_all/edger_analysis
├── hk_genes_normals_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_normals_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_normals_dge_ruvg_edger_clustering.pdf
├── hk_genes_normals_dge_ruvg_edger_histogram_controls.pdf
└── hk_genes_normals_dge_ruvg_edger_histogram_full_transcriptome.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_target_all/edger_analysis
├── hk_genes_tumor_normals_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_tumor_normals_dge_ruvg_edger_clustering.pdf
├── hk_genes_tumor_normals_dge_ruvg_edger_histogram_controls.pdf
└── hk_genes_tumor_normals_dge_ruvg_edger_histogram_full_transcriptome.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass edgeR-based DGE analysis (p-adjust < 0.05)
output/match_target_all/edger_analysis
├── dge_empirical_hk_genes_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_hk_genes_dge_ruvg_edger_clustering.pdf
├── dge_empirical_hk_genes_dge_ruvg_edger_histogram_controls.pdf
└── dge_empirical_hk_genes_dge_ruvg_edger_histogram_full_transcriptome.pdf

## edgeR + RUVr analysis
# using residuals identified in first pass edgeR-based DEG analysis
output/match_target_all/edger_analysis
├── dge_ruvr_edger_chisq_pvalues.tsv
├── dge_ruvr_edger_ks_pvalues.tsv
├── dge_ruvr_edger_clustering.pdf
├── dge_ruvr_edger_histogram_controls.pdf
└── dge_ruvr_edger_histogram_full_transcriptome.pdf
```

#### ComBat and ComBat_seq analysis 

In this section, we will generate some QC plots with and without `sva::ComBat` and `sva::ComBat_seq` batch correction using `RNA_library` as the batch variable.

Example run:

* ComBat

```sh
# matched samples in PBTA cohort
# 4 CNS Embryonal tumor, 4 Diffuse midline glioma and 10 High-grade glioma/astrocytoma
# sample size is okay in HGG so we will use that as a test
Rscript code/04-test_combat.R \
--dataset match_pbta_hgg 

# matched samples in TARGET cohort
# 24 Acute Lymphoblastic Leukemia samples have been sequenced with > 1 technology
Rscript code/04-test_combat.R \
--dataset match_target_all

# full run script
bash run_combat.sh
```

* ComBat_seq

```sh
# matched samples in PBTA cohort
# 4 CNS Embryonal tumor, 4 Diffuse midline glioma and 10 High-grade glioma/astrocytoma
# sample size is okay in HGG so we will use that as a test
Rscript code/04-test_combat_seq.R \
--dataset match_pbta_hgg 

# matched samples in TARGET cohort
# 24 Acute Lymphoblastic Leukemia samples have been sequenced with > 1 technology
Rscript code/04-test_combat_seq.R \
--dataset match_target_all 

# full run script
bash run_combat_seq.sh
```

Output files:

1. `{expressed_genes, hk_genes_normals, hk_genes_tumors_normals}_boxplots_with_and_without_bc.pdf`: Boxplots of matched samples using either all expressed genes, housekeeping genes in Normals from HRT, housekeeping genes in Tumors + Normals using HRT protocol. 

2. `{expressed_genes, hk_genes_normals, hk_genes_tumors_normals}_clustering_with_and_without_bc.pdf`: PCA and UMAP of matched samples using either all expressed genes, housekeeping genes in Normals from HRT, housekeeping genes in Tumors + Normals using HRT protocol.

* ComBat outputs

```sh
# matched samples from PBTA HGG
output/combat_output/match_pbta_hgg
├── expressed_genes_boxplots_with_and_without_bc.pdf 
├── expressed_genes_clustering_with_and_without_bc.pdf
├── hk_genes_normals_boxplots_with_and_without_bc.pdf
├── hk_genes_normals_clustering_with_and_without_bc.pdf
├── hk_genes_tumors_normals_boxplots_with_and_without_bc.pdf
└── hk_genes_tumors_normals_clustering_with_and_without_bc.pdf

# matched samples from TARGET ALL
output/combat_output/match_target_all
├── expressed_genes_boxplots_with_and_without_bc.pdf
├── expressed_genes_clustering_with_and_without_bc.pdf
├── hk_genes_normals_boxplots_with_and_without_bc.pdf
├── hk_genes_normals_clustering_with_and_without_bc.pdf
├── hk_genes_tumors_normals_boxplots_with_and_without_bc.pdf
└── hk_genes_tumors_normals_clustering_with_and_without_bc.pdf
```

* ComBat_seq outputs

```sh
# matched sampless from PBTA HGG
output/combatseq_output/match_pbta_hgg
├── expressed_genes_boxplots_with_and_without_bc.pdf
├── expressed_genes_clustering_with_and_without_bc.pdf
├── hk_genes_normals_boxplots_with_and_without_bc.pdf
├── hk_genes_normals_clustering_with_and_without_bc.pdf
├── hk_genes_tumors_normals_boxplots_with_and_without_bc.pdf
└── hk_genes_tumors_normals_clustering_with_and_without_bc.pdf

# matched samples from TARGET ALL
output/combatseq_output/match_target_all
├── expressed_genes_boxplots_with_and_without_bc.pdf
├── expressed_genes_clustering_with_and_without_bc.pdf
├── hk_genes_normals_boxplots_with_and_without_bc.pdf
├── hk_genes_normals_clustering_with_and_without_bc.pdf
├── hk_genes_tumors_normals_boxplots_with_and_without_bc.pdf
└── hk_genes_tumors_normals_clustering_with_and_without_bc.pdf
```
