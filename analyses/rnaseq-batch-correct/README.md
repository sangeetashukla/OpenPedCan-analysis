## Batch Correction of RNA-seq expression matrices

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to identify and correct the source of technical variation in RNA-seq datasets. 

### Analysis scripts

#### Identify housekeeping genes

The goal of this section is to generate a list of house keeping genes identified collectively in tumors as well as normal samples using the HRT Atlas protocol: https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv 


Example run:

```sh
Rscript 00-tumor_specific_hk_genes.R
```

Output files in .rds format are generated under `input` folder:

```sh
input
├── hk_genes_normals.rds # downloaded from HRT atlas
└── hk_genes_tumor_normals.rds # output of 00-tumor_specific_hk_genes.R
```

#### Initial QC plots using PCA and UMAP

Here the goal is to create initial QC clustering plots using PCA and UMAP with the most variable genes (> 90% quantile) as well as housekeeping genes. We are using three datasets for testing:

1. Tumor + Normal samples
2. Tumor samples only
3. Matched tumor samples i.e. patients with samples sequenced with different `RNA_library` protocol.

Example run:

```sh
# compute PCA and UMAP from count data
Rscript 01-full_dataset_compute_pca_umap_counts.R

# create PCA plots 
Rscript 02-full_dataset_plot_pca_counts.R

# create UMAP plots
Rscript 02-full_dataset_plot_umap_counts.R

# full run script
bash run_pca_umap_clustering.sh
```

Output files of PCA and UMAP output (.rds) as well as plots (.pdf):

```sh
# PCA
# tumor + normal samples across all cohorts
output/QC_clustering
├── pca_output_tumors_normals.rds # most variable genes > 90% quantile 
├── pca_output_tumors_normals_hk.rds # housekeeping genes
├── pca_output_tumors_normals_composition_combined.pdf # combined plot with composition type
└── pca_output_tumors_normals_library_combined.pdf # combined plot with RNA library type

# tumor samples across all cohorts
output/QC_clustering
├── pca_output_tumors.rds # most variable genes > 90% quantile 
├── pca_output_tumors_hk.rds # housekeeping genes
├── pca_output_tumors_composition_combined.pdf # combined plot with composition type 
└── pca_output_tumors_library_combined.pdf # combined plot with RNA library type

# matched samples from PBTA and TARGET
output/QC_clustering
├── pca_output_matched_samples.rds # most variable genes > 90% quantile
├── pca_output_matched_samples_hk.rds # housekeeping genes
├── pca_output_matched_samples_composition_combined.pdf # combined plot with composition type 
└── pca_output_matched_samples_library_combined.pdf # combined plot with RNA library type

# UMAP
# tumor + normal samples across all cohorts
output/QC_clustering
├── umap_output_tumors_normals.rds
├── umap_output_tumors_normals_hk.rds
├── umap_output_tumors_normals_composition_combined.pdf
└── umap_output_tumors_normals_library_combined.pdf

# tumor samples across all cohorts
output/QC_clustering
├── umap_output_tumors.rds
├── umap_output_tumors_hk.rds
├── umap_output_tumors_composition_combined.pdf
└── umap_output_tumors_library_combined.pdf

# matched samples from PBTA and TARGET
output/QC_clustering
├── umap_output_matched_samples.rds
├── umap_output_matched_samples_hk.rds
├── umap_output_matched_samples_composition_combined.pdf
└── umap_output_matched_samples_library_combined.pdf
```

#### RUVSeq analysis 

The goal of this section is to identify and remove factors of unwanted variation. We evaluate the following methods: 

1. `RUVg`: Estimating the factors of unwanted variation using control genes. RUVg was tested with both DESeq2 as well as edgeR. Three sets of control genes were evaluated:
a. Housekeeping genes identified in Normal tissues by the HRT protocol
b. Housekeeping genes identified collectively in Tumors + Normal tissues using the HRT protocol
c. Empirical control genes from the first-pass DGE analysis. These are “in-silico empirical” negative controls, e.g., least significantly DE genes (p-adjust > 0.05) based on a first-pass DE analysis performed prior to RUVg normalization.
d. Housekeeping genes identified in Normal tissues by the HRT protocol that are differentially expressed with p-adjust < 0.05 from the first-pass DGE analysis.

2. `RUVr`: use the residuals from a DGE analysis (without batch correction) as a source of unwanted variation to correct for in a second pass DGE. RUVr was tested only using edgeR.

For each of these methods, multiple Ks are evaluated from 1 to 5. A chisq test as well as a KS test is performed on the p-values obtained using second pass DEG analysis to identify if the distribution of the p-values is significantly different than uniform distribution.   

Example run:

```sh
# matched samples in PBTA cohort
# 4 CNS Embryonal tumor, 4 Diffuse midline glioma and 10 High-grade glioma/astrocytoma
# sample size is okay in HGG so we will use that as a test
Rscript code/03-test_ruvseq.R \
--dataset match_pbta_hgg \
--k_value 5

# matched samples in TARGET cohort
# 24 Acute Lymphoblastic Leukemia samples have been sequenced with > 1 technology
Rscript code/03-test_ruvseq.R \
--dataset match_target_all \
--k_value 5

# full run script
bash run_ruvseq.sh
```

General description of output files:

* RUVg analysis:

1. `{dge_empirical_genes, hk_genes_normals, hk_genes_tumor_normals, dge_empirical_hk_genes}_stranded_vs_polya_dge_ruvg_{edger, deseq2}_{chisq, ks}_pvalues.tsv`: Chisq/KS test of uniformity on p-values obtained using second pass DEG analysis for all Ks. 

2. `{dge_empirical_genes, hk_genes_normals, hk_genes_tumor_normals, dge_empirical_hk_genes}_stranded_vs_polya_dge_ruvg_{edger, deseq2}_clustering.pdf`: PCA and UMAP after RUVg for all Ks

3. `{dge_empirical_genes, hk_genes_normals, hk_genes_tumor_normals, dge_empirical_hk_genes}_stranded_vs_polya_dge_ruvg_{edger, deseq2}_histogram.pdf`: Histogram of p-values for second-pass DESeq2 analysis for all Ks

* RUVr analysis:

1. `stranded_vs_polya_dge_ruvr_edger_{chisq, ks}_pvalues.tsv`: Chisq/KS test of uniformity on p-values obtained using second pass DEG analysis for all Ks. 

2. `stranded_vs_polya_dge_ruvr_edger_clustering.pdf`: PCA and UMAP after applying RUVg for all Ks.

3. `stranded_vs_polya_dge_ruvr_edger_histogram.pdf`: Histogram of p-values for second-pass DESeq2 analysis for all Ks.

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
└── stranded_vs_polya_dge_deseq2_histogram.pdf 

# using empirical genes identified in first pass DESeq2-based DEG analysis (p-adjust > 0.05)
output/match_pbta_hgg/deseq2_analysis
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv 
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv 
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf 
└── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf 

# using housekeeping genes in normals (HRT atlas)
output/match_pbta_hgg/deseq2_analysis
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_pbta_hgg/deseq2_analysis
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass DESeq2-based DGE analysis (p-adjust < 0.05)
output/match_pbta_hgg/deseq2_analysis
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

## edgeR + RUVg analysis
# histogram of p-values in first-pass edgeR-based DEG analysis
output/match_pbta_hgg/edger_analysis
└── stranded_vs_polya_dge_edger_histogram.pdf

# using empirical genes identified in first pass edgeR-based DEG analysis (p-adjust > 0.05)
output/match_pbta_hgg/edger_analysis
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

# using housekeeping genes in normals (HRT atlas)
output/match_pbta_hgg/edger_analysis
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_pbta_hgg/edger_analysis
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass edgeR-based DGE analysis (p-adjust < 0.05)
output/match_pbta_hgg/edger_analysis
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

## edgeR + RUVr analysis
# using residuals identified in first pass edgeR-based DEG analysis
output/match_pbta_hgg/edger_analysis
├── stranded_vs_polya_dge_ruvr_edger_chisq_pvalues.tsv
├── stranded_vs_polya_dge_ruvr_edger_ks_pvalues.tsv
├── stranded_vs_polya_dge_ruvr_edger_clustering.pdf
└── stranded_vs_polya_dge_ruvr_edger_histogram.pdf
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
└── stranded_vs_polya_dge_deseq2_histogram.pdf

# using empirical genes identified in first pass DESeq2-based DEG analysis (p-adjust > 0.05)
output/match_target_all/deseq2_analysis
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── dge_empirical_genes_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

# using housekeeping genes in normals (HRT atlas)
output/match_target_all/deseq2_analysis
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── hk_genes_normals_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_target_all/deseq2_analysis
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass DESeq2-based DGE analysis (p-adjust < 0.05)
output/match_target_all/deseq2_analysis
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_chisq_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_ks_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_clustering.pdf
└── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_deseq2_histogram.pdf

## edgeR + RUVg analysis
# histogram of p-values in first-pass edgeR-based DEG analysis
output/match_target_all/edger_analysis
└── stranded_vs_polya_dge_edger_histogram.pdf

# using empirical genes identified in first pass edgeR-based DEG analysis (p-adjust > 0.05)
output/match_target_all/edger_analysis
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── dge_empirical_genes_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

# using housekeeping genes in normals (HRT atlas)
output/match_target_all/edger_analysis
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── hk_genes_normals_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

# using housekeeping genes identified in tumors + normals (following HRT protocol)
output/match_target_all/edger_analysis
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── hk_genes_tumor_normals_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

# differentially expressed housekeeping genes in normals (HRT atlas) from the first-pass edgeR-based DGE analysis (p-adjust < 0.05)
output/match_target_all/edger_analysis
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_chisq_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_ks_pvalues.tsv
├── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_clustering.pdf
└── dge_empirical_hk_genes_stranded_vs_polya_dge_ruvg_edger_histogram.pdf

## edgeR + RUVr analysis
# using residuals identified in first pass edgeR-based DEG analysis
output/match_target_all/edger_analysis
├── stranded_vs_polya_dge_ruvr_edger_chisq_pvalues.tsv
├── stranded_vs_polya_dge_ruvr_edger_ks_pvalues.tsv
├── stranded_vs_polya_dge_ruvr_edger_clustering.pdf
└── stranded_vs_polya_dge_ruvr_edger_histogram.pdf
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
