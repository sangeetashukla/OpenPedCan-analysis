# script to compute PCA and UMAP from log2-normalized count data

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(uwot)
})

# output directory
output_dir <- file.path('output', 'QC_clustering')
dir.create(output_dir, showWarnings = F, recursive = T)

# histology file
hist_file <- read.delim('../../data/histologies.tsv')
hist_file <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         !cohort %in% "TCGA")

# housekeeping genes
hk_genes_tumors_normals <- readRDS('input/hk_genes_tumor_normals.rds')

# read expression data
dat <- readRDS('../../data/gene-counts-rsem-expected_count-collapsed.rds')
dat <- dat[rowSums(dat) > 0,] # remove all genes with zero counts across all samples
log_expression <- log2(dat + 1)

# 1. TUMOR + NORMALS
# 1.1 filter to high variance genes above 90% quantile (3666 x 19995)
row_variances <- matrixStats::rowVars(as.matrix(log_expression))
high_var_exp <- log_expression[which(row_variances > quantile(row_variances, 0.9)), ]

# pca
fname <- file.path(output_dir, 'pca_output_tumors_normals.rds')
if(!file.exists(fname)){
  pca_results <- prcomp(high_var_exp)
  saveRDS(pca_results, file = fname)
}

# umap
fname <- file.path(output_dir, 'umap_output_tumors_normals.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(high_var_exp), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(high_var_exp)
  saveRDS(umap_out, file = fname)
}

# 1.2. only house keeping genes found in tumor + normals (HRT protocol) (25 x 19995)
hk_genes <- log_expression[rownames(log_expression) %in% hk_genes_tumors_normals,] 

# pca
fname <- file.path(output_dir, 'pca_output_tumors_normals_hk.rds')
if(!file.exists(fname)){
  pca_results_hk <- prcomp(hk_genes)
  saveRDS(pca_results_hk, file = fname)
}

# umap
fname <- file.path(output_dir, 'umap_output_tumors_normals_hk.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(hk_genes), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(hk_genes)
  saveRDS(umap_out, file = fname)
}

# 2. TUMORS only
tumors <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         sample_type == "Tumor",
         !cohort %in% c("GTEx", "TCGA"))
log_expression <- log_expression %>%
  dplyr::select(tumors$Kids_First_Biospecimen_ID)

# 2.1. filter to high variance genes above 90% quantile (3666 x 2601)
row_variances <- matrixStats::rowVars(as.matrix(log_expression))
high_var_exp <- log_expression[which(row_variances > quantile(row_variances, 0.9)), ]

# pca
fname <- file.path(output_dir, 'pca_output_tumors.rds')
if(!file.exists(fname)){
  pca_results <- prcomp(high_var_exp)
  saveRDS(pca_results, file = fname)
}

# umap
fname <- file.path(output_dir, 'umap_output_tumors.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(high_var_exp), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(high_var_exp)
  saveRDS(umap_out, file = fname)
}

# 2.2. only house keeping genes (25 x 2601)
hk_genes <- log_expression[rownames(log_expression) %in% hk_genes_tumors_normals,]

# pca
fname <- file.path(output_dir, 'pca_output_tumors_hk.rds')
if(!file.exists(fname)){
  pca_results_hk <- prcomp(hk_genes)
  saveRDS(pca_results_hk, file = fname)
}

# umap
fname <- file.path(output_dir, 'umap_output_tumors_hk.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(hk_genes), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(hk_genes)
  saveRDS(umap_out, file = fname)
}

# 3. TUMORS matched samples (keep only matched participant + sample ids)
matched_samples_hist <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         sample_type == "Tumor") %>%
  group_by(Kids_First_Participant_ID, sample_id) %>%
  mutate(n = length(sort(unique(RNA_library)))) %>% 
  filter(n > 1)

# 3.1. filter to high variance genes above 90% quantile (3666 x 42)
log_expression_matched_samples <- log_expression %>%
  dplyr::select(matched_samples_hist$Kids_First_Biospecimen_ID)
row_variances <- matrixStats::rowVars(as.matrix(log_expression_matched_samples))
matched_samples_high_var <- log_expression_matched_samples[which(row_variances > quantile(row_variances, 0.9)), ]

# pca
fname <- file.path(output_dir, 'pca_output_matched_samples.rds')
if(!file.exists(fname)){
  pca_results <- prcomp(matched_samples_high_var)
  saveRDS(pca_results, file = fname)
}

# umap
fname <- file.path(output_dir, 'umap_output_matched_samples.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(matched_samples_high_var), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(matched_samples_high_var)
  saveRDS(umap_out, file = fname)
}

# 3.2. only house keeping genes (25 x 42)
matched_samples_hk <- log_expression_matched_samples[rownames(log_expression_matched_samples) %in% hk_genes_tumors_normals,]

# pca
fname <- file.path(output_dir, 'pca_output_matched_samples_hk.rds')
if(!file.exists(fname)){
  pca_results <- prcomp(matched_samples_hk)
  saveRDS(pca_results, file = fname)
}

# umap
fname <- file.path(output_dir, 'umap_output_matched_samples_hk.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(matched_samples_hk), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(matched_samples_hk)
  saveRDS(umap_out, file = fname)
}
