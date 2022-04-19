# This script performs the following functions:
# 1. DESeq2 analysis with RUVg 
# 2. edgeR analysis with RUVg and RUVr

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(RUVSeq)
})

# parse parameters
option_list <- list(
  make_option(opt_str = "--dataset", type = "character",
              help = "Dataset for running differential gene expression analysis: match_pbta_hgg or match_target_all"),
  make_option(opt_str = "--k_value", type = "character",
              help = "number of K values to evaluate")
)
opt <- parse_args(OptionParser(option_list = option_list))
dataset <- opt$dataset
k_value <- as.numeric(opt$k_value)
output_dir <- file.path('output', dataset)
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions 
source('util/deseq2_pvals_histogram.R') # DESeq2 pval histograms
source('util/edaseq_plot.R') # PCA and UMAP clustering
source('util/box_plots.R') # boxplots for samples
source('util/ruvg_test.R') # function to run RUVg
source('util/ruvr_test.R') # function to run RUVr (only works with edgeR)

# read histology
htl_df <- readr::read_tsv('../../data/histologies.tsv')

# read expected counts
cnt_df <- readRDS('../../data/gene-counts-rsem-expected_count-collapsed.rds')

# filter histology
if (dataset == "match_target_all") {
  # 24 Acute Lymphoblastic Leukemia samples have been sequenced with > 1 technology
  selected_htl_df <- htl_df %>%
    filter(experimental_strategy == "RNA-Seq",
           sample_type == "Tumor",
           cohort == "TARGET") %>%
    group_by(Kids_First_Participant_ID, sample_id) %>%
    mutate(n = length(sort(unique(RNA_library)))) %>% 
    filter(n > 1)
} else if (dataset == "match_pbta_hgg") {
  # 4 CNS Embryonal tumor, 4 Diffuse midline glioma and 10 High-grade glioma/astrocytoma
  # sample size is okay in HGG so we will use that as a test
  selected_htl_df <- htl_df %>%
    filter(experimental_strategy == "RNA-Seq",
           sample_type == "Tumor",
           cohort == "PBTA",
           cancer_group == "High-grade glioma/astrocytoma") %>%
    group_by(Kids_First_Participant_ID, sample_id) %>%
    mutate(n = length(sort(unique(RNA_library)))) %>% 
    filter(n > 1)
} else {
  stop(paste0('unknown dataset', dataset))
}

# filter expression
cnt_df <- cnt_df %>%
  dplyr::select(selected_htl_df$Kids_First_Biospecimen_ID)
cnt_df <- cnt_df[rowSums(cnt_df) > 0,]

# filter by expression and convert to DGElist
# for paired analysis, we need to consider patient id (e.g. https://www.biostars.org/p/437279/)
bs_id <- gsub('TARGET-[0-9]{2}-', '', selected_htl_df$Kids_First_Biospecimen_ID) # shorten bs ids for TARGET samples
patient_id <- factor(as.character(selected_htl_df$Kids_First_Participant_ID))
rna_library <- factor(as.character(selected_htl_df$RNA_library))
design <- model.matrix(~ 0 + patient_id + rna_library)
counts_object <- edgeR::DGEList(counts = cnt_df, group = rna_library)
counts_object_filtered <- edgeR::filterByExpr(counts_object, group = rna_library)
counts_object <- counts_object[counts_object_filtered, , keep.lib.sizes = FALSE]

# create new expression set
seq_expr_set <- newSeqExpressionSet(counts = round(counts_object$counts), phenoData = data.frame(patient_id, rna_library, bs_id, row.names = colnames(round(counts_object$counts))))
pca_p <- edaseq_plot(object = seq_expr_set, title = "Before normalization", type = "PCA")
umap_p <- edaseq_plot(object = seq_expr_set, title = "Before normalization", type = "UMAP")
boxplot_p <- box_plots(object = seq_expr_set, title = "Before normalization")

# from https://support.bioconductor.org/p/95805/#95808
# The "betweenLaneNormalization" is just a poorly named function that perform between-sample normalization, 
# independently of the fact that they are organized in lanes, plates, etc. 
# It is used in RUVSeq just to make sure that the factors of unwanted variation don't capture the difference in library size that should be captured by the size factors. 
# normalize the data using upper-quartile (UQ) normalization
seq_expr_set <- EDASeq::betweenLaneNormalization(seq_expr_set, which = "upper")
pca_p_norm <- edaseq_plot(object = seq_expr_set, title = "After UQ normalization", type = "PCA")
umap_p_norm <- edaseq_plot(object = seq_expr_set, title = "After UQ normalization", type = "UMAP")
boxplot_p_norm <- box_plots(object = seq_expr_set, title = "After UQ normalization")

# save both UMAP and PCA in a single file
p <- ggpubr::ggarrange(pca_p, pca_p_norm, umap_p, umap_p_norm, common.legend = T, legend = "bottom")
fname <- file.path(output_dir, 'clustering_with_and_without_norm.pdf')
ggsave(filename = fname, plot = p, width = 6, height = 6, device = "pdf", bg = "white")

# save boxplots in a single file
p <- ggpubr::ggarrange(boxplot_p, boxplot_p_norm, nrow = 2, common.legend = T, legend = "bottom")
fname <- file.path(output_dir, 'boxplots_with_and_without_norm.pdf')
ggsave(filename = fname, plot = p, width = 12, height = 6, device = "pdf", bg = "white")

########################################### DESeq2 analysis ##############################################################

# 1.1. RUVg using hk genes in tumor + normals
emp_neg_ctrl_genes_hk_tumor_normal <- readRDS('input/hk_genes_tumor_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_hk_tumor_normal, prefix = "hk_genes_tumor_normals", 
          diff_type = "deseq2", 
          output_dir = output_dir)

# 1.2. RUVg using hk genes in normals only (full HRT atlas)
emp_neg_ctrl_genes_normals <- readRDS('input/hk_genes_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_normals, prefix = "hk_genes_normals", 
          diff_type = "deseq2", 
          output_dir = output_dir)

# 1.3. RUVg using empirical control genes (dge from DESeq2 analysis)
# From DESeq2 documentation:
# 1. DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts_object$counts), colData =  data.frame(patient_id, rna_library), design = design)
dds <- DESeq2::DESeq(dds)
dge_output <- DESeq2::results(dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
dge_output <- dge_output %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>%
  arrange(padj) 

# export results (commenting out to reduce output files)
# filename <- file.path(output_dir, 'stranded_vs_polya_dge_deseq2_result.tsv')
# write_tsv(dge_output, file = filename)

# plot and save p-value histogram
p <- deseq2_pvals_histogram(res_df = dge_output,
                            xlab = 'stranded vs poly-A RNA-seq DGE RLE nbinomWaldTest p-value',
                            ylab = 'Gene count', title = paste0('Histogram of stranded vs poly-A RNA-seq\n', 'DESeq2 nbinomWaldTest p-values'))
filename <- file.path(output_dir, 'deseq2_analysis', 'stranded_vs_polya_dge_deseq2_histogram.pdf')
ggsave(filename = filename, plot = p, width = 8, height = 7, bg = "white")

# pull empirical control genes to be used in RUVg (not differential expressed)
emp_neg_ctrl_genes_dge_deseq2 <- dge_output %>% 
  filter(padj > 0.05) %>%
  pull(gene)
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value,
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_dge_deseq2,  prefix = "dge_empirical_genes", 
          diff_type = "deseq2", 
          output_dir = output_dir)

# 1.4 HK genes (HRT atlas) differentially expressed in DEG analysis
hk_genes_normals <- readRDS('input/hk_genes_normals.rds')
emp_neg_ctrl_hk_genes_dge_deseq2 <- dge_output %>% 
  filter(padj < 0.05,
         gene %in% hk_genes_normals) %>%
  pull(gene)
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value,
          emp_neg_ctrl_genes = emp_neg_ctrl_hk_genes_dge_deseq2,  prefix = "dge_empirical_hk_genes", 
          diff_type = "deseq2", 
          output_dir = output_dir)

print("DESeq2 analysis finished.")
########################################### DESeq2 analysis ##############################################################

########################################### edgeR analysis ##############################################################

# 2.1. RUVg using hk genes in tumor + normals
emp_neg_ctrl_genes_hk_tumor_normal <- readRDS('input/hk_genes_tumor_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_hk_tumor_normal, prefix = "hk_genes_tumor_normals", 
          diff_type = "edger", 
          output_dir = output_dir)

# 2.2. RUVg using hk genes in normals only (full HRT atlas)
emp_neg_ctrl_genes_normals <- readRDS('input/hk_genes_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_normals, prefix = "hk_genes_normals", 
          diff_type = "edger", 
          output_dir = output_dir)

# 2.3 RUVg using empirical control genes (dge from edgeR)
y <- DGEList(counts = counts(seq_expr_set), group = rna_library)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
edger_fit <- glmFit(y, design)

# export results
lrt <- glmLRT(edger_fit, coef = grep('rna_library', colnames(edger_fit$coefficients))) # for RNA library
dge_output <- topTags(lrt, n = Inf)$table %>%
  rownames_to_column('gene') %>%
  dplyr::rename("pvalue" = "PValue", "padj" = "FDR")
# (commenting out to reduce output files)
# filename <- file.path(output_dir, 'edger_analysis', 'stranded_vs_polya_dge_edger_result.csv')
# write_tsv(dge_output, file = filename)

# plot and save p-value histogram
p <- deseq2_pvals_histogram(res_df = dge_output,
                            xlab = 'stranded vs poly-A RNA-seq edgeR p-value',
                            ylab = 'Gene count',
                            title = paste0('Histogram of stranded vs poly-A RNA-seq\n',
                                           'differential gene expression edgeR p-values'))
filename <- file.path(output_dir, 'edger_analysis', 'stranded_vs_polya_dge_edger_histogram.pdf')
ggsave(filename = filename, plot = p, width = 8, height = 7, bg = "white")

# pull empirical control genes to be used in RUVg (not differential expressed)
emp_neg_ctrl_genes_dge_edgeR <- dge_output %>% 
  filter(padj > 0.05) %>%
  pull(gene)
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_dge_edgeR, prefix = "dge_empirical_genes", 
          diff_type = "edger", 
          output_dir = output_dir)

# 2.4 HK genes (HRT atlas) differentially expressed in DEG analysis
hk_genes_normals <- readRDS('input/hk_genes_normals.rds')
emp_neg_ctrl_hk_genes_dge_edgeR <- dge_output %>% 
  filter(padj < 0.05,
         gene %in% hk_genes_normals) %>%
  pull(gene)
print(paste0("Differentially expressed HK genes:", length(emp_neg_ctrl_hk_genes_dge_edgeR)))
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_hk_genes_dge_edgeR, prefix = "dge_empirical_hk_genes", 
          diff_type = "edger", 
          output_dir = output_dir)

# 2.5 RUVr using empirical control genes (dge from edgeR)
# Estimating the factors of unwanted variation using residuals
# uses residuals, e.g., from a first-pass GLM regression of the counts on the covariates of interest.
res <- residuals(edger_fit, type = "deviance") # edger_fit from line 169

# we can use all the genes to estimate the factors of unwanted variation.
all_genes <- rownames(seq_expr_set@assayData$counts)
ruvr_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = all_genes, 
          residuals = res, 
          output_dir = output_dir)

print("edgeR analysis finished.")
########################################### edgeR analysis ##############################################################
