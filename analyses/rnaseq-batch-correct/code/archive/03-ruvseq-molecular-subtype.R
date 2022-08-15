# This script performs the following functions:
# 1. DESeq2 analysis with RUVg 
# 2. edgeR analysis with RUVg and RUVr

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(msigdbr)
  library(data.table)
})

# parse parameters
option_list <- list(
  make_option(opt_str = "--dataset", type = "character",
              help = "Give any name to the dataset to create output folder"),
  make_option(opt_str = "--cancer_group_value", type = "character",
              help = "cancer group value"),
  make_option(opt_str = "--cohort_value", type = "character",
              help = "cohort value"),
  make_option(opt_str = "--k_value", type = "character",
              help = "number of K values to evaluate")
)
opt <- parse_args(OptionParser(option_list = option_list))
dataset <- opt$dataset
cancer_group_value <- opt$cancer_group
cohort_value <- opt$cohort
k_value <- as.numeric(opt$k_value)

output_dir <- file.path('output', dataset)
dir.create(output_dir, showWarnings = F, recursive = T)

input_dir <- file.path('input')

dataset = 'target_nbl'
cancer_group_value <- 'Neuroblastoma'
cohort_value <- 'TARGET'
k_value <- 5

# source functions 
source('util/deseq2_pvals_histogram.R') # DESeq2 pval histograms
source('util/edaseq_plot.R') # PCA and UMAP clustering
source('util/box_plots.R') # boxplots for samples
source('util/ruvg_test.R') # function to run RUVg
source('util/ruvr_test.R') # function to run RUVr (only works with edgeR)

# read histology
htl_df <- readr::read_tsv('../../data/v10/histologies.tsv')

# read expected counts
cnt_df <- readRDS('../../data/v10/gene-counts-rsem-expected_count-collapsed.rds')

# filter histology
selected_htl_df <- htl_df %>%
  filter(experimental_strategy == "RNA-Seq",
         sample_type == "Tumor",
         cancer_group %in% cancer_group_value,
         !is.na(molecular_subtype)) 

# filter expression
cnt_df <- cnt_df %>%
  dplyr::select(selected_htl_df$Kids_First_Biospecimen_ID)
cnt_df <- cnt_df[rowSums(cnt_df) > 0,]

# filter by expression and convert to DGElist
bs_id <- gsub('TARGET-[0-9]{2}-', '', selected_htl_df$Kids_First_Biospecimen_ID) # shorten bs ids for TARGET samples
molecular_subtype <- factor(as.character(selected_htl_df$molecular_subtype))
RNA_library <- factor(as.character(selected_htl_df$RNA_library))
design <- model.matrix(~molecular_subtype)
counts_object <- edgeR::DGEList(counts = cnt_df, group = molecular_subtype)
counts_object_filtered <- edgeR::filterByExpr(counts_object, group = molecular_subtype)
counts_object <- counts_object[counts_object_filtered, , keep.lib.sizes = FALSE]

# create new expression set
seq_expr_set <- EDASeq::newSeqExpressionSet(counts = round(counts_object$counts), phenoData = data.frame(bs_id, molecular_subtype, RNA_library, row.names = colnames(round(counts_object$counts))))
pca_p <- edaseq_plot(object = seq_expr_set, title = "Before normalization", type = "PCA", color_var = "RNA_library", shape_var = "RNA_library")
umap_p <- edaseq_plot(object = seq_expr_set, title = "Before normalization", type = "UMAP", color_var = "RNA_library", shape_var = "RNA_library")
boxplot_p <- box_plots(object = seq_expr_set, title = "Before normalization", facet_var = "RNA_library", color_var = "RNA_library")

# from https://support.bioconductor.org/p/95805/#95808
# The "betweenLaneNormalization" is just a poorly named function that perform between-sample normalization, 
# independently of the fact that they are organized in lanes, plates, etc. 
# It is used in RUVSeq just to make sure that the factors of unwanted variation don't capture the difference in library size that should be captured by the size factors. 
# normalize the data using upper-quartile (UQ) normalization
seq_expr_set <- EDASeq::betweenLaneNormalization(seq_expr_set, which = "upper")
pca_p_norm <- edaseq_plot(object = seq_expr_set, title = "After UQ normalization", type = "PCA", color_var = "RNA_library", shape_var = "RNA_library")
umap_p_norm <- edaseq_plot(object = seq_expr_set, title = "After UQ normalization", type = "UMAP", color_var = "RNA_library", shape_var = "RNA_library")
boxplot_p_norm <- box_plots(object = seq_expr_set, title = "After UQ normalization", facet_var = "RNA_library", color_var = "RNA_library")

# save both UMAP and PCA in a single file
p <- ggpubr::ggarrange(pca_p, pca_p_norm, umap_p, umap_p_norm, common.legend = T, legend = "bottom")
fname <- file.path(output_dir, 'clustering_with_and_without_norm_rna_lib.pdf')
ggsave(filename = fname, plot = p, width = 6, height = 6, device = "pdf", bg = "white")

# save boxplots in a single file
p <- ggpubr::ggarrange(boxplot_p, boxplot_p_norm, nrow = 2, common.legend = T, legend = "bottom")
fname <- file.path(output_dir, 'boxplots_with_and_without_norm_rna_lib.pdf')
ggsave(filename = fname, plot = p, width = 12, height = 6, device = "pdf", bg = "white")

########################################### DESeq2 analysis ##############################################################

# 1.1. RUVg using hk genes in tumor + normals
emp_neg_ctrl_genes_hk_tumor_normal <- readRDS('input/hk_genes_tumor_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_hk_tumor_normal, prefix = "hk_genes_tumor_normals", 
          diff_type = "deseq2", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "RNA_library", shape_var = "RNA_library")

# 1.2. RUVg using hk genes in normals only (full HRT atlas)
emp_neg_ctrl_genes_normals <- readRDS('input/hk_genes_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_normals, prefix = "hk_genes_normals", 
          diff_type = "deseq2", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "RNA_library", shape_var = "RNA_library")

# Retrieve gene symbols of neuroblastoma and hgg use cases
if(cancer_group_value == 'Neuroblastoma'){
  
  msig.gs <- msigdbr::msigdbr(species = "Homo sapiens", category = 'C2')
  msig.gs <- msig.gs %>%
    dplyr::filter(gs_name %in% grep('KIM_MYCN_AMPLIFICATION_TARGETS', gs_name, value = TRUE)) %>%
    data.table::setDT(.)
  
  readr::write_rds(msig.gs$human_gene_symbol, file.path(input_dir, 'MYCN_targets_M2919_M18532.rds'))
  
} else if(cancer_group_value %in% c('	
High-grade glioma/astrocytoma', 'Diffuse midline glioma')){
  
  genes.hgg <- data.table::fread(file.path(input_dir, '12915_2022_1324_MOESM4_ESM.txt'))
  genes.hgg <- genes.hgg[,c('External Gene Name','K27M  vs WT day5  adj p-value'), with = FALSE] %>%
    .[`K27M  vs WT day5  adj p-value` < 0.05]
  genes.hgg <- genes.hgg$`External Gene Name`
}

# 1.3. RUVg using empirical control genes (dge from DESeq2 analysis)
# From DESeq2 documentation:
# 1. DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts(seq_expr_set)), colData =  data.frame(bs_id, molecular_subtype), design = design)
dds <- DESeq2::DESeq(dds)
dge_output <- DESeq2::results(dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
dge_output <- dge_output %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>%
  arrange(padj) 

readr::write_tsv(x = dge_output, file = file.path(output_dir, paste0('DESeq2_noBatch_', dataset, '_dge.tsv')))

# plot and save p-value histogram
p <- deseq2_pvals_histogram(res_df = dge_output,
                            xlab = 'DGE RLE nbinomWaldTest p-value',
                            ylab = 'Gene count', title = paste0('Histogram of DESeq2 nbinomWaldTest p-values'))
filename <- file.path(output_dir, 'deseq2_analysis', 'dge_deseq2_histogram.pdf')
ggsave(filename = filename, plot = p, width = 8, height = 7, bg = "white")


# pull empirical control genes to be used in RUVg (not differential expressed)
emp_neg_ctrl_genes_dge_deseq2 <- dge_output %>% 
  filter(padj > 0.05) %>%
  pull(gene)
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value,
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_dge_deseq2,  prefix = "dge_empirical_genes", 
          diff_type = "deseq2", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

# 1.4 HK genes (HRT atlas) differentially expressed in DEG analysis
hk_genes_normals <- readRDS('input/hk_genes_normals.rds')
emp_neg_ctrl_hk_genes_dge_deseq2 <- dge_output %>% 
  filter(padj < 0.05,
         gene %in% hk_genes_normals) %>%
  pull(gene)
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value,
          emp_neg_ctrl_genes = emp_neg_ctrl_hk_genes_dge_deseq2,  prefix = "dge_empirical_hk_genes", 
          diff_type = "deseq2", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

print("DESeq2 analysis finished.")
########################################### DESeq2 analysis ##############################################################

########################################### edgeR analysis ##############################################################

# 2.1. RUVg using hk genes in tumor + normals
emp_neg_ctrl_genes_hk_tumor_normal <- readRDS('input/hk_genes_tumor_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_hk_tumor_normal, prefix = "hk_genes_tumor_normals", 
          diff_type = "edger", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

# 2.2. RUVg using hk genes in normals only (full HRT atlas)
emp_neg_ctrl_genes_normals <- readRDS('input/hk_genes_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_normals, prefix = "hk_genes_normals", 
          diff_type = "edger", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

# 2.3 RUVg using empirical control genes (dge from edgeR)
y <- edgeR::DGEList(counts = counts(seq_expr_set))
y <- edgeR::calcNormFactors(y, method = "upperquartile")
y <- edgeR::estimateGLMCommonDisp(y, design)
y <- edgeR::estimateGLMTagwiseDisp(y, design)
edger_fit <- edgeR::glmFit(y, design)

# export results
lrt <- edgeR::glmLRT(edger_fit, coef = grep('molecular_subtype', colnames(edger_fit$coefficients))) # for RNA library
dge_output <- edgeR::topTags(lrt, n = Inf)$table %>%
  rownames_to_column('gene') %>%
  dplyr::rename("pvalue" = "PValue", "padj" = "FDR")

# plot and save p-value histogram
p <- deseq2_pvals_histogram(res_df = dge_output,
                            xlab = 'edgeR p-value',
                            ylab = 'Gene count',
                            title = paste0('Histogram of edgeR p-values'))
filename <- file.path(output_dir, 'edger_analysis', 'dge_edger_histogram.pdf')
ggsave(filename = filename, plot = p, width = 8, height = 7, bg = "white")

# pull empirical control genes to be used in RUVg (not differential expressed)
emp_neg_ctrl_genes_dge_edgeR <- dge_output %>% 
  filter(padj > 0.05) %>%
  pull(gene)
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_dge_edgeR, prefix = "dge_empirical_genes", 
          diff_type = "edger", 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

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
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

# 2.5 RUVr using empirical control genes (dge from edgeR)
# Estimating the factors of unwanted variation using residuals
# uses residuals, e.g., from a first-pass GLM regression of the counts on the covariates of interest.
res <- residuals(edger_fit, type = "deviance") # edger_fit from line 169

# we can use all the genes to estimate the factors of unwanted variation.
all_genes <- rownames(seq_expr_set@assayData$counts)
ruvr_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = all_genes, 
          residuals = res, 
          output_dir = output_dir,
          design_variable = "molecular_subtype",
          color_var = "molecular_subtype", shape_var = "molecular_subtype")

print("edgeR analysis finished.")
########################################### edgeR analysis ##############################################################
