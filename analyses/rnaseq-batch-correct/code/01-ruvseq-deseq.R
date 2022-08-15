# This script performs the following functions:
# 1. DESeq2 tumor-only analysis with RUVg by molecular subtype
# Authors: Komal Rathi, updated by Adam Kraya

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(DESeq2)
  library(RUVSeq)
  library(EDASeq)
  library(edgeR)
  library(stringr)
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

# extract parameters
opt <- parse_args(OptionParser(option_list = option_list))
dataset <- opt$dataset
cancer_group_value <- unlist(stringr::str_split(opt$cancer_group, ','))
cohort_value <- opt$cohort
k_value <- as.numeric(opt$k_value)

# establish directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, 'analyses', 'rnaseq-batch-correct')
plots_dir <- file.path(analysis_dir, 'plots', dataset)
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

output_dir <- file.path(analysis_dir, 'output', dataset)
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive=TRUE)
}

# source functions 
source('util/deseq2_pvals_histogram.R') # DESeq2 pval histograms
source('util/edaseq_plot.R') # PCA and UMAP clustering
source('util/ruvg_test.R') # function to run RUVg

# read histology
message("Reading in histologies file")
htl_df <- readr::read_tsv('../../data/v10/histologies.tsv')

# read expected counts
message("Reading RSEM expected counts file")
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
design <- model.matrix(~ molecular_subtype)
counts_object <- edgeR::DGEList(counts = cnt_df, group = molecular_subtype)
counts_object_filtered <- edgeR::filterByExpr(counts_object, group = molecular_subtype)
counts_object <- counts_object[counts_object_filtered, , keep.lib.sizes = FALSE]

# create new expression set
seq_expr_set <- EDASeq::newSeqExpressionSet(counts = round(counts_object$counts), phenoData = data.frame(bs_id, molecular_subtype, RNA_library, row.names = colnames(round(counts_object$counts))))

########################################### DESeq2 analysis ##############################################################

# 1.1. RUVg using hk genes in normals only (full HRT atlas)
message("Running RUVg-DESeq2 analysis for k = 1-5")
emp_neg_ctrl_genes_normals <- readRDS('input/hk_genes_normals.rds')
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_normals, prefix = "hk_genes_normals", 
          diff_type = "deseq2", 
          output_dir = output_dir,
          plot_dir = plots_dir,
          design_variable = "molecular_subtype",
          color_var = "RNA_library", shape_var = "RNA_library")

message("Batch correction DGE complete")
# 1.2 Non-batch corrected DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
message("Running non-batch corrected DESeq2 analysis")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts_object$counts), colData =  data.frame(bs_id, molecular_subtype), design = design)
dds <- DESeq2::DESeq(dds)
dge_output <- DESeq2::results(dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
dge_output <- dge_output %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>%
  arrange(padj) 

message("Writing DESeq2 non batch-corrected results")
readr::write_tsv(x = dge_output, file = file.path(output_dir, 'deseq2_analysis', paste0('DESeq2_noBatchC_', dataset, '_dge.tsv')))
readr::write_rds(x = dds, file = file.path(output_dir, 'deseq2_analysis', paste0('DESeq2_noBatchC_', dataset, '_dge.rds')))

# plot and save p-value histogram
message("Saving p-value histogram")
p <- deseq2_pvals_histogram(res_df = dge_output,
                            xlab = 'DGE RLE nbinomWaldTest p-value',
                            ylab = 'Gene count', title = paste0('Histogram of DESeq2 nbinomWaldTest p-values'))
filename <- file.path(plots_dir, 'dge_deseq2_histogram.pdf')
ggsave(filename = filename, plot = p, width = 8, height = 7, bg = "white")


print("DESeq2 analysis finished.")
########################################### DESeq2 analysis ##############################################################
