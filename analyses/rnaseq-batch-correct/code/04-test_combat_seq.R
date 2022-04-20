# function to perform ComBat_seq batch correction
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(sva)
  library(edgeR)
  library(EDASeq)
})

# parse parameters
option_list <- list(
  make_option(opt_str = "--dataset", type = "character",
              help = "Dataset for running differential gene expression analysis: match_pbta_hgg or match_target_all")
)
opt <- parse_args(OptionParser(option_list = option_list))
dataset <- opt$dataset
output_dir <- file.path('output', 'combatseq_output', dataset)
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source('util/edaseq_plot.R') # PCA and UMAP clustering
source('util/box_plots.R') # boxplots

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
colnames(cnt_df) <- gsub('TARGET-[0-9]{2}-', '', colnames(cnt_df))

# filter by expression and convert to DGElist
bs_id <- gsub('TARGET-[0-9]{2}-', '', selected_htl_df$Kids_First_Biospecimen_ID) # shorten bs ids for TARGET samples
patient_id <- factor(as.character(selected_htl_df$Kids_First_Participant_ID))
rna_library <- factor(as.character(selected_htl_df$RNA_library))
design <- model.matrix(~ 0 + patient_id + rna_library)
counts_object <- edgeR::DGEList(counts = cnt_df, 
                                samples = data.frame(patient_id, rna_library, bs_id), 
                                group = rna_library)
counts_object_filtered <- edgeR::filterByExpr(counts_object, group = rna_library)
counts_object <- counts_object[counts_object_filtered, , keep.lib.sizes = FALSE]

run_combat <- function(counts_object, genes_to_filter = NULL, prefix = ""){
  
  # extract counts from counts object
  if(is.null(genes_to_filter)){
    count_data <- counts_object$counts
  } else {
    count_data <- counts_object$counts[which(rownames(counts_object$counts) %in% genes_to_filter), ]
  }
  seq_expr_set <- EDASeq::newSeqExpressionSet(counts = round(count_data), phenoData = counts_object$samples)
  
  # before batch correction
  pca_p <- edaseq_plot(object = seq_expr_set, title = "Before ComBat_seq", type = "PCA")
  umap_p <- edaseq_plot(object = seq_expr_set, title = "Before ComBat_seq", type = "UMAP")
  boxplot_p <- box_plots(object = seq_expr_set, title = "Before ComBat_seq")
  
  # after batch correction
  counts_batch_corrected <- sva::ComBat_seq(counts = count_data, batch = rna_library, group = rep(1, ncol(count_data)), full_mod = FALSE)
  seq_expr_set_bc <- EDASeq::newSeqExpressionSet(counts = round(counts_batch_corrected), phenoData = counts_object$samples)
  pca_p_bc <- edaseq_plot(object = seq_expr_set_bc, title = "After ComBat_seq", type = "PCA")
  umap_p_bc <- edaseq_plot(object = seq_expr_set_bc, title = "After ComBat_seq", type = "UMAP")
  boxplot_p_bc <- box_plots(object = seq_expr_set_bc, title = "After ComBat_seq")
  
  # save both UMAP and PCA in a single file
  p <- ggpubr::ggarrange(pca_p, pca_p_bc, umap_p, umap_p_bc, common.legend = T, legend = "bottom")
  fname <- file.path(output_dir, paste0(prefix, '_clustering_with_and_without_bc.pdf'))
  ggsave(filename = fname, plot = p, width = 6, height = 6, device = "pdf", bg = "white")
  
  # save boxplots in a single file
  p <- ggpubr::ggarrange(boxplot_p, boxplot_p_bc, nrow = 2, common.legend = T, legend = "bottom")
  fname <- file.path(output_dir, paste0(prefix, '_boxplots_with_and_without_bc.pdf'))
  ggsave(filename = fname, plot = p, width = 12, height = 6, device = "pdf", bg = "white")
}

# all expressed genes
run_combat(counts_object = counts_object, genes_to_filter = NULL, prefix = "expressed_genes")

# subset to housekeeping genes from HRT protocol
hk_genes_normals <- readRDS('input/hk_genes_normals.rds')
run_combat(counts_object = counts_object, genes_to_filter = hk_genes_normals, prefix = "hk_genes_normals")

# subset to housekeeping genes identified in tumor samples per the HRT protocol
hk_genes_tumors_normals <- readRDS('input/hk_genes_tumor_normals.rds')
run_combat(counts_object = counts_object, genes_to_filter = hk_genes_tumors_normals, prefix = "hk_genes_tumors_normals")
