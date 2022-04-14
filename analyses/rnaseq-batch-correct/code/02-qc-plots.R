# Author: Komal S. Rathi
# Function: QC plots for batch correction

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(ggplot2))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "rnaseq-batch-correct")
plots_dir <- file.path(analyses_dir, "plots")

# source functions
source(file.path(analyses_dir, "util", "pubTheme.R"))
source(file.path(analyses_dir, "util", "clustering_plots.R"))
source(file.path(analyses_dir, "util", "density_plots.R"))

# parameters
option_list <- list(
  make_option(c("--uncorrected_mat"), type = "character",
              help = "combined expression matrix with multiple batches (.rds)"),
  make_option(c("--corrected_mat"), type = "character",
              help = "corrected expression matrix with multiple batches (.rds)"),
  make_option(c("--combined_metadata"), type = "character",
              help = "combined metadata file with batch information (.tsv)"),
  make_option(c("--var_prop"), type = "numeric",
              help = "proportion of most variable genes to be used"),
  make_option(c("--sample_id"), type = "character",
              help = "sample identifier column in metadata file matching column names in expression datasets"),
  make_option(c("--hk_genes"), type = "character",
              help = "path to housekeeping genes from Housekeeping Transcript Atlas (.csv)"),
  make_option(c("--clustering_type"), type = "character",
              help = "type of clustering to use: umap or tsne"),
  make_option(c("--plots_prefix"), type = "character",
              help = "prefix for clustering plots"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
uncorrected_mat <- opt$uncorrected_mat
corrected_mat <- opt$corrected_mat
combined_metadata <- opt$combined_metadata
var_prop <- opt$var_prop
var_prop <- as.numeric(var_prop)
sample_id <- opt$sample_id
hk_genes <- opt$hk_genes
clustering_type <- opt$clustering_type
plots_prefix <- opt$plots_prefix

# output files
clustering_plots <- file.path(plots_dir, paste0(plots_prefix, '_', clustering_type, '.pdf'))
density_plots <- file.path(plots_dir, paste0(plots_prefix, '_density.pdf'))

# read input data
uncorrected_mat <- readRDS(uncorrected_mat)
corrected_mat <- readRDS(corrected_mat)
combined_metadata <- read.delim(combined_metadata, stringsAsFactors = F, check.names = F)

# housekeeping genes (from Housekeeping Transcript Atlas)
hk_genes <- read.csv(hk_genes, sep = ";")
hk_genes <- unique(hk_genes$Gene.name)

# check if everything is lined-up between expression matrices and metadata
if(identical(rownames(combined_metadata), colnames(uncorrected_mat)) & 
   identical(rownames(combined_metadata), colnames(corrected_mat))){
  print("Matching dimensions")
} else {
  print("Check inputs")
  break
}

# t-SNE (uncorrected matrix)
# most variable genes
n <- round(var_prop/100*nrow(uncorrected_mat))
vars <- apply(uncorrected_mat, 1, var)
varorder <- order(vars, decreasing = T)
uncorrected_mat_mostvar <- uncorrected_mat[varorder,] %>% head(n)
p <- clustering_plot(mat = uncorrected_mat_mostvar, 
                     metadata = combined_metadata, 
                     color_var = 'group', 
                     shape_var = 'batch',
                     type = clustering_type,
                     title = paste0('Clustering before batch correction\n', var_prop, '% most variable genes'))

# house keeping genes
uncorrected_mat_hk_genes <- uncorrected_mat %>%
  rownames_to_column('gene') %>%
  filter(gene %in% hk_genes) %>%
  column_to_rownames('gene')
r <- clustering_plot(mat = uncorrected_mat_hk_genes, 
                     metadata = combined_metadata, 
                     color_var = 'group', 
                     shape_var = 'batch',
                     type = clustering_type,
                     title = 'Clustering before batch correction\nHousekeeping genes')

# t-SNE (corrected matrix)
# most variable genes
n <- round(var_prop/100*nrow(corrected_mat))
vars <- apply(uncorrected_mat, 1, var)
varorder <- order(vars, decreasing = T)
corrected_mat_mostvar <- uncorrected_mat[varorder,] %>% head(n)
q <- clustering_plot(mat = corrected_mat_mostvar, 
                     metadata = combined_metadata, 
                     color_var = 'group', 
                     shape_var = 'batch',
                     type = clustering_type,
                     title = paste0('Clustering after batch correction\n', var_prop, '% most variable genes'))

# house keeping genes
corrected_mat_hk_genes <- corrected_mat %>%
  rownames_to_column('gene') %>%
  filter(gene %in% hk_genes) %>%
  column_to_rownames('gene')
s <- clustering_plot(mat = corrected_mat_hk_genes, 
                     metadata = combined_metadata, 
                     color_var = 'group', 
                     shape_var = 'batch',
                     type = clustering_type,
                     title = 'Clustering after batch correction\nHousekeeping genes')

# save t-SNE plots
ggarrange(p, q, r, s, common.legend = T, legend = "right") %>%
  ggexport(filename = clustering_plots, width = 14, height = 10)

# distribution of housekeeping genes
# uncorrected  mat
uncorrected_mat_hk_genes <- uncorrected_mat_hk_genes %>% 
  rownames_to_column('gene') %>% 
  gather(-gene, key = !!sample_id, value = 'value') %>%
  inner_join(combined_metadata, by = sample_id)
p <- density_plot(mat = uncorrected_mat_hk_genes, 
                  color_var = 'batch', 
                  title = 'House Keeping Genes (Before ComBat correction)',
                  xlab = 'log2(Uncorrected value + 1)')

# corrected mat
corrected_mat_hk_genes <- corrected_mat_hk_genes %>% 
  as.data.frame() %>%
  rownames_to_column('gene') %>% 
  gather(-gene, key = !!sample_id, value = 'value') %>%
  inner_join(combined_metadata, by = sample_id)
q <- density_plot(mat = corrected_mat_hk_genes, 
                  color_var = 'batch', 
                  title = 'House Keeping Genes (After ComBat correction)',
                  xlab = 'log2(ComBat corrected value + 1)')

# save plots
ggarrange(p, q, ncol = 2, common.legend = T, legend = "right") %>%
  ggexport(filename = density_plots, width = 14, height = 6)
