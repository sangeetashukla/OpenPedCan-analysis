# Author: Komal S. Rathi
# Function: Batch correction using sva::ComBat and sva::ComBat_seq

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sva))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "rnaseq-batch-correct")
output_dir <- file.path(analyses_dir, "output")

# parameters
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "comma-separated list of expression matrices to combine (TPM, FPKM or expected counts) (.rds)"),
  make_option(c("--type"), type = "character",
              help = "Type of expression data: TPM, FPKM or expected_count"),
  make_option(c("--metadata"), type = "character",
              help = "comma-separated list of sample metadata/histologies files to combine (.tsv)"),
  make_option(c("--id_col"), type = "character",
              help = "identifier column from metadata file that matches with expression matrix columns"),
  make_option(c("--batch_cols"), type = "character",
              help = "comma-separated list of columns from meta file to be used to create batch variable"),
  make_option(c("--other_cols_to_keep"), type = "character",
              help = "comma-separated list of columns to keep in the combined metadata in addition to the sample identifier and batch variables"),
  make_option(c("--output_prefix"), type = "character",
              help = "prefix for output files"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
mat <- opt$mat
mat <- trimws(strsplit(mat,",")[[1]]) 
metadata <- opt$metadata
metadata <- trimws(strsplit(metadata,",")[[1]]) 
id_col <- opt$id_col
batch_cols <- opt$batch_cols
batch_cols <- trimws(strsplit(batch_cols,",")[[1]]) 
other_cols_to_keep <- opt$other_cols_to_keep
other_cols_to_keep <- trimws(strsplit(other_cols_to_keep,",")[[1]]) 
type <- opt$type
output_prefix <- opt$output_prefix

# output files
uncorrected_outfile <- file.path(output_dir, paste0(output_prefix, '_uncorrected.rds'))
corrected_outfile <- file.path(output_dir, paste0(output_prefix, '_corrected.rds'))
metadata_with_batch_info <- file.path(output_dir, paste0(output_prefix, '_metadata.tsv'))

# source functions
source(file.path(analyses_dir, "util", "combine_mat.R"))
source(file.path(analyses_dir, "util", "combine_meta.R"))

# combine matrices
uncorrected_mat <- combine_mat(mat)

# combine metadata file 
# the output combined file will have the sample identifier and batch columns along with any other specified columns 
combined_metadata <- combine_meta(metadata, cols = c(id_col, batch_cols, other_cols_to_keep))

# combine the batch columns into a batch variable 
combined_metadata <- combined_metadata %>%
  unite('batch', batch_cols, remove = FALSE) %>%
  mutate(tmp = get(id_col)) %>%
  column_to_rownames('tmp')

# take an intersection if metadata and expression have different number of samples
samples_to_use <- intersect(rownames(combined_metadata), colnames(uncorrected_mat))
print(paste0("Samples to use:", length(samples_to_use)))

# subset matrix columns and metadata rows 
combined_metadata <- combined_metadata %>%
  filter(get(id_col) %in% samples_to_use)
uncorrected_mat <- uncorrected_mat %>%
  dplyr::select(all_of(samples_to_use))

# double-check if everything is lined-up between expression matrix and metadata
if(identical(rownames(combined_metadata), colnames(uncorrected_mat))){
  print("Matching dimensions")
  print("Save uncorrected file...")
  saveRDS(uncorrected_mat, file = uncorrected_outfile)
  print("Save metadata with batch info...")
  write.table(combined_metadata, file = metadata_with_batch_info, sep = "\t", quote = F)
} else {
  print("Check inputs")
  break
}

# batch correct using ComBat
print("Batch-correct uncorrected file..")
if(type != "expected_count"){
  # for fpkm/tpm we will use log2(x+1)
  corrected_mat <- ComBat(dat = log2(uncorrected_mat + 1), batch = combined_metadata$batch)
  # back-transform for down-stream analyses
  corrected_mat <- 2^(corrected_mat) 
  corrected_mat <- as.data.frame(corrected_mat)
} else {
  # for expected counts, ComBat_seq must be used
  corrected_mat <- ComBat_seq(counts = as.matrix(combined_mat), batch = combined_metadata$batch)
}

print("Save corrected file...")
saveRDS(corrected_mat, file = corrected_outfile)
