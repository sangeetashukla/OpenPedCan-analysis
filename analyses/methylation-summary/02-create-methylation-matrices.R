# Create combined methylation `Beta-values` and `M-values` matrices for all
# histologies (cancer types)

# Eric Wafula for Pediatric OpenTargets
# 03/18/2022

# Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to module and results directories
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "analyses", "methylation-summary")
results_dir <- file.path(module_dir, "results")

# Get primary tumor methylation sample IDs
required_cols <- c("Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", 
                   "experimental_strategy", "sample_type", "tumor_descriptor",
                   "cohort", "cancer_group")
primary_tumors <- data.table::fread(file.path(data_dir, "histologies.tsv"), 
                                 sep = "\t", 
                                 select = required_cols, 
                                 showProgress = FALSE) %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(sample_type == "Tumor" & tumor_descriptor == "Primary Tumor" &  
                  experimental_strategy == "Methylation") %>% 
  dplyr::pull(Kids_First_Biospecimen_ID)


# Get array annotation probes corresponding to GENCODE version 38 (Ensembl 104) 
# gene symbols
annot_file <- file.path(results_dir, 
                        "methylation-array-probes-annotations.tsv.gz")
probe_ids <- data.table::fread(annot_file, select = "Probe_ID", 
                               showProgress = FALSE) %>% 
  dplyr::distinct() %>% 
  dplyr::pull(Probe_ID)

# Function to create methylation matrix
create_matrix <- function(meth_file, probe_ids) {
  required_cols <- c("Probe", "Sample_Name", "Meth_Value")
  meth_table <- data.table::fread(meth_file, 
                                  select = required_cols, 
                                  showProgress = FALSE) %>% 
    tibble::tibble() %>% 
    dplyr::rename(Probe_ID = Probe) %>% 
    tidyr::pivot_wider(names_from = Sample_Name, values_from = Meth_Value) %>% 
    dplyr::filter(Probe_ID %in% probe_ids)
  return(meth_table)
}

# Creating Beta-values matrix for all  preprocessed methylation samples
message("===================================================================")
message(c("Creating merged Beta-values matrix for all methylation samples"))
message("===================================================================\n")

# NBL beta-values samples
message("Creating Beta-values matrix for Neuroblastoma (NBL) samples...\n")
nbl_matrix <- create_matrix(
  file.path(data_dir, "NBL-beta-values-methylation.tsv.gz"), probe_ids)

# OS beta-values samples
message("Creating Beta-values matrix for Osteosarcoma (OS) samples...\n")
os_matrix <- create_matrix(
  file.path(data_dir, "OS-beta-values-methylation.tsv.gz"), probe_ids)
beta_matrix <- nbl_matrix %>% dplyr::full_join(os_matrix, by = "Probe_ID")
rm(nbl_matrix, os_matrix)

# CCSK beta-values samples
message("Creating Beta-values matrix for Clear Cell Sarcoma of the Kidney (CCSK) samples...\n")
ccsk_matrix <- create_matrix(
  file.path(data_dir, "CCSK-beta-values-methylation.tsv.gz"), probe_ids)
beta_matrix <- beta_matrix %>% dplyr::full_join(ccsk_matrix, by = "Probe_ID")
rm(ccsk_matrix)

# WT beta-values samples 
message("Creating Beta-values matrix for Wilms Tumor (WT) samples...\n")
wt_matrix <- create_matrix(
  file.path(data_dir, "WT-beta-values-methylation.tsv.gz"), probe_ids)
beta_matrix <- beta_matrix %>% dplyr::full_join(wt_matrix, by = "Probe_ID")
rm(wt_matrix)

# AML beta-values samples 
message("Creating Beta-values matrix for Acute Myeloid Leukemia (AML) samples...\n")
aml_matrix <- create_matrix(
  file.path(data_dir, "AML450k-beta-values-methylation.tsv.gz"), probe_ids)
beta_matrix <- beta_matrix %>% dplyr::full_join(aml_matrix, by = "Probe_ID")
rm(aml_matrix)

# write merged beta-values to file
message("Writing merged Beta-values matrix to methylation-beta-values-matrix.tsv file...\n")
beta_matrix %>% dplyr::select(tidyselect::any_of(c("Probe_ID", primary_tumors))) %>% 
  data.table::setDT() %>%
  data.table::fwrite(file.path(results_dir,
                               "methylation-beta-values-matrix.tsv.gz"), 
                     sep="\t", compress = "auto")
rm(beta_matrix)

# Creating M-values matrix for all methylation samples
message("=================================================================")
message("Creating merged M-values matrix for all methylation samples")
message("=================================================================\n")

# NBL m-values samples
message("Creating M-values matrix for Neuroblastoma (NBL) samples...\n")
nbl_matrix <- create_matrix(
  file.path(data_dir, "NBL-m-values-methylation.tsv.gz"), probe_ids)

# OS m-values samples
message("Creating M-values matrix for Osteosarcoma (OS) samples...\n")
os_matrix <- create_matrix(
  file.path(data_dir, "OS-m-values-methylation.tsv.gz"), probe_ids)
m_matrix <- nbl_matrix %>% dplyr::full_join(os_matrix, by = "Probe_ID")
rm(nbl_matrix, os_matrix)

# CCSK m-values samples
message("Creating M-values matrix for Clear Cell Sarcoma of the Kidney (CCSK) samples...\n")
ccsk_matrix <- create_matrix(
  file.path(data_dir, "CCSK-m-values-methylation.tsv.gz"), probe_ids)
m_matrix <- m_matrix %>% dplyr::full_join(ccsk_matrix, by = "Probe_ID")
rm(ccsk_matrix)

# WT m-values samples 
message("Creating M-values matrix for Wilms Tumor (WT) samples...\n")
wt_matrix <- create_matrix(
  file.path(data_dir, "WT-m-values-methylation.tsv.gz"), probe_ids)
m_matrix <- m_matrix %>% dplyr::full_join(wt_matrix, by = "Probe_ID")
rm(wt_matrix)

# AML m-values samples 
message("Creating M-values matrix for Acute Myeloid Leukemia (AML) samples...\n")
aml_matrix <- create_matrix(
  file.path(data_dir, "AML450k-m-values-methylation.tsv.gz"), probe_ids)
m_matrix <- m_matrix %>% dplyr::full_join(aml_matrix, by = "Probe_ID")
rm(aml_matrix)

# write merged m-values to file
message("Writing merged M-values matrix to methylation-m-values-matrix.tsv file...\n")
m_matrix %>% dplyr::select(tidyselect::any_of(c("Probe_ID", primary_tumors))) %>% 
  data.table::setDT() %>%
  data.table::fwrite(file.path(results_dir,
                               "methylation-m-values-matrix.tsv.gz"), 
                     sep="\t", compress = "auto")
rm(m_matrix)

message("Analysis Done..\n")