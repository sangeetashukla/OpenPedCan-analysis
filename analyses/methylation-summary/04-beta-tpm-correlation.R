# Calculate representative probe-level correlations between RNA-Seq (TPM) and 
# methylation (Beta) for patients who have samples in both datasets 

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

# Get required columns from histologies file for RNA-Seq and methyalation 
# sample IDs
required_cols <- c("Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", 
                   "experimental_strategy", "sample_type", "tumor_descriptor",
                   "cohort", "cancer_group")
histologies <- data.table::fread(file.path(data_dir, "histologies.tsv"), 
                                 sep = "\t", 
                                 select = required_cols, 
                                 showProgress = FALSE) %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(sample_type == "Tumor" & 
                  (tumor_descriptor == "Primary Tumor" | 
                     tumor_descriptor == "Initial CNS Tumor") & 
                  (experimental_strategy == "RNA-Seq" | 
                     experimental_strategy == "Methylation"))

# Get methylation beta values
beta <- readr::read_rds(file.path(results_dir, "methyl-beta-values.rds"))

# Get RNA-Seq expression TPM values for array probes in GENCODE version 38 
# (Ensembl 104) gene symbols
tpm <- readr::read_rds(file.path(data_dir, 
                                 "gene-expression-rsem-tpm-collapsed.rds")) %>% 
  tibble::rownames_to_column(var = "Gene_symbol")
tpm <- data.table::fread(
  file.path(results_dir, "methyl-probe-annotations.tsv.gz"), 
  select = c("Probe_ID", "Gene_symbol"), showProgress = FALSE)  %>% 
  tibble::as_tibble() %>% 
  dplyr::distinct()  %>%
  dplyr::inner_join(tpm, by = "Gene_symbol") %>% 
  dplyr::select(-Gene_symbol)

# Calculating representative probe-level correlations between RNA-Seq (TPM) and 
# preprocessed methylation (Beta) samples
message("========================================================================================")
message("Calculating probe-level correlations between methylation beta and expression tpm values")
message("========================================================================================\n")
# Compute cancer type methylation to expression correlation
rnaseq_histologies <- histologies %>% 
  dplyr::filter(experimental_strategy == "RNA-Seq")
meth_histologies <- histologies %>% 
  dplyr::filter(experimental_strategy == "Methylation")
rm(histologies)
beta_tmp_correlation <- tibble::tibble()
for (cancer_type in unique(meth_histologies$cancer_group)) {
  message(c("Calculating probe-level correlations between beta and tpm values, ", cancer_type, " samples...\n"))
  # Get cancer type patients with RNA-Seq data
  cancer_type_rnaseq_ids <-  rnaseq_histologies %>% 
    dplyr::filter(cancer_group == cancer_type) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, 
                  Kids_First_Participant_ID, cohort) %>% 
    dplyr::rename(RNASeq_ID = Kids_First_Biospecimen_ID)
  cancer_type_ids <- meth_histologies %>% 
    dplyr::filter(cancer_group == cancer_type) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, 
                  Kids_First_Participant_ID,  cohort) %>% 
    dplyr::rename(Meth_ID = Kids_First_Biospecimen_ID) %>% 
    dplyr::inner_join(cancer_type_rnaseq_ids, 
                      by = c("Kids_First_Participant_ID", "cohort")) %>% 
    dplyr::rename(Patient_ID = Kids_First_Participant_ID) %>%
    dplyr::select(-cohort)
 
  # Get cancer type tpm values for patients with both RNA-Seq and 
  # methylation data and only keep array probes present in cancer type
  #  that have standard deviation == 0
  cancer_type_tpm <- tpm %>% 
    dplyr::select(tidyselect::any_of(c("Probe_ID",
                                       unique(cancer_type_ids$RNASeq_ID)))) %>%
    tidyr::pivot_longer(-Probe_ID, names_to = "RNASeq_ID", values_to = "TPM") %>% 
    dplyr::left_join(cancer_type_ids, by = "RNASeq_ID") %>% 
    dplyr::select(-c("Meth_ID", "RNASeq_ID")) %>%
    tidyr::pivot_wider(names_from = Patient_ID, 
                       values_from = TPM, values_fn = median) %>% 
    tidyr::pivot_longer(-Probe_ID, 
                        names_to = "Patient_ID", values_to = "TPM") %>%
    tidyr::pivot_wider(names_from = Probe_ID, values_from = TPM) %>% 
    purrr::discard(~ all(is.numeric(.x) & sd(.x) == 0))

  # Get cancer type beta values for patients with both RNA-Seq and 
  # methylation data and only keep array probes present in cancer type
  # tmp values that have standard deviation > 0
  cancer_type_beta <- beta %>% 
    dplyr::select(tidyselect::any_of(c("Probe_ID", 
                                       unique(cancer_type_ids$Meth_ID)))) %>%
    dplyr::filter(!is.na(.)) %>% 
    dplyr::filter(Probe_ID %in% unique(colnames(cancer_type_tpm))) %>% 
    tidyr::pivot_longer(-Probe_ID, names_to = "Meth_ID", values_to = "Beta") %>%
    dplyr::left_join(cancer_type_ids, by = "Meth_ID") %>% 
    dplyr::select(-c("Meth_ID", "RNASeq_ID")) %>%
    tidyr::pivot_wider(names_from = Patient_ID, 
                       values_from = Beta, values_fn = median) %>% 
    tidyr::pivot_longer(-Probe_ID, 
                        names_to = "Patient_ID", values_to = "Beta") %>%
    tidyr::pivot_wider(names_from = Probe_ID, values_from = Beta) %>% 
    purrr::discard(~ all(is.numeric(.x) & sd(.x) == 0))

  # Only keep cancer type tpm values for patient IDs and array probes that are 
  # in beta values then sort patient IDs and array probes to be in the same
  # order in both matrices. 
  # Once again drop any array probes with standard deviation == 0 as result of
  # filtering 
  cancer_type_tpm <- cancer_type_tpm %>% 
    dplyr::filter(Patient_ID %in% cancer_type_beta$Patient_ID) %>%
    purrr::discard(~ all(is.numeric(.x) & sd(.x) == 0)) %>% 
    dplyr::arrange(Patient_ID)
  cancer_type_beta <- cancer_type_beta %>% 
    dplyr::select(tidyselect::any_of(colnames(cancer_type_tpm))) %>%
    dplyr::arrange(Patient_ID) %>%
    tibble::column_to_rownames(var = "Patient_ID") %>% 
    dplyr::select(sort(tidyselect::peek_vars()))
  cancer_type_tpm <- cancer_type_tpm %>% 
    tibble::column_to_rownames(var = "Patient_ID") %>%
    dplyr::select(tidyselect::all_of(colnames(cancer_type_beta))) %>%
    dplyr::select(sort(tidyselect::peek_vars()))
  
  # calculate probe correlation between methylation beta values
  # RNA-Seq expression tpm values
  beta_tmp_correlation <-  beta_tmp_correlation %>%
    dplyr::bind_rows(
      tibble::tibble(
        "Probe_ID" = names(cancer_type_beta),
        "RNA_Correlation" = 
          sapply(1:ncol(cancer_type_beta),
                 function(i) cor(cancer_type_beta[,i], cancer_type_tpm[,i])),
      ) %>% 
        dplyr::mutate(Disease = cancer_type)
    )
  rm(cancer_type_rnaseq_ids, cancer_type_ids, cancer_type_beta, cancer_type_tpm)
}
rm(beta, tpm, rnaseq_histologies, meth_histologies) 

# Write probe-level correlations between methylation beta values
# RNA-Seq expression tpm values
message("Writing probe-level correlations to methyl-probe-beta-tpm-correlations.tsv file...\n")
beta_tmp_correlation %>% data.table::setDT() %>%
  data.table::fwrite(file.path(results_dir,
                               "methyl-probe-beta-tpm-correlations.tsv.gz"), 
                     sep="\t", compress = "auto")

message("Analysis Done..\n")