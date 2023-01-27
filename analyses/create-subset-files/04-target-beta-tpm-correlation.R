# Calculate representative probe-level correlations between RNA-Seq (TPM) and 
# methylation (Beta) for TARGET patients who have samples in both datasets 

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
                  tumor_descriptor == "Primary Tumor" & 
                  cohort == "TARGET" & 
                  (experimental_strategy == "RNA-Seq" | 
                     experimental_strategy == "Methylation"))

# Get methylation beta values
beta <- data.table::fread(file.path(results_dir, 
                                    "methylation-beta-values-matrix.tsv.gz"),
                          sep = "\t", 
                          showProgress = FALSE) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(Probe_ID, tidyselect::starts_with("TARGET-"))

# Get RNA-Seq expression TPM values for array probes in GENCODE version 38 
# (Ensembl 104) gene symbols
tpm <- readr::read_rds(file.path(data_dir, 
                                 "gene-expression-rsem-tpm-collapsed.rds")) %>% 
  dplyr::select(dplyr::starts_with("TARGET-")) %>% 
  tibble::rownames_to_column(var = "Gene_symbol")
tpm <- data.table::fread(
  file.path(results_dir, "methylation-probe-annotations.tsv.gz"), 
  select = c("Probe_ID", "Gene_symbol"), showProgress = FALSE)  %>% 
  tibble::as_tibble() %>% 
  dplyr::distinct()  %>%
  dplyr::inner_join(tpm, by = "Gene_symbol") %>% 
  dplyr::select(-Gene_symbol)

# Calculating representative probe-level correlations between RNA-Seq (TPM) and 
# preprocessed methylation (Beta) samples for the TARGET cohort
message("=====================================================================================")
message("Calculating probe-level correlations between beta and tpm values for TARGET samples")
message("=====================================================================================\n")
# Compute cancer type methylation to expression correlation for TARGET cohort
rnaseq_histologies <- histologies %>% 
  dplyr::filter(experimental_strategy == "RNA-Seq")
meth_histologies <- histologies %>% 
  dplyr::filter(experimental_strategy == "Methylation")
rm(histologies)
beta_tmp_correlation <- tibble::tibble()
for (cancer_type in unique(meth_histologies$cancer_group)) {
  message(c("Calculating probe-level correlations between beta and tpm values for TARGET, ", cancer_type, " samples...\n"))
  # Get cancer type patients with RNA-Seq data
  rnaseq_usi <-  rnaseq_histologies %>% 
    dplyr::filter(cancer_group == cancer_type) %>% 
    dplyr::pull(Kids_First_Participant_ID)
  rnaseq_usi <- unique(rnaseq_usi) 
  sample_ids <- meth_histologies %>% 
    dplyr::filter(cancer_group == cancer_type,
                  Kids_First_Participant_ID %in% rnaseq_usi) %>% 
    dplyr::pull(Kids_First_Biospecimen_ID)
  # Get cancer type beta values for patients with both RNA-Seq and 
  # methylation data
  cancer_type_beta <- beta %>% 
    dplyr::select(tidyselect::any_of(c("Probe_ID", sample_ids))) %>%
    tidyr::pivot_longer(-Probe_ID, names_to = "Sample", values_to = "Beta") %>% 
    tidyr::pivot_wider(names_from = Probe_ID, values_from = Beta) %>%
    dplyr::mutate(USI = stringr::str_split(Sample, "-", simplify = TRUE)[,3]) %>%
    dplyr::select(-Sample) %>%
    tidyr::pivot_longer(-USI, names_to = "Probe_ID", values_to = "Beta") %>%
    tidyr::pivot_wider(names_from = USI, values_from = Beta, values_fn = median)
  # Get cancer type tpm values for patients with both RNA-Seq and 
  # methylation data
  cancer_type_tpm <- tpm %>% 
    dplyr::select(tidyselect::contains(colnames(cancer_type_beta))) %>% 
    tidyr::pivot_longer(-Probe_ID, names_to = "Sample", values_to = "TPM") %>%
    tidyr::pivot_wider(names_from = Probe_ID, values_from = TPM, 
                       values_fn = median) %>%
    dplyr::mutate(USI = stringr::str_split(Sample, "-", simplify = TRUE)[,3]) %>%
    dplyr::select(-Sample) %>%
    tidyr::pivot_longer(-USI, names_to = "Probe_ID", values_to = "TPM") %>%
    tidyr::pivot_wider(names_from = USI, values_from = TPM, values_fn = median)
  # Keep array probes in cancer type beta values to that only present 
  # in cancer type tmp values
  cancer_type_beta <- cancer_type_beta %>% 
    dplyr::filter(Probe_ID %in% unique(cancer_type_tpm$Probe_ID))
  # calculate probe correlation between methylation beta values
  # RNA-Seq expression tpm values
  cancer_type_beta <- cancer_type_beta %>%
    tidyr::pivot_longer(-Probe_ID, names_to = "Sample", values_to = "Beta") %>%
    tidyr::pivot_wider(names_from = Probe_ID, values_from = Beta) %>%
    dplyr::arrange(Sample) %>%
    dplyr::select(-Sample)
  cancer_type_tpm <- cancer_type_tpm %>%
    tidyr::pivot_longer(-Probe_ID, names_to = "Sample", values_to = "TPM") %>%
    tidyr::pivot_wider(names_from = Probe_ID, values_from = TPM) %>%
    dplyr::arrange(Sample) %>%
    dplyr::select(-Sample)
  for (probe in colnames(cancer_type_tpm)) {
    if (sd(unlist(cancer_type_beta[,probe])) == 0 || 
        sd(unlist(cancer_type_tpm[,probe])) == 0) {
      next 
    }
    beta_tmp_correlation <-  beta_tmp_correlation %>%
      dplyr::bind_rows(
        tibble("Probe_ID" = probe,
               "RNA_Correlation" = cor(cancer_type_beta[,probe], 
                                       cancer_type_tpm[,probe][1]),
               "Disease" = cancer_type)
        )
  }
  rm(rnaseq_usi, sample_ids, cancer_type_beta, cancer_type_tpm)
}
rm(beta, tpm, rnaseq_histologies, meth_histologies) 

# Write probe-level correlations between methylation beta values
# RNA-Seq expression tpm values
message("Writing probe-level correlations to methylation-probe-beta-tpm-correlations.tsv file...\n")
beta_tmp_correlation %>% data.table::setDT() %>%
  data.table::fwrite(file.path(results_dir,
                               "methylation-probe-beta-tpm-correlations.tsv.gz"), 
                     sep="\t", compress = "auto")

message("Analysis Done..\n")