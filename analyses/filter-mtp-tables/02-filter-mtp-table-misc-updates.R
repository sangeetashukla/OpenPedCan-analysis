
# Author: Sangeeta Shukla for Pediatric OpenTargets
# Purpose: Remove Ensembl (ESNG) gene identifier in the mutation frequency tables, including SNV, CNV and fusion that are not in GENCODE v39 and Ensembl package 104.


## Usage

## To run this from the command line, use:
## Rscript --vanilla filter-mpt-table-misc-updates.R
## This assumes you are in the modules directory of the repository, OpenPedCan-analysis/analyses/filter-mutation-frequency-tables._


# Load R analysis packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))

# Magrittr pipe
`%>%` <- dplyr::`%>%`




# Set up directories for input and output files
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses")
module_dir <- file.path(analyses_dir, "filter-mtp-tables")
input_dir <- file.path(module_dir, "results") ##Input data to be updated is in the results directory


# Read files to be updated
fusion_freq <- read_tsv(file.path(input_dir, "putative-oncogene-fusion-freq.tsv.gz"))
fused_gene_freq <- read_tsv(file.path(input_dir,"putative-oncogene-fused-gene-freq.tsv.gz"))
var_snv_freq <- read_tsv(file.path(input_dir,"variant-level-snv-consensus-annotated-mut-freq.tsv.gz"))
gene_snv_freq <- read_tsv(file.path(input_dir,"gene-level-snv-consensus-annotated-mut-freq.tsv.gz"))
gene_cnv_freq <- read_tsv(file.path(input_dir,"gene-level-cnv-consensus-annotated-mut-freq.tsv.gz"))
tpm_gene_zscore <- read_tsv(file.path(input_dir,"long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz"))
tpm_group_zscore <- read_tsv(file.path(input_dir,"long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz"))


# Replace EFO code for Wilms tumor
fusion_freq <- fusion_freq %>%
  mutate(diseaseFromSourceMappedId = replace(diseaseFromSourceMappedId, diseaseFromSourceMappedId == "MONDO_0006058", "MONDO_0019004"))

fused_gene_freq <- fused_gene_freq %>%
  mutate(diseaseFromSourceMappedId = replace(diseaseFromSourceMappedId, diseaseFromSourceMappedId == "MONDO_0006058", "MONDO_0019004"))

var_snv_freq <- var_snv_freq %>%
  mutate(diseaseFromSourceMappedId = replace(diseaseFromSourceMappedId, diseaseFromSourceMappedId == "MONDO_0006058", "MONDO_0019004"))

gene_snv_freq <- gene_snv_freq %>%
  mutate(diseaseFromSourceMappedId = replace(diseaseFromSourceMappedId, diseaseFromSourceMappedId == "MONDO_0006058", "MONDO_0019004"))

gene_cnv_freq <- gene_cnv_freq %>%
  mutate(diseaseFromSourceMappedId = replace(diseaseFromSourceMappedId, diseaseFromSourceMappedId == "MONDO_0006058", "MONDO_0019004"))

tpm_gene_zscore <- tpm_gene_zscore %>%
  mutate(MONDO = replace(MONDO, MONDO == "MONDO_0006058", "MONDO_0019004"))

tpm_group_zscore <- tpm_group_zscore %>%
  mutate(MONDO = replace(MONDO, MONDO == "MONDO_0006058", "MONDO_0019004"))



# Remove all DGD samples
fusion_freq <- fusion_freq %>%
  filter(!Dataset == "CHOP P30 Panel")

fused_gene_freq <- fused_gene_freq %>%
  filter(!Dataset == "CHOP P30 Panel")

var_snv_freq <- var_snv_freq %>%
  filter(!Dataset == "CHOP P30 Panel")

gene_snv_freq <- gene_snv_freq %>%
  filter(!Dataset == "CHOP P30 Panel")

gene_cnv_freq <- gene_cnv_freq %>%
  filter(!Dataset == "CHOP P30 Panel")

tpm_gene_zscore <- tpm_gene_zscore %>%
  filter(!cohort == "CHOP P30 Panel")

tpm_group_zscore <- tpm_group_zscore %>%
  filter(!cohort == "CHOP P30 Panel")




# Validate if the replaced code value is not for Wilms tumor
stopifnot(identical("Wilms tumor", fusion_freq %>% filter(diseaseFromSourceMappedId=="MONDO_0019004") %>% pull(Disease) %>% unique()))
stopifnot(identical("Wilms tumor", fused_gene_freq %>% filter(diseaseFromSourceMappedId=="MONDO_0019004") %>% pull(Disease) %>% unique()))
stopifnot(identical("Wilms tumor", var_snv_freq %>% filter(diseaseFromSourceMappedId=="MONDO_0019004") %>% pull(Disease) %>% unique()))
stopifnot(identical("Wilms tumor", gene_snv_freq %>% filter(diseaseFromSourceMappedId=="MONDO_0019004") %>% pull(Disease) %>% unique()))
stopifnot(identical("Wilms tumor", gene_cnv_freq %>% filter(diseaseFromSourceMappedId=="MONDO_0019004") %>% pull(Disease) %>% unique()))
stopifnot(identical("Wilms tumor", tpm_gene_zscore %>% filter(MONDO=="MONDO_0019004") %>% pull(Disease) %>% unique()))
stopifnot(identical("Wilms tumor", tpm_group_zscore %>% filter(MONDO=="MONDO_0019004") %>% pull(Disease) %>% unique()))

# Validate if the removed Dataset still exists
# Note: In some datasets, cohort="CHOP P30 Panel" may be labelled as cohort="DGD" which also needs to be removed
stopifnot(identical(FALSE,any(fusion_freq$Dataset=="CHOP P30 Panel")))
stopifnot(identical(FALSE,any(fused_gene_freq$Dataset=="CHOP P30 Panel")))
stopifnot(identical(FALSE,any(var_snv_freq$Dataset=="CHOP P30 Panel")))
stopifnot(identical(FALSE,any(gene_snv_freq$Dataset=="CHOP P30 Panel")))
stopifnot(identical(FALSE,any(gene_cnv_freq$Dataset=="CHOP P30 Panel")))
stopifnot(identical(FALSE,any(tpm_gene_zscore$cohort=="CHOP P30 Panel")))
stopifnot(identical(FALSE,any(tpm_group_zscore$cohort=="CHOP P30 Panel")))
