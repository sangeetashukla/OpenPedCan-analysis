#This script assigns ATRT into three known subtypes using methylation result.
#Subtypiong esults is saved as ATRT-molecular-subtypes.tsv

# Set up library
library(tidyverse)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to result
results_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "results")


# Read in histologies
histo <- 
  readr::read_tsv(file.path(root_dir, "data", "histologies-base.tsv"), guess_max = 10000)

# Filter histo, 
# select all ATRT biospecimens from PBTA and/DGD
atrt_df <- histo %>%
  dplyr::filter(short_histology == "ATRT",
                sample_type == "Tumor",
                composition == "Solid Tissue",
                cohort == "PBTA" | cohort == "DGD")

# Create a dataframe with sample_id and matched biospecimens id for RNA-Seq, WGS and methylation
atrt_df_meth <- atrt_df %>% 
  filter(experimental_strategy == "Methylation") %>%
  select(sample_id, Kids_First_Biospecimen_ID, 
         cns_methylation_subclass, cns_methylation_subclass_score)

atrt_df_WGS <- atrt_df %>% 
  filter(experimental_strategy == "WGS") %>%
  select(sample_id,  Kids_First_Biospecimen_ID)

atrt_df_RNA <- atrt_df %>% 
  filter(experimental_strategy == "RNA-Seq") %>%
  select(sample_id,  Kids_First_Biospecimen_ID) 

atrt_subtype <- atrt_df_meth %>%
  full_join(atrt_df_WGS, by = "sample_id", suffix = c("_meth", "_DNA")) %>% 
  full_join(atrt_df_RNA, by = "sample_id") %>% 
  rename(Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID)

# create a list for ATRT subtypes
ATRT_subtype_list <- c("ATRT, MYC", "ATRT, SHH", "ATRT, TYR")

# for the samples, whose cns_methlation_subclass_score >= 0.8 and cns_methylation_subclass is one of the three types in ATRT_subtype_list, their molecular subtype are same as cns_methylation_subclass
# for the samples without methylation sequencing, their molecular subtype are "ATRT, To be classified."
# For the samples fit all the other situations, their molecular subtype are "ATRT, To be classified."
atrt_subtype <- atrt_subtype %>% 
  mutate(molecular_subtype = case_when(cns_methylation_subclass_score >= 0.8 & cns_methylation_subclass %in% ATRT_subtype_list ~ cns_methylation_subclass, 
                                       is.na(Kids_First_Biospecimen_ID_meth) ~ "ATRT, To be classified",
                                       TRUE ~ "ATRT, To be classified")) %>%
  select(sample_id, Kids_First_Biospecimen_ID_meth, Kids_First_Biospecimen_ID_DNA, Kids_First_Biospecimen_ID_RNA, molecular_subtype) %>%
  
# write result
  readr::write_tsv(file.path(results_dir, "ATRT-molecular-subtypes.tsv"))

