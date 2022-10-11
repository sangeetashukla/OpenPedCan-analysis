setwd("~/OpenPedCan-analysis/analyses/molecular-subtyping-NBL")
source('~/OpenPedCan-analysis/analyses/molecular-subtyping-NBL/independent_rna_samples_NBL.R')
library(dplyr)
library(data.table)

# Load the files from the ticket
consensus_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/consensus_wgs_plus_cnvkit_wxs.tsv.gz"))
hist_df<- as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/histologies.tsv"))
cnv_cnvkit_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/cnv-cnvkit.seg.gz"))
cnv_controlfreec_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/cnv-controlfreec.tsv.gz"))
gene_expression_df<-as.data.frame(readRDS("/home/rstudio/OpenPedCan-analysis/data/gene-expression-rsem-tpm-collapsed.rds"))

# should be like independent-specimens.wgswxspanel.primary.
#test_df <- as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/independent-specimens.wgswxspanel.primary.prefer.wgs.tsv"))



# For each NBL sample, gather MYCN amplification status 

## Step 1: Filter the neuroblastoma, ganglioneuroblastoma, and ganglioneuroma samples from hist file first
hist_filtered_df <- hist_df %>% filter(pathology_diagnosis=="Neuroblastoma"|
                                pathology_diagnosis=="Ganglioneuroblastoma"|
                                pathology_diagnosis=="Ganglioneuroblastoma,nodular"|
                                pathology_diagnosis=="Ganglioneuroblastoma, intermixed"|
                                pathology_diagnosis=="Ganglioneuroma, maturing subtype OR Ganglioneuroblastoma, well differentiated")


## Step 2: Filter the MCYN records from consensus_wgs_plus_cnvkit_wxs.tsv.gz
consensus_filtered_df <- consensus_df %>% filter(gene_symbol=="MYCN")

## Join the two data frames to get MYCN details with respect to the diagnosis
Myc_df <- inner_join(hist_filtered_df,consensus_filtered_df,by= c("Kids_First_Biospecimen_ID" = "biospecimen_id"))


# Create alterations/output table for matched DNA/RNA biospecimens. Match ID can be adapted from the independent-samples module

## First create a df like A data frame of samples, with columnscorresponding to those in `independent-specimens.wgswxspanel.primary.tsv
# independent-specimens.wgswxspanel.relapse.tsv independent-specimens.wgswxspanel.primary-plus.tsv
# columns that I need here : "Kids_First_Participant_ID" "Kids_First_Biospecimen_ID" "cohort" "cancer_group" "experimental_strategy""tumor_descriptor"

independent_dna_samples_df <- Myc_df %>% select(Kids_First_Participant_ID,Kids_First_Biospecimen_ID,
                                             cohort,cancer_group,experimental_strategy,tumor_descriptor)

matched_independent_NBL_df <- independent_rna_samples_NBL(independent_dna_samples_df,hist_filtered_df) ####

## Step 3: Add MYCN TPM and MYCN copy number and status to table. It may also be useful to see a quick ascending barplot of MYCN TPM (y-axis) 
## and sample (x-axis) colored by clinical status (see 4 below). 


