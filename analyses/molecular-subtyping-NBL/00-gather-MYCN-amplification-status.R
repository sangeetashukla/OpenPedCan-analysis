setwd("~/OpenPedCan-analysis/analyses/molecular-subtyping-NBL")
library(dplyr)
library(data.table)

# Load the files from the ticket
consensus_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/consensus_wgs_plus_cnvkit_wxs.tsv.gz"))
hist_df<- as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/histologies.tsv"))
cnv_cnvkit_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/cnv-cnvkit.seg.gz"))
cnv_controlfreec_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/cnv-controlfreec.tsv.gz"))
gene_expression_df<-as.data.frame(readRDS("/home/rstudio/OpenPedCan-analysis/data/gene-expression-rsem-tpm-collapsed.rds"))






# For each NBL sample, gather MYCN amplification status 

## Step 1: Filter the neuroblastoma, ganglioneuroblastoma, and ganglioneuroma samples from hist file first
hist_df <- hist_df %>% filter(pathology_diagnosis=="Neuroblastoma"|
                                pathology_diagnosis=="Ganglioneuroblastoma"|
                                pathology_diagnosis=="Ganglioneuroblastoma,nodular"|
                                pathology_diagnosis=="Ganglioneuroblastoma, intermixed"|
                                pathology_diagnosis=="Ganglioneuroma, maturing subtype OR Ganglioneuroblastoma, well differentiated")


## Step 2: Filter the MCYN records from consensus_wgs_plus_cnvkit_wxs.tsv.gz
consensus_df <- consensus_df %>% filter(gene_symbol=="MYCN")

## Join the two data frames to get MYCN details with respect to the diagnosis
Myc_df <- inner_join(hist_df,consensus_df,by= c("Kids_First_Biospecimen_ID" = "biospecimen_id"))
