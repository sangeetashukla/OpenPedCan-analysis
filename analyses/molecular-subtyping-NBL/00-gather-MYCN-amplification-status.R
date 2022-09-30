library(dplyr)
library(data.table)

# Load the Files in ticket
consensus_wgs_plus_cnvkit_wxs_df<-as.data.frame(fread("/home/rstudio/OpenPedCan-analysis/data/consensus_wgs_plus_cnvkit_wxs.tsv.gz"))
