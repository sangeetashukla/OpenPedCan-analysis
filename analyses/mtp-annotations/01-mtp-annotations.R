library(purrr)
library(tidyverse)
library(jsonlite)
library(data.table)

#Load in all the csv files and add them to a dataframe
target_data_all <- list.files(path = "~/volume/OpenPedCan-analysis/scratch/mtp-csv/targets/",  # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  bind_rows                                         # Combine data sets into one data set 

#Pulls out the hugo gene symbol, ensemble gene id, and ensemble transcript id and renames the cols
target_formated_data_all <- target_data_all %>%
  select(approvedSymbol, id, transcriptIds) %>% 
  mutate(transcriptIds = gsub("\\[|\\]", "", transcriptIds)) %>% 
  separate_rows(transcriptIds) %>% 
  rename(gene_symbol = approvedSymbol, gene_id = id, transcript_id = transcriptIds) %>% distinct()

#writes a tsv file to the mtp-annotation results folder
write.table(target_formated_data_all, "~/volume/OpenPedCan-analysis/analyses/mtp-annotations/results/mtp-target-mapping.tsv", sep = "\t")

#Load in all the csv files and add them to a dataframe
print("Loading in disease csv files....")
disease_data_all <- list.files(path = "~/volume/OpenPedCan-analysis/scratch/mtp-csv/diseases",  # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>%
  lapply(read_csv) %>%                              # Store all files in list
  bind_rows                                         # Combine data sets into one data set


#Pulls out the EFO and Mondo disease identifiers and renames the cols
disease_formated_data_all <- disease_data_all %>%
  select(id,name,description,dbXRefs) %>%
  rename(id = id,name= name,description = description, dbXRefs = dbXRefs) %>% distinct()

#Write a tsv file to the results folder
write.table(disease_formated_data_all, "~/volume/OpenPedCan-analysis/analyses/mtp-annotations/results/mtp-disease-mapping.tsv", sep = "\t")
