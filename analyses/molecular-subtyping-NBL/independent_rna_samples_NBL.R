independent_rna_samples_NBL <- function(independent_dna_sample_df, histology_df){
  
#### ***** START ********* ####
histology_df_processed <- hist_filtered_df %>%mutate(match_id = paste(Kids_First_Participant_ID, sample_id, sep = "_"))
  
  
independent_dna <- histology_df_processed %>%filter(Kids_First_Participant_ID %in% independent_dna_samples_df$Kids_First_Participant_ID)
  
matched_rna <- histology_df_processed %>%filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                                                  RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                                                  match_id %in% independent_dna$match_id)
#################################################################################
  
# primary_descs <- c("Initial CNS Tumor", "Primary Tumor")
#   
#   
# only_rna_primary <- histology_df_processed %>% 
#     # keep rna from histology_df
#     dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
#                   RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
#                   tumor_descriptor %in% primary_descs,
#                   # find and remove participants which have 
#                   # matching dna samples in independent_wgswxspanel
#                   !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
#   only_rna_plus <- histology_df_processed %>% 
#     # keep rna from histology_df
#     dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
#                   RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
#                   # find and remove participants which have 
#                   # matching dna samples in independent_wgswxspanel
#                   !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
#                   # and participant not in only_rna_primary sample set
#                   !Kids_First_Participant_ID %in% only_rna_primary$Kids_First_Participant_ID
#     )
#   # has rna samples which match the independent samples provided plus rna only sample which are primary tumors plus rna samples where no primary primaries exists
#   sample_df <- bind_rows(matched_rna, only_rna_primary, only_rna_plus)
#   
#   
  
  #################################################################################
  # cohort_cancer_group_combo <- sample_df %>%
  #   dplyr::select(cohort, cancer_group) %>%
  #   unique() 
  # independent_each <- data.frame(Kids_First_Participant_ID = character(), 
  #                                Kids_First_Biospecimen_ID = character(), 
  #                                stringsAsFactors = FALSE)
  # 
  # # loop through each cohort and cancer group and combine the results together
  # for(i in 1:nrow(cohort_cancer_group_combo)){
  #   # filter to the specific cancer group and cohort
  #   cohort_name <- cohort_cancer_group_combo[i,1] %>% as.character()
  #   cancer_group_name <- cohort_cancer_group_combo[i,2] %>% as.character()
  #   
  #   # deal with cancer group is NA to avoid missing samples
  #   if(is.na(cancer_group_name)){
  #     filtered_df <- sample_df %>% 
  #       dplyr::filter(cohort == cohort_name) %>%
  #       dplyr::filter(is.na(cancer_group))
  #   }else{
  #     filtered_df <- sample_df %>% 
  #       dplyr::filter(cohort == cohort_name) %>%
  #       dplyr::filter(cancer_group == cancer_group_name)
  #   }
  #   
  #   # split filtered_df for each cohort into rnaseq and panel sample and rbind 
  #   # with panel samples at the end to allow dplyr::distinct to preferentially select  
  #   # rnaseq samples whenever there are samples for a participant in both rna libraries
  #   rnaseq <- filtered_df %>% filter(RNA_library != "exome_capture")
  #   exome_capture <- filtered_df %>% filter(RNA_library == "exome_capture")
  #   filtered_df <- rbind(rnaseq, exome_capture)
  #   
  #   # find the independent samples for the specific cancer group and cohort
  #   # "If there are multiple rows for a given combination of inputs, only the first
  #   # row will be preserved. If omitted, will use all variables." -- distinct in dplyr 0.8.3
  #   independent_filtered <- filtered_df %>%
  #     dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>%
  #     dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, 
  #                   cohort, cancer_group, experimental_strategy, tumor_descriptor)
  #   
  #   
  #   # merge the independent samples together
  #   independent_each <- rbind(independent_each, independent_filtered)
  # }
  # #### ***** EOF   ********* ####
  #return(list(independent_each,sample_df,matched_rna))
  return(matched_rna)
}