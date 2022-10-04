independent_rna_samples <- function(independent_dna_sample_df = NULL, 
                                    histology_df,
                                    independent_level = c("each-cohort", "all-cohorts", "all-cohorts-pre-release"),
                                    match_type = c("independent_dna", "independent_dna_plus_only_rna", "none"),
                                    tumor_description_rna_only = c("primary", "relapse", "primary_plus"),
                                    seed){
  match_type <- match.arg(match_type)
  tumor_description_rna_only <- match.arg(tumor_description_rna_only)
  independent_level <- match.arg(independent_level)
  if(!missing(seed)){set.seed(seed)}
  primary_descs <- c("Initial CNS Tumor", "Primary Tumor")
  relapse_descs <- c("Recurrence", "Progressive", "Progressive Disease Post-Mortem")
  
  # add dna-rna match_id column comprising of "Kids_First_Participant_ID" and "sample_id" to histology_df
  # because some Kids_First_Participant_ID might have more than one sample_id
  histology_df <- histology_df %>%
    dplyr::mutate(match_id = paste(Kids_First_Participant_ID, sample_id, sep = "_"))
  
  # Find sample set for the dna independent samples 
  # This will always be the included since in both the following
  # conditions "independent_dna" "independent_dna_plus_only_rna"
  # 
  independent_dna <- histology_df %>%
    # include matched independent_dna samples
    dplyr::filter(Kids_First_Participant_ID %in%
                    independent_dna_sample_df$Kids_First_Participant_ID)
  matched_rna <- histology_df %>%
    # keep rna from histology_df
    dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                  # exome capture and RNA-Seq
                  RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                  # keep match_id that only present in independent_dna
                  match_id %in% independent_dna$match_id)
  
  matched_rna_primary <- histology_df %>%
    # keep rna from histology_df
    dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                  RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                  tumor_descriptor %in% primary_descs,
                  match_id %in% independent_dna$match_id)
  
  matched_rna_relapse <- histology_df %>%
    # keep rna from histology_df
    dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                  RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                  tumor_descriptor %in% relapse_descs,
                  match_id %in% independent_dna$match_id)
  
  # Here we are adding only initial only-RNA-Seq samples
  # since this will always be part of independent_dna_plus_only_rna
  # regardless tumor_description_rna_only is "primary" OR "primary_plus"
  #
  if( match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "primary") {
    # find sample set where we initial only-RNA-Seq samples
    only_rna_primary <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                    RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                    tumor_descriptor %in% primary_descs,
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors
    sample_df <- bind_rows(matched_rna_primary, only_rna_primary)
  }
  
  if( match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "relapse") {
    # find sample set where we initial only-RNA-Seq samples
    only_rna_relapse <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                    RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                    tumor_descriptor %in% relapse_descs,
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors
    sample_df <- bind_rows(matched_rna_relapse,only_rna_relapse)
  }
  
  # Here we are adding only-RNA-Seq samples which are not initial
  # if tumor_description_rna_only == "primary_plus"
  #
  if(match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "primary_plus"){
    # find sample set where we only find rna samples
    only_rna_primary <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                    RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                    tumor_descriptor %in% primary_descs,
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    only_rna_plus <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq" | experimental_strategy == "Targeted Sequencing",
                    RNA_library %in% c("exome_capture", "stranded", "poly-A", "poly-A stranded"),
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                    # and participant not in only_rna_primary sample set
                    !Kids_First_Participant_ID %in% only_rna_primary$Kids_First_Participant_ID
      )
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors plus rna samples where no primary primaries exists
    sample_df <- bind_rows(matched_rna, only_rna_primary, only_rna_plus)
  } 
  
  if(independent_level == "each-cohort"){
    # find out the list of cohorts and cancer groups and loop through them to run the process for each cancer group and cohort
    # combination of cohort and cancer groups
    cohort_cancer_group_combo <- sample_df %>%
      dplyr::select(cohort, cancer_group) %>%
      unique() 
    
    independent_each <- data.frame(Kids_First_Participant_ID = character(), 
                                   Kids_First_Biospecimen_ID = character(), 
                                   stringsAsFactors = FALSE)
    
    # loop through each cohort and cancer group and combine the results together
    for(i in 1:nrow(cohort_cancer_group_combo)){
      # filter to the specific cancer group and cohort
      cohort_name <- cohort_cancer_group_combo[i,1] %>% as.character()
      cancer_group_name <- cohort_cancer_group_combo[i,2] %>% as.character()
      
      # deal with cancer group is NA to avoid missing samples
      if(is.na(cancer_group_name)){
        filtered_df <- sample_df %>% 
          dplyr::filter(cohort == cohort_name) %>%
          dplyr::filter(is.na(cancer_group))
      }else{
        filtered_df <- sample_df %>% 
          dplyr::filter(cohort == cohort_name) %>%
          dplyr::filter(cancer_group == cancer_group_name)
      }
      
      # split filtered_df for each cohort into rnaseq and panel sample and rbind 
      # with panel samples at the end to allow dplyr::distinct to preferentially select  
      # rnaseq samples whenever there are samples for a participant in both rna libraries
      rnaseq <- filtered_df %>% filter(RNA_library != "exome_capture")
      exome_capture <- filtered_df %>% filter(RNA_library == "exome_capture")
      filtered_df <- rbind(rnaseq, exome_capture)
      
      # find the independent samples for the specific cancer group and cohort
      # "If there are multiple rows for a given combination of inputs, only the first
      # row will be preserved. If omitted, will use all variables." -- distinct in dplyr 0.8.3
      independent_filtered <- filtered_df %>%
        dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>%
        dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, 
                      cohort, cancer_group, experimental_strategy, tumor_descriptor)
      
      
      # merge the independent samples together
      independent_each <- rbind(independent_each, independent_filtered)
    }
    return(independent_each)
  }
  
  if(independent_level == "all-cohorts"){
    rnaseq <- sample_df %>% filter(RNA_library != "exome_capture")
    exome_capture <- sample_df %>% filter(RNA_library == "exome_capture")
    sample_df <- rbind(rnaseq, exome_capture)
    
    independent_all <- sample_df %>%
      dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>%
      dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, 
                    cohort, cancer_group, experimental_strategy, tumor_descriptor)
    
    return(independent_all)
  }
  
  if(match_type == "none" & independent_level == "all-cohorts-pre-release"){
    if(tumor_description_rna_only %in% c("primary_plus")){
      # find cases where non-primary is the only option
      no_primary <- histology_df %>% 
        dplyr::group_by(Kids_First_Participant_ID) %>%
        dplyr::summarize(n_primary = sum(tumor_descriptor %in% primary_descs)) %>%
        dplyr::filter(n_primary == 0) %>%
        dplyr::pull(Kids_First_Participant_ID)
    } else {
      no_primary <- c()
    }
    
    if(tumor_description_rna_only %in% c("primary", "primary_plus")){
      primary_df <- histology_df %>%
        dplyr::filter(tumor_descriptor %in% primary_descs)
      noprimary_df <- histology_df %>%
        dplyr::filter(Kids_First_Participant_ID %in% no_primary)
      sample_df <- dplyr::bind_rows(primary_df, noprimary_df)
    } 
    
    if(tumor_description_rna_only == "primary"){
      sample_df <- histology_df %>%
        dplyr::filter(tumor_descriptor %in% primary_descs)
    } 
    
    if(tumor_description_rna_only == "relapse"){
      sample_df <- histology_df %>%
        dplyr::filter(tumor_descriptor %in% relapse_descs)
    } 
    
    independent_all <- sample_df %>%
      dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>%
      dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, cohort, tumor_descriptor) %>%
      arrange(Kids_First_Biospecimen_ID)
    
    return(independent_all)
  }
}