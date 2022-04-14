# script to visualize the PCA results from 01-full_dataset_compute_pca_umap_counts.R

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(plotly)  
  library(ggpubr)
})

# generalized function to call for plotting
plot_data <- function(pca_output, title = "", color_var, shape_var, label_var = NULL){
  # % variability
  var_explained = pca_output$sdev^2 / sum(pca_output$sdev^2)
  
  # first two components explain majority of the variation 
  pca_output_for_plot <- pca_output$rotation
  pca_output_for_plot <- data.frame(pca_output_for_plot)[1:2]
  pca_output_for_plot <- readr::read_tsv('../../data/histologies.tsv') %>%
    dplyr::select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id, sample_type, cohort, RNA_library, cancer_group, gtex_group, composition) %>%
    inner_join(pca_output_for_plot %>%
                 rownames_to_column('Kids_First_Biospecimen_ID')) %>% 
    mutate(label = paste0(cohort, "_", RNA_library),
           label2 = paste0(cohort, "_", composition),
           group = ifelse(cohort == "GTEx", gtex_group, cancer_group))
  
  # plot
  plot_out <- ggplot(pca_output_for_plot, 
                     aes(PC1, PC2, color = get(color_var),
                         shape = get(shape_var),
                         text = paste0("Sample type:", sample_type,
                                       "\nGroup:", group,
                                       "\nID:", sample_id))) +
    geom_point(size = 2, alpha = 0.5) +
    xlab(paste0("PC1: ", round(var_explained[1]*100, 2), "%")) +
    ylab(paste0("PC2: ", round(var_explained[2]*100, 2), "%")) +
    ggtitle(title) +
    ggpubr::theme_pubr() +
    guides(color=guide_legend(title=color_var),
           shape=guide_legend(title=shape_var))
  
  if(!is.null(label_var)){
    plot_out <- plot_out +
      ggrepel::geom_text_repel(aes(label = get(label_var)), size = 2, show.legend = FALSE) 
  }
  
  return(plot_out)
}

# output directory
output_dir <- file.path('output', 'QC_clustering')
dir.create(output_dir, showWarnings = F, recursive = T)

# 1. TUMORS + NORMALS
# 1.1. log2 + highly variable genes > 90% quantile
pca_tumors_normals_top_var <- readRDS(file.path(output_dir, 'pca_output_tumors_normals.rds'))
p1 <- plot_data(pca_output = pca_tumors_normals_top_var, 
                title = paste0("Tumors + Normals: Top variable genes (> 90% quantile) (n = ", nrow(pca_tumors_normals_top_var$x), ")"), 
                color_var = "cohort", 
                shape_var = "composition") + theme(legend.position = "right")
p2 <- plot_data(pca_output = pca_tumors_normals_top_var, 
                title = paste0("Tumors + Normals: Top variable genes (> 90% quantile) (n = ", nrow(pca_tumors_normals_top_var$x), ")"), 
                color_var = "cohort", 
                shape_var = "RNA_library") + theme(legend.position = "right")

# plot_data(pca_output = pca_tumors_normals_top_var, title = "", color_var = "cohort", shape_var = "composition") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors + Normals: Top variable genes (> 90% quantile)",
#          legend = list(title = list(text='Cohort_Composition\n')))
# plot_data(pca_output = pca_tumors_normals_top_var, title = "", color_var = "cohort", shape_var = "RNA_library") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors + Normals: Top variable genes (> 90% quantile)",
#          legend = list(title = list(text='Cohort_Library\n')))

# 1.2. log2 + filter house keeping genes
pca_tumors_normals_hk <- readRDS(file.path(output_dir, 'pca_output_tumors_normals_hk.rds'))
q1 <- plot_data(pca_output = pca_tumors_normals_hk, 
                title = paste0("Tumors + Normals: House keeping genes (n = ", nrow(pca_tumors_normals_hk$x), ")"), 
                color_var = "cohort", 
                shape_var = "composition") + theme(legend.position = "right")
q2 <- plot_data(pca_output = pca_tumors_normals_hk, 
                title = paste0("Tumors + Normals: House keeping genes (n = ", nrow(pca_tumors_normals_hk$x), ")"), 
                color_var = "cohort", 
                shape_var = "RNA_library") + theme(legend.position = "right")

# plot_data(pca_output = pca_tumors_normals_hk, title = "", color_var = "cohort", shape_var = "composition") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors + Normals: House keeping genes (n = 25)",
#          legend = list(title = list(text='Cohort_Composition\n')))
# plot_data(pca_output = pca_tumors_normals_hk, title = "", color_var = "cohort", shape_var = "RNA_library") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors + Normals: House keeping genes (n = 25)",
#          legend = list(title = list(text='Cohort_Library\n')))

# combine both plots
# by composition type
tumors_normals_combined_plot <- ggarrange(p1, q1, common.legend = T, legend = "bottom")
ggsave(tumors_normals_combined_plot, 
       filename = file.path(output_dir, 'pca_output_tumors_normals_composition_combined.pdf'), 
       width = 15, height = 6)

# by library type
tumors_normals_combined_plot <- ggarrange(p2, q2, common.legend = T, legend = "bottom")
ggsave(tumors_normals_combined_plot, 
       filename = file.path(output_dir, 'pca_output_tumors_normals_library_combined.pdf'), 
       width = 15, height = 6)

# 2. TUMORS only
# 2.1. log2 + highly variable genes > 90% quantile
pca_tumors_top_var <- readRDS(file.path(output_dir, 'pca_output_tumors.rds'))
p1 <- plot_data(pca_output = pca_tumors_top_var, 
                title = paste0("Tumors: Top variable genes (> 90% quantile) (n = ", nrow(pca_tumors_top_var$x), ")"), 
                color_var = "cohort", 
                shape_var = "composition") + theme(legend.position = "right")
p2 <- plot_data(pca_output = pca_tumors_top_var, 
                title = paste0("Tumors: Top variable genes (> 90% quantile) (n = ", nrow(pca_tumors_top_var$x), ")"), 
                color_var = "cohort", 
                shape_var = "RNA_library") + theme(legend.position = "right")

# plot_data(pca_output = pca_tumors_top_var, title = "", color_var = "cohort", shape_var = "composition") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors: Top variable genes (> 90% quantile)",
#          legend = list(title = list(text='Cohort_Composition\n')))
# plot_data(pca_output = pca_tumors_top_var, title = "", color_var = "cohort", shape_var = "RNA_library") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors: Top variable genes (> 90% quantile)",
#          legend = list(title = list(text='Cohort_Library\n')))

# 2.2. log2 + filter house keeping genes
pca_tumors_hk <- readRDS(file.path(output_dir, 'pca_output_tumors_hk.rds'))
q1 <- plot_data(pca_output = pca_tumors_hk, 
                title = paste0("Tumors: House keeping genes (n = ", nrow(pca_tumors_hk$x), ")"), 
                color_var = "cohort", 
                shape_var = "composition") + theme(legend.position = "right")
q2 <- plot_data(pca_output = pca_tumors_hk, 
                title = paste0("Tumors: House keeping genes (n = ", nrow(pca_tumors_hk$x), ")"), 
                color_var = "cohort", 
                shape_var = "RNA_library") + theme(legend.position = "right")

# plot_data(pca_output = pca_tumors_hk, title = "", color_var = "cohort", shape_var = "composition") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors: House keeping genes (n = 25)",
#          legend = list(title = list(text='Cohort_Composition\n')))
# plot_data(pca_output = pca_tumors_hk, title = "", color_var = "cohort", shape_var = "RNA_library") %>%
#   plotly::ggplotly() %>%
#   layout(title = "Tumors: House keeping genes (n = 25)",
#          legend = list(title = list(text='Cohort_Library\n')))

# combine both plots
# by composition type
tumors_combined_plot <- ggarrange(p1, q1, common.legend = T, legend = "bottom")
ggsave(tumors_combined_plot, 
       filename = file.path(output_dir, 'pca_output_tumors_composition_combined.pdf'), 
       width = 15, height = 6)

# by library type
tumors_combined_plot <- ggarrange(p2, q2, common.legend = T, legend = "bottom")
ggsave(tumors_combined_plot, 
       filename = file.path(output_dir, 'pca_output_tumors_library_combined.pdf'), 
       width = 15, height = 6)

# 3. TUMORS matched samples
# 3.1. log2 + highly variable genes > 90% quantile
pca_matched_tumors_top_var <- readRDS(file.path(output_dir, 'pca_output_matched_samples.rds'))
p1 <- plot_data(pca_output = pca_matched_tumors_top_var, 
                title = paste0("Matched tumor samples: Top variable genes (> 90% quantile) (n =", nrow(pca_matched_tumors_top_var$x),")"), 
                color_var = "cohort", shape_var = "composition", label_var = "Kids_First_Participant_ID") +
  theme_pubr(legend = "right")
p2 <- plot_data(pca_output = pca_matched_tumors_top_var, 
                title = paste0("Matched tumor samples: Top variable genes (> 90% quantile) (n =", nrow(pca_matched_tumors_top_var$x),")"), 
                color_var = "cohort", shape_var = "RNA_library", label_var = "Kids_First_Participant_ID") +
  theme_pubr(legend = "right")

# 3.2. log2 + filter house keeping genes
pca_matched_tumors_hk <- readRDS(file.path(output_dir, 'pca_output_matched_samples_hk.rds'))
q1 <- plot_data(pca_output = pca_matched_tumors_hk, 
                title = paste0("Matched tumor samples: House keeping genes (n =", nrow(pca_matched_tumors_hk$x), ")"), 
                color_var = "cohort", shape_var = "composition", label_var = "Kids_First_Participant_ID") +
  theme_pubr(legend = "right")
q2 <- plot_data(pca_output = pca_matched_tumors_hk, 
                title = paste0("Matched tumor samples: House keeping genes (n =", nrow(pca_matched_tumors_hk$x), ")"), 
                color_var = "cohort", shape_var = "RNA_library", label_var = "Kids_First_Participant_ID") +
  theme_pubr(legend = "right")

# combine both plots
# by composition type
match_samples_combined_plot <- ggarrange(p1, q1, common.legend = T, legend = "bottom")
ggsave(match_samples_combined_plot, 
       filename = file.path(output_dir, 'pca_output_matched_samples_composition_combined.pdf'), 
       width = 15, height = 6)

# by library type
match_samples_combined_plot <- ggarrange(p2, q2, common.legend = T, legend = "bottom")
ggsave(match_samples_combined_plot, 
       filename = file.path(output_dir, 'pca_output_matched_samples_library_combined.pdf'), 
       width = 15, height = 6)

# Seems like tissue type is the primary driver of expression differences (solid vs non-solid tissues)
# so the major differences are biological rather than technical
