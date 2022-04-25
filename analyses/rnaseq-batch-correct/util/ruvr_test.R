# function to evaluate RUVr: Estimating the factors of unwanted variation using residuals
# NOTE: this function can only be run with edgeR
ruvr_test <- function(seq_expr_set, emp_neg_ctrl_genes, residuals, k_val = 1:2, output_dir){
  
  # if no genes left in emp_neg_ctrl_genes, then don't run function
  if(length(emp_neg_ctrl_genes) == 0){
    return(NULL)
  } 
  
  # if length of emp_neg_ctrl_genes is less than k_val but greater than 0, then set max k_val number of genes
  if(length(emp_neg_ctrl_genes) < max(k_val)) {
    k_val <- c(1:length(emp_neg_ctrl_genes))
  }
  
  # create output directory
  output_dir <- file.path(output_dir, "edger_analysis")
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # loop through all Ks
  cluster_plot <- list()
  pval_hist_plot <- list()
  dge_output_neg_control_genes <- list()
  chisq_out <- list()
  ks_out <- list()
  for(i in 1:length(k_val)){
    print(paste0("K = ", i))
    print(paste0('Run differential gene expression with RUVSeq estimated batch effect on poly-A vs stranded RNA-seq...'))
    
    # run RUVr assuming there are k_val factors of unwanted variation
    ruvr_set <- RUVr(x = seq_expr_set, cIdx = emp_neg_ctrl_genes, k = k_val[i], residuals = residuals)
    
    # pca and umap after ruvr
    ruvr_pca <- edaseq_plot(object = ruvr_set, title = paste0("PCA: RUVg output (k = ", i, ")"), type = "PCA")
    ruvr_umap <- edaseq_plot(object = ruvr_set, title = paste0("UMAP: RUVg output (k = ", i, ")"), type = "UMAP")
    cluster_plot[[i]] <- ggpubr::ggarrange(ruvr_pca, ruvr_umap, common.legend = T, legend = "bottom")
    
    # differential expression after RUVg
    # W corresponds to the factors of "unwanted variation"
    # factor for unwanted variation comes last for edgeR
    design <- model.matrix(as.formula(paste0('~0 + patient_id + rna_library +', 'W_', i)), data = pData(ruvr_set))
    y <- DGEList(counts = counts(ruvr_set), group = rna_library)
    y <- calcNormFactors(y, method = "upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = grep('rna_library', colnames(fit$coefficients)))
    dge_output <- topTags(lrt, n = Inf)$table %>%
      rownames_to_column('gene') %>%
      dplyr::rename("pvalue" = "PValue", "padj" = "FDR")
    
    # save dge output (commenting out to reduce output files)
    # filename <- file.path(output_dir, paste0('stranded_vs_polya_dge_ruvr_k', k_val[i], '_edger_result.csv'))
    # write_tsv(dge_output, file = filename)
    
    # plot and save p-value histogram
    # evaluate the distribution of p-values among negative control genes only
    dge_output_neg_control_genes[[i]] <- dge_output %>%
      filter(gene %in% emp_neg_ctrl_genes)
    pval_hist_plot[[i]] <- deseq2_pvals_histogram(res_df = dge_output_neg_control_genes[[i]],
                                                  xlab = 'stranded vs poly-A RUVr p-value (negative control genes)',
                                                  ylab = 'Gene count',
                                                  title = paste0('Histogram of stranded vs poly-A paired analysis (k = ', i, ')'))
    
    # chisq test for p-values (negative control genes only)
    chisq_out[[i]] <- chisq.test(x = dge_output_neg_control_genes[[i]]$pvalue)
    chisq_out[[i]] <- broom::tidy(chisq_out[[i]])
    chisq_out[[i]]$k <- k_val[i]
    
    # ks test for p-values
    ks_out[[i]] <- ks.test(x = dge_output_neg_control_genes[[i]]$pvalue, punif)
    ks_out[[i]] <- broom::tidy(ks_out[[i]])
    ks_out[[i]]$k <- k_val[i]
  }
  
  # save the plots for all k values in a multi-page pdf file
  # clustering output (PCA/UMAP)
  pdf(file = file.path(output_dir, 'stranded_vs_polya_dge_ruvr_edger_clustering.pdf'), width = 6, height = 4)
  print(cluster_plot)
  dev.off()
  
  # p-value histogram
  pdf(file = file.path(output_dir, 'stranded_vs_polya_dge_ruvr_edger_histogram.pdf'), width = 8, height = 7)
  print(pval_hist_plot)
  dev.off()
  
  # rbind and save chisq values
  data.table::rbindlist(chisq_out) %>%
    write_tsv(file = file.path(output_dir, 'stranded_vs_polya_dge_ruvr_edger_chisq_pvalues.tsv'))
  
  # rbind and save ks values
  data.table::rbindlist(ks_out) %>%
    write_tsv(file = file.path(output_dir, 'stranded_vs_polya_dge_ruvr_edger_ks_pvalues.tsv'))
}
