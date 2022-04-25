# function to evaluate RUVg: Estimating the factors of unwanted variation using control genes
# uses negative control genes, assumed to have constant expression across samples
ruvg_test <- function(seq_expr_set, emp_neg_ctrl_genes, k_val = 1:2, prefix, diff_type = c("deseq2", "edger"), output_dir){
  
  # create output directory
  if(diff_type == "edger"){
    output_dir <- file.path(output_dir, "edger_analysis")
  } else {
    output_dir <- file.path(output_dir, "deseq2_analysis")
  }
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # to estimate the factors of unwanted variation,
  # subset to negative control genes (genes that can be assumed not to be influenced by the covariates of interest)
  emp_neg_ctrl_genes <- intersect(emp_neg_ctrl_genes, rownames(seq_expr_set@assayData$counts))
  emp_neg_ctrl_genes <- emp_neg_ctrl_genes[emp_neg_ctrl_genes %in% rownames(seq_expr_set@assayData$counts)]
  print(paste("Number of negative control genes:", length(emp_neg_ctrl_genes)))
  
  # if no genes left in emp_neg_ctrl_genes, then don't run function
  if(length(emp_neg_ctrl_genes) == 0){
    return(NULL)
  } 
  
  # if length of emp_neg_ctrl_genes is less than k_val but greater than 0, then set max k_val number of genes
  if(length(emp_neg_ctrl_genes) < max(k_val)) {
    k_val <- c(1:length(emp_neg_ctrl_genes))
  }
  
  # loop through all Ks
  cluster_plot <- list()
  pval_hist_plot <- list()
  dge_output_neg_control_genes <- list()
  chisq_out <- list()
  ks_out <- list()
  for(i in 1:length(k_val)){
    print(paste0("K = ", i))
    print(paste0('Run differential gene expression with RUVSeq estimated batch effect on poly-A vs stranded RNA-seq...'))
    
    # run RUVg assuming there are k_val factors of unwanted variation
    ruvg_set <- RUVg(x = seq_expr_set, cIdx = emp_neg_ctrl_genes, k = k_val[i])
    
    # pca and umap after ruvg
    ruvg_pca <- edaseq_plot(object = ruvg_set, title = paste0("PCA: RUVg output (k = ", i, ")"), type = "PCA")
    ruvg_umap <- edaseq_plot(object = ruvg_set, title = paste0("UMAP: RUVg output (k = ", i, ")"), type = "UMAP")
    cluster_plot[[i]] <- ggpubr::ggarrange(ruvg_pca, ruvg_umap, common.legend = T, legend = "bottom")
    
    # differential expression after RUVg
    if(diff_type == "deseq2"){
      # W_i corresponds to the factors of "unwanted variation"
      # factor for unwanted variation comes first for deseq2
      design <- model.matrix(as.formula(paste0('~0 + W_', i, '+ patient_id + rna_library')), data = pData(ruvg_set))
      ruv_dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts(ruvg_set), colData = pData(ruvg_set), design = design)
      ruv_dds <- DESeq2::DESeq(ruv_dds)
      
      # Export results and save as CSV
      dge_output <- DESeq2::results(ruv_dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
      dge_output <- dge_output %>% 
        as.data.frame() %>% 
        rownames_to_column('gene') %>%
        arrange(padj) 
    } else if(diff_type == "edger") {
      # W corresponds to the factors of "unwanted variation"
      # factor for unwanted variation comes last for edgeR
      design <- model.matrix(as.formula(paste0('~0 + patient_id + rna_library +', 'W_', i)), data = pData(ruvg_set))
      y <- DGEList(counts = counts(ruvg_set), group = rna_library)
      y <- calcNormFactors(y, method = "upperquartile")
      y <- estimateGLMCommonDisp(y, design)
      y <- estimateGLMTagwiseDisp(y, design)
      fit <- glmFit(y, design)
      lrt <- glmLRT(fit, coef = grep('rna_library', colnames(fit$coefficients)))
      dge_output <- topTags(lrt, n = Inf)$table %>%
        rownames_to_column('gene') %>%
        dplyr::rename("pvalue" = "PValue", "padj" = "FDR")
    }
    
    # save dge output (commenting out to reduce output files)
    # filename <- file.path(output_dir, paste0(prefix, '_stranded_vs_polya_dge_ruvg_k', k_val[i], '_', diff_type, '_result.csv'))
    # write_tsv(dge_output, file = filename)
    
    # plot and save p-value histogram
    # evaluate the distribution of p-values for full transcriptome
    pval_hist_plot[[i]] <- deseq2_pvals_histogram(res_df = dge_output,
                                xlab = 'stranded vs poly-A RUVg p-value (negative control genes)',
                                ylab = 'Gene count',
                                title = paste0('Histogram of stranded vs poly-A paired analysis (k = ', i, ')'))
    
    # test for uniformity (negative control genes only)
    # chisq test for p-values 
    dge_output_neg_control_genes[[i]] <- dge_output %>%
      filter(gene %in% emp_neg_ctrl_genes)
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
  pdf(file = file.path(output_dir, paste0(prefix, '_stranded_vs_polya_dge_ruvg_', diff_type, '_clustering.pdf')), width = 6, height = 4)
  print(cluster_plot)
  dev.off()

  # p-value histogram (negative control genes only)
  pdf(file = file.path(output_dir, paste0(prefix, '_stranded_vs_polya_dge_ruvg_', diff_type, '_histogram.pdf')), width = 8, height = 7)
  print(pval_hist_plot)
  dev.off()
  
  # rbind and save chisq values
  data.table::rbindlist(chisq_out) %>%
    write_tsv(file = file.path(output_dir, paste0(prefix, '_stranded_vs_polya_dge_ruvg_', diff_type, '_chisq_pvalues.tsv')))
  
  # rbind and save ks values
  data.table::rbindlist(ks_out) %>%
    write_tsv(file = file.path(output_dir, paste0(prefix, '_stranded_vs_polya_dge_ruvg_', diff_type, '_ks_pvalues.tsv')))
}
