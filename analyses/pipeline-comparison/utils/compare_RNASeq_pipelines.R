# Invoke libraries
suppressPackageStartupMessages(
  {
    library(tidyr)
    library(dplyr)
  }
)

# This function performs quality check on expected RNA-Seq counts input  matrices for which correlation coefficient can be calculated
# Input: 2 matrices with RNA-seq counts processed by different pipelines
# Output: R list object containing 2 matrices
# Output_note: Each result matrix can be accessed as calc_cor_input[[mat1]]
# Output_note: , calc_cor_input[[mat2]]
qc_matrix_input <- function(matrix1, matrix2){
  #rownames must refer to gene names
  common_genes <- intersect(rownames(matrix1),rownames(matrix2)) %>%
    sort()
  mat1_all_genes <- matrix1[common_genes,]
  mat2_all_genes <- matrix2[common_genes,]
  
  #colnames must refer to sample ids
  common_samples <- intersect(colnames(mat1_all_genes), colnames(mat2_all_genes)) %>%
    sort()
  mat1_final <- mat1_all_genes[, common_samples]
  mat2_final <- mat2_all_genes[, common_samples]
  
  #Validate all rows/genes and all columns/sample_ids match across matrices
  mat1_col_extra <- setdiff(colnames(mat1_final),colnames(mat2_final))
  # character(0)
  mat2_col_extra <- setdiff(colnames(mat2_final),colnames(mat1_final))
  # character(0)
  
  calc_cor_input <- list(mat1_final,mat2_final)
  return(calc_cor_input)
  
}


#
calculate_matrix_cor <- function(matrix1, matrix2, wf1_name, wf2_name){
  
  result_compare_all_genes <- data.frame()
  for(i in 1:nrow(matrix1)){
  
    p1_p2_coef = cor(matrix1[i,],matrix2[i,], method = "pearson", use = "complete.obs")
    result_compare_all_genes <- rbind(result_compare_all_genes, p1_p2_coef)
    
  }
  
  row_names <- rownames(matrix1)
  rownames(result_compare_all_genes) <- row_names
  colnames(result_compare_all_genes) <- "coef_p1_p2"
  na.omit(result_compare_all_genes)
  
  
  serr_r <- sapply(result_compare_all_genes, function(r) (1 - (r^2))/length(colnames(result_compare_all_genes)))
  perr_r <- sapply(result_compare_all_genes, function(r) 0.6745*((1 - (r^2))/sqrt(length(colnames(result_compare_all_genes)))))
  na.omit(result_compare_all_genes)
  
  final_result <- data.frame()
  final_result <- cbind(result_compare_all_genes, serr_r, perr_r)
  #print(c(paste("coeff",wf1_name,wf2_name,sep = "_"), "stderr_coef", "proberr_coef"))
  colnames(final_result) <- c(paste("coeff",wf1_name,wf2_name,sep = "_"), "stderr_coef", "proberr_coef")
  #print(class(final_result))
  return(final_result)
  
}