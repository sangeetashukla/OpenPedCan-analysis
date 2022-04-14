deseq2_pvals_histogram <- function(res_df, xlab, ylab, title) {
  stopifnot(all(c('pvalue', 'padj') %in% colnames(res_df)))
  
  p <- ggplot(res_df, aes(x = pvalue)) +
    geom_histogram(aes(y = ..density..), alpha=0.7, fill="#33AADE", color="black") +
    geom_density(alpha=0.3, fill="red") +
    geom_vline(aes(xintercept=mean(pvalue)), color="black", linetype="dashed", size=1) +
    theme_bw() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(paste0(title, '\n',
                   'Total number of genes: ', nrow(res_df), '\n', 
                   sum(res_df$pvalue < 0.05),
                   ' genes have p-value < 0.05\n',
                   sum(res_df$pvalue >= 0.05),
                   ' genes have p-value >= 0.05\n',
                   sum(res_df$padj < 0.05),
                   ' genes have BH FDR < 0.05\n',
                   sum(res_df$padj >= 0.05),
                   ' genes have BH FDR >= 0.05'))
  
  # p <- ggplot(res_df, aes(x = pvalue)) +
  #   geom_histogram(binwidth = 0.05, center = 0.025) +
  #   theme_classic() +
  #   scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
  #   scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  #   xlab(xlab) +
  #   ylab(ylab) +
  #   ggtitle(paste0(title, '\n',
  #                  'Total number of genes: ', nrow(res_df), '\n', 
  #                  sum(res_df$pvalue < 0.05),
  #                  ' genes have p-value < 0.05\n',
  #                  sum(res_df$pvalue >= 0.05),
  #                  ' genes have p-value >= 0.05\n',
  #                  sum(res_df$padj < 0.05),
  #                  ' genes have BH FDR < 0.05\n',
  #                  sum(res_df$padj >= 0.05),
  #                  ' genes have BH FDR >= 0.05')) +
  #   theme(text = element_text(size=15))
  return(p)
}
