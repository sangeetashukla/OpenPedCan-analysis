setwd('~/Projects/PediatricOpenTargets/OpenPedCan-analysis/analyses/rnaseq-batch-correct/')

suppressPackageStartupMessages({
  library(tidyverse)
})

# histology
hist_file <- read.delim('../../data/histologies.tsv')

# subset to tumor samples
hist_file <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         cohort %in% c("GMKF", "PBTA", "TARGET"),
         sample_type == "Tumor")

# counts
expr_counts <- readRDS('../../data/gene-counts-rsem-expected_count-collapsed.rds')
expr_counts <- expr_counts %>%
  dplyr::select(hist_file$Kids_First_Biospecimen_ID)

# compute gene lengths
gencode_gtf <- rtracklayer::import(con = file.path('../../data/gencode.v27.primary_assembly.annotation.gtf.gz'))
gencode_gtf <- as.data.frame(gencode_gtf)
gene_lengths <- rowsum(gencode_gtf$end-gencode_gtf$start+1, gencode_gtf$gene_name)
gene_lengths <- gene_lengths[rownames(expr_counts),]

# create DGEList object
y <- edgeR::DGEList(counts = expr_counts, genes = data.frame(Length = gene_lengths))
y <- edgeR::calcNormFactors(y, method = "TMM")
rpkm_mat <- rpkm(y)

# apply filters
# 1. the gene must be expressed at non-zero level in all tissue and cell types
rpkm_mat <- rpkm_mat %>%
  as.data.frame()
rpkm_mat <- rpkm_mat[rowSums(rpkm_mat) > 0,]

# 2. the variability of transcript expression should be low within all tissues and cell types, 
# as evidence by a standard deviation of the log2 RPKM < 1;
row_sd <- apply(log2(rpkm_mat + 1), MARGIN = 1, sd)
rpkm_mat <- rpkm_mat[row_sd < 1,]

# 3. the maximum fold change (MFC), represented by the 
# ratio between maximum and average log2 RPKM of the transcript ((maximum log2 RPKM)/(average log2 RPKM)), must be lower than 2.
row_mfc <- apply(log2(rpkm_mat + 1), MARGIN = 1, FUN = function(x) max(x)/mean(x))
rpkm_mat <- rpkm_mat[row_mfc < 2,]

# total tumor-specific HK genes
dim(rpkm_mat)
# 37 2601

# overlap with HRT atlas
load('input/Housekeeping_GenesHuman.RData')
hk_tumors_normals <- intersect(rownames(rpkm_mat), Housekeeping_Genes$Gene.name) # 25
saveRDS(hk_tumors_normals, file = 'output/hk_genes_tumor_normals.rds')

