# script to use HGG WT vs DMG H3K28M PBTA data as a use case for RUVSeq batch correction

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(uwot)
  library(tweeDEseq)
  library(DESeq2)
  library(tximport)
  library(tximportData)
})

# source functions
source('util/umap_plot.R')
source('util/edaseq_plot.R')
source('util/deseq2_pvals_histogram.R') # DESeq2 pval histograms
source('util/box_plots.R') # boxplots for samples
source('util/ruvg_test.R') # function to run RUVg
source('util/F_betaScoreFuns.R')
source('util/F_glm.nb2.R')
source('util/F_mgfNB.r')
source('util/F_odScore.R')
source('util/F_rowMultiply.R')
source('util/F_scoreNB.R')
source('util/F_smoothTestMat.R')
source('util/F_smoothTestPhylo.R')

# define directories
## data dir

dir <- system.file("extdata", package = "tximportData")
root_dir = rprojroot::find_root(rprojroot::has_dir('.git'))
data_dir = file.path(root_dir, 'data', 'v11')

## output dirs

umap_output_dir <- file.path('output', 'pbta_hgg_test', 'umap')
dge_output_dir <- file.path('output', 'pbta_hgg_test', 'dge')
plot_dir <- file.path('output', 'pbta_hgg_test', 'plots')

dir.create(dge_output_dir, showWarnings = F, recursive = T)
dir.create(umap_output_dir, showWarnings = F, recursive = T)
dir.create(plot_dir, showWarnings = F, recursive = T)

# histology file
hist_file <- read.delim(file.path(data_dir, 'histologies.tsv'))
hist.hgg <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         cohort %in% "PBTA", 
         harmonized_diagnosis %in% c('High-grade glioma/astrocytoma, H3 wildtype', 'Diffuse midline glioma, H3 K28-mutant'))

# HRT atlas genes
hk_genes_hrt <- readRDS('input/hk_genes_normals.rds')

# read expression data
# counts <- readRDS(file.path(data_dir, 'gene-counts-rsem-expected_count-collapsed.rds'))
# counts.hgg <- counts %>%
#   dplyr::select(hist.hgg$Kids_First_Biospecimen_ID)

# read transcript count data
file.metadata <- readr::read_tsv('~/Documents/temp_hgg_dmg_openpedcan/manifest_20220801_120351.tsv')
files <- file.metadata$name
files <- file.path('~/Documents/temp_hgg_dmg_openpedcan', files)
names(files) <- file.metadata$`Kids First Biospecimen ID`
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
rownames(txi.rsem[['counts']]) = gsub('_.+', '', rownames(txi.rsem[['counts']]))
rownames(txi.rsem[['length']]) = gsub('_.+', '', rownames(txi.rsem[['length']]))
rownames(txi.rsem[['abundance']]) = gsub('_.+', '', rownames(txi.rsem[['abundance']]))
txi.rsem <- tximport::summarizeToGene(txi.rsem, tx2gene = tx2gene)

# filter for zero count and protein coding genes from gencode v27
counts.hgg <- txi.rsem[['counts']]
counts.hgg = counts.hgg[rowSums(counts.hgg) > 0,] # remove all genes with zero counts across all samples
counts.hgg = as.data.frame(counts.hgg)

gencode_gtf <- rtracklayer::import(con = file.path(root_dir, 'data', 'gencode.v27.primary_assembly.annotation.gtf.gz'))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_pc <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

counts.hgg <- counts.hgg %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% gencode_pc$gene_id) %>%
  column_to_rownames("gene_id")

# diagnostic umap of HGG vs DMG
## high variance genes above 90% quantile (3666 x 19995)

counts.hgg.log = log2(counts.hgg +1)
row_variances <- matrixStats::rowVars(as.matrix(counts.hgg))
high_var_exp <- counts.hgg[which(row_variances > quantile(row_variances, 0.9)), ]
# high_var_exp.round <- na.omit(round(high_var_exp))
# high_var_exp.round <- as.matrix(high_var_exp.round)

fname <- file.path(umap_output_dir, 'umap_output_hgg_dmg.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(high_var_exp), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(high_var_exp)
  saveRDS(umap_out, file = fname)
}

p1 <- plot_data(umap_output = umap_out, 
                title = paste0("HGG WT vs DMG H3 K28"), 
                color_var = "harmonized_diagnosis", 
                shape_var = "RNA_library") + theme(legend.position = "right")

ggsave(p1, 
       filename = file.path(plot_dir, 'umap_output_hgg_dmg.pdf'), 
       width = 15, height = 6)

## housekeeping gene umap

hk_genes <- counts.hgg.log[rownames(counts.hgg.log) %in% hk_genes_hrt,]

fname <- file.path(umap_output_dir, 'umap_output_hgg_dmg_hk.rds')
if(!file.exists(fname)){
  set.seed(100)
  umap_out <- uwot::umap(X = t(hk_genes), n_components = 2, n_sgd_threads = 1, metric = "correlation")
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- colnames(hk_genes)
  saveRDS(umap_out, file = fname)
}

p2 <- plot_data(umap_output = umap_out, 
                title = paste0("HGG WT vs DMG H3 K28, Housekeeping genes only"), 
                color_var = "harmonized_diagnosis", 
                shape_var = "RNA_library") + theme(legend.position = "right")

ggsave(p2, 
       filename = file.path(plot_dir, 'umap_output_hgg_dmg_hk.pdf'), 
       width = 15, height = 6)

# test RUVseq vs no RUVseq DESeq2 models
#bs_id <- gsub('TARGET-[0-9]{2}-', '', selected_htl_df$Kids_First_Biospecimen_ID) # shorten bs ids for TARGET samples

## no RUVseq DGE
### build design matrix

harmonized_diagnosis <- factor(as.character(hist.hgg$harmonized_diagnosis))
design <- model.matrix(~ 0 + harmonized_diagnosis)
bs_id <- hist.hgg$Kids_First_Biospecimen_ID
RNA_library = hist.hgg$RNA_library
dds <- DESeqDataSetFromTximport(txi.rsem, hist.hgg, design = design)

### filter lowly expressed genes
counts_object <- edgeR::DGEList(counts = dds@assays@data@listData$counts, group = harmonized_diagnosis)
counts_object <- calcNormFactors(counts_object)
keep <- edgeR::filterByExpr(counts_object, group = harmonized_diagnosis)
counts_object_filtered <- counts_object[keep, , keep.lib.sizes = F]
counts_object_filtered <- estimateDisp(counts_object_filtered, design)
fit <- glmQLFit(counts_object_filtered, design)
qlf <- glmQLFTest(fit, contrast=c(1,0))
res <- topTags(qlf, n = Inf)

### build EDA expressionset
seq_expr_set <- EDASeq::newSeqExpressionSet(counts = dds@assays@data@listData$counts, phenoData = data.frame(bs_id, RNA_library, harmonized_diagnosis, row.names = colnames(dds@assays@data@listData$counts)), featureData = as.data.frame(dds@assays@data@listData$avgTxLength))

# from https://support.bioconductor.org/p/95805/#95808
# The "betweenLaneNormalization" is just a poorly named function that perform between-sample normalization, 
# independently of the fact that they are organized in lanes, plates, etc. 
# It is used in RUVSeq just to make sure that the factors of unwanted variation don't capture the difference in library size that should be captured by the size factors. 
# normalize the data using upper-quartile (UQ) normalization
# seq_expr_set <- EDASeq::betweenLaneNormalization(seq_expr_set, which = "upper")
umap_p_norm <- edaseq_plot(object = seq_expr_set, title = "After UQ normalization", type = "UMAP", color_var = "harmonized_diagnosis", shape_var = "RNA_library")

p <- ggpubr::ggarrange(umap_p_norm, common.legend = T, legend = "bottom")
fname <- file.path(umap_output_dir, 'umap_hgg_uq_norm.pdf')
ggsave(filename = fname, plot = p, width = 6, height = 6, device = "pdf", bg = "white")

# 1.2. RUVg using hk genes in normals only (full HRT atlas)
emp_neg_ctrl_genes_normals <- readRDS('input/hk_genes_normals.rds')
k_value = 5
ruvg_test(seq_expr_set = seq_expr_set, k_val = 1:k_value, 
          emp_neg_ctrl_genes = emp_neg_ctrl_genes_normals, prefix = "hk_genes_normals", 
          diff_type = "deseq2", 
          output_dir = dge_output_dir,
          design_variable = "harmonized_diagnosis",
          color_var = "harmonized_diagnosis", shape_var = "RNA_library")

# 1. DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts(seq_expr_set)), colData =  data.frame(bs_id, harmonized_diagnosis), design = design)
dds <- DESeq2::DESeq(dds)
dge_output <- DESeq2::results(dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
dge_output <- dge_output %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>%
  arrange(padj) 

ind = which(dge_output$padj < 0.05)
dge_output.f = dge_output[ind,]
dge_output.f = dge_output.f[which(dge_output.f$gene %in% gencode_pc$gene_id),]
dge_output.f = merge(dge_output.f, gencode_pc, by.x = 'gene', by.y = 'gene_id')


# plot and save p-value histogram
p <- deseq2_pvals_histogram(res_df = dge_output,
                            xlab = 'DGE RLE nbinomWaldTest p-value',
                            ylab = 'Gene count', title = paste0('Histogram of DESeq2 nbinomWaldTest p-values'))
filename <- file.path(output_dir, 'deseq2_analysis', 'dge_deseq2_histogram.pdf')
ggsave(filename = filename, plot = p, width = 8, height = 7, bg = "white")


test = tweeDEseq::gofTest(high_var_exp.round, a = 0)
test.sort = sort(pchisq(test, df=1, lower.tail=FALSE))
res.test = data.table('test.res' = test.sort, 'adj.p' = p.adjust(test.sort, method = 'BH'))
length(which(res.test$adj.p < 0.05))

qq <- qqchisq(test, main="Chi2 Q-Q Plot", ylim = c(0, 60))

mat = t(high_var_exp)
test2 = smoothTestMat(mat = t(high_var_exp)[1:10,], x = design)
