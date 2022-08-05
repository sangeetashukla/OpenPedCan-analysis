# script to use HGG WT vs DMG H3K28M PBTA data as a use case for RUVSeq batch correction

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(uwot)
  library(DESeq2)
  library(tximport)
  library(tximportData)
  library(MASS)
  library(stats)
})

# source functions
source('util/umap_plot.R')
source('util/edaseq_plot.R')
source('util/deseq2_pvals_histogram.R') # DESeq2 pval histograms
source('util/box_plots.R') # boxplots for samples
source('util/ruvg_test.R') # function to run RUVg

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

# gencode reference

gencode_gtf <- rtracklayer::import(con = file.path(root_dir, 'data', 'gencode.v27.primary_assembly.annotation.gtf.gz'))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_pc <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()
# read transcript count data

file.metadata <- readr::read_tsv('~/Documents/temp_hgg_dmg_openpedcan/manifest_20220801_120351.tsv')
files <- file.metadata$name
files <- file.path('~/Documents/temp_hgg_dmg_openpedcan', files)
names(files) <- file.metadata$`Kids First Biospecimen ID`

# obtain transcript to gene mappings
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)

# remove _geneName and retain only transcript IDs for tximport to map transcripts to genes
rownames(txi.rsem[['counts']]) = gsub('_.+', '', rownames(txi.rsem[['counts']]))
rownames(txi.rsem[['length']]) = gsub('_.+', '', rownames(txi.rsem[['length']]))
rownames(txi.rsem[['abundance']]) = gsub('_.+', '', rownames(txi.rsem[['abundance']]))
txi.rsem <- tximport::summarizeToGene(txi.rsem, tx2gene = tx2gene)

# build DESeq2 dataset from tximport and specify design. 

harmonized_diagnosis <- factor(as.character(hist.hgg$harmonized_diagnosis))
design <- model.matrix(~harmonized_diagnosis)
bs_id <- hist.hgg$Kids_First_Biospecimen_ID
RNA_library = hist.hgg$RNA_library
dds <- DESeqDataSetFromTximport(txi.rsem, hist.hgg, design = design)

# build expressionset for EDA analysis and downstream DGE analysis
seq_expr_set <- EDASeq::newSeqExpressionSet(counts = dds@assays@data@listData$counts, phenoData = data.frame(bs_id, RNA_library, harmonized_diagnosis, row.names = colnames(dds@assays@data@listData$counts)), featureData = as.data.frame(dds@assays@data@listData$avgTxLength))

# 1. DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
# dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts(seq_expr_set), colData =  data.frame(bs_id, harmonized_diagnosis), design = design)

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

# test for goodness of fit
counts.hgg = txi.rsem$counts
counts.hgg = counts.hgg[which(rownames(counts.hgg) %in% gencode_pc$gene_id),]

