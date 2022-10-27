# Create probe annotations for Illumina infinium methylation arrays using 
# Current GENCODE: GENCODE version 39 (Ensembl 105) gene symbols

# Eric Wafula for Pediatric OpenTargets
# 10/14/2022

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressWarnings(
  suppressPackageStartupMessages(library(rtracklayer))
)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# set up optparse options
option_list <- list(
  make_option(opt_str = "--probes_manifest", type = "character", default = NULL,
              help = "The latest Illumina Infinuim array probe manifest with 
              cpg annotation metadata.",  
              metavar = "character"),
  
  make_option(opt_str = "--gencode_gtf", type = "character", default = NULL,
              help = "The current GENCODE GTF utilized in OpenPedCan analyses
              modules.",
              metavar = "character")
)

# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
probes_manifest <- opt$probes_manifest
gencode_gtf <- opt$gencode_gtf

# establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to module and results directories
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "analyses", "methylation-summary")
results_dir <- file.path(module_dir, "results")

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Annotate array probes with GENCODE gene symbols and Ensembl IDs
# for the version currently utilized in the OpenPedCan analyses modules
message("===========================================================")
message("Creating annotations for methylation array probes")
message("===========================================================\n")

# using on-disk RSQLite database for to join tables; dplyr requires utilizes
# lots memory for table joins.
db_connect <- DBI::dbConnect(RSQLite::SQLite(), path = "")

# load array probe manifest 
DBI::dbWriteTable(db_connect, "manifest",
                  readr::read_csv(probes_manifest, 
                                  skip = 7, 
                                  guess_max = 10000) %>%
                    dplyr::select(Name, CHR, MAPINFO, 
                                  GencodeCompV12_Accession)  %>% 
                    dplyr::rename(Probe_ID = Name, 
                                  Chromosome = CHR, 
                                  Location = MAPINFO,
                                  transcript_id = GencodeCompV12_Accession) %>% 
                    tidyr::separate_rows(transcript_id, sep = ";", convert = FALSE) %>%
                    dplyr::mutate(transcript_id =
                                    stringr::str_extract(transcript_id, "ENST\\d+")) %>%
                    tidyr::drop_na() %>% 
                    dplyr::distinct()
)
                  
# load GENCODE GTF file
DBI::dbWriteTable(db_connect, "gencode",
                  rtracklayer::import(con = gencode_gtf) %>% 
                    tibble::as_tibble() %>%
                    dplyr::select(gene_id, transcript_id, gene_name) %>%
                    dplyr::rename(targetFromSourceId = gene_id, 
                                  Gene_symbol = gene_name) %>% 
                    dplyr::mutate(transcript_id = 
                                    stringr::str_extract(transcript_id, "ENST\\d+"), 
                                  targetFromSourceId = 
                                    stringr::str_extract(targetFromSourceId, "ENSG\\d+")
                                  )
)

# merge array probe annotations with GENCODE gene symbols and Ensembl IDs
probe_annotations <- dplyr::inner_join(dplyr::tbl(db_connect, "manifest"), 
                                       dplyr::tbl(db_connect, "gencode"), 
                                       by = "transcript_id") %>% 
  dplyr::distinct() %>% 
  tibble::as_tibble()

# remove database connection
dbDisconnect(db_connect)

# write array probe GENCODE annotations to file
message("Writing array probes GENCODE annotations to methyl-probe-annotations.tsv file...\n")
probe_annotations %>% 
  readr::write_tsv(file.path(results_dir, "methyl-probe-annotations.tsv.gz")) 

message("Analysis Done..\n")

