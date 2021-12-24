suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
})

process_annotate_overlaps <- function(cnv_df, exon_granges) {
  # This function takes a standardized data.frame that contains genomic range
  # information (cnv_df) and finds the overlaps with a GRanges object (exon_granges). 
  #
  # Args:
  #   cnv_df: standardized data.frame that contains the segments used in the CNV caller
  #   exon_granges: exons to be merged with the CNV data.frame;
  #
  # Returns:
  #   A data.frame with the following columns: 
  #   biospecimen_id, status, copy_number, ploidy, ensembl, gene_symbol, cytoband
  
  # make CNV data.frame into GRanges and sort by coordinates
  cnv_gr <- cnv_df %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                                            starts.in.df.are.0based = FALSE) %>%
    sort()
  
  # Create a data.frame with the overlaps between the CNV file and Gencode exon coordinates
  query <- cnv_gr
  subject <- exon_granges[exon_granges@seqnames %in% query@seqnames] # subset to only sequences in query
  overlapping_ranges <- findOverlaps(query = query, subject = subject)
  overlapping_ranges <- data.frame(query[queryHits(overlapping_ranges),], subject[subjectHits(overlapping_ranges),]) 
  
  # return object with specific columns
  overlapping_ranges <- overlapping_ranges %>%
    dplyr::rename("biospecimen_id" = "Kids_First_Biospecimen_ID",
                  "ploidy" = "tumor_ploidy",
                  "ensembl" = "gene_id",
                  "gene_symbol" = "gene_name") %>%
    dplyr::select(biospecimen_id, status, copy_number, ploidy, ensembl, gene_symbol, cytoband) %>%
    dplyr::distinct()
  
  return(overlapping_ranges)
}