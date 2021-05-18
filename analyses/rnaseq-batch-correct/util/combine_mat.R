# function to merge input matrices (input must be collapsed)
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyverse))
combine_mat <- function(...){
  x <- lapply(..., FUN = function(x) readRDS(x))
  x <- lapply(x, FUN = function(x) x %>% rownames_to_column('gene'))
  x <- plyr::join_all(x, by = 'gene', type = 'inner')
  x <- x %>%
    column_to_rownames('gene')
  return(x)
}