# function to rbind all input metadata
suppressPackageStartupMessages(library(dplyr))
combine_meta <- function(..., cols){
  x <- lapply(..., FUN = function(x) read.delim(x, stringsAsFactors = F, check.names = F))
  x <- lapply(x, FUN = function(x) x %>% dplyr::select(cols))
  x <- do.call("rbind", x)
  return(x)
}
