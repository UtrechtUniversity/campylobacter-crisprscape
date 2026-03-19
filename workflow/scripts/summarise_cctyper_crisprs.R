#!/usr/bin/env Rscript

suppressPackageStartupMessages(
  library(tidyverse)
)

sink(
  file = file(snakemake@log[[1]], open = "wt"),
  type = "message"
)

cctyper_tables <- snakemake@input[["table"]]

read_cctyper_table <- function(cctyper_file) {
  batch <- snakemake@params[["batch"]]

  return(read_delim(cctyper_file, show_col_types = F) %>%
    mutate(batch = batch))
}


concatenated_cctyper_table <- do.call(
  rbind,
  lapply(X = cctyper_tables, FUN = read_cctyper_table)
)

write_delim(
  x = concatenated_cctyper_table, file = snakemake@output[[1]],
  delim = "\t"
)
