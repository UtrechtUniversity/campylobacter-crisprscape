#!/usr/bin/env Rscript

suppressPackageStartupMessages(
  library(tidyverse)
)

sink(
  file = file(snakemake@log[[1]], open = "wt"),
  type = "message"
)

tables <- snakemake@input[["table"]]

read_table <- function(file) {
  batch <- dirname(file) %>%
    basename()

  return(read_delim(file, show_col_types = F) %>%
    mutate(batch = batch))
}


concatenated_table <- do.call(
  rbind,
  lapply(X = tables, FUN = read_table)
)

write_delim(
  x = concatenated_table, file = snakemake@output[[1]],
  delim = "\t"
)
