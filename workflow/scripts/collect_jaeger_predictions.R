#!/usr/bin/env Rscript

sink(
  file = file(snakemake@log[[1]], open = "wt"), type = "message"
)

suppressPackageStartupMessages(
  library(tidyverse)
)

# Read all output files in a batch (one per genome)
jaeger_files <- snakemake@input[["default"]]
phage_files <- snakemake@input[["phages"]]

# Concatenate them all in one dataframe
jaeger_df <- do.call(
  rbind,
  lapply(X = jaeger_files, FUN = read_delim, show_col_types = FALSE)
)

phage_df <- do.call(
  rbind,
  lapply(X = phage_files, FUN = read_delim, show_col_types = FALSE)
)

# Simplify the dataframe by extracting essential info
simplify_df <- function(df) {
  return(
    df %>%
      mutate(
        contig = gsub(
          pattern = " .*",
          replacement = "",
          x = contig_id
        ),
        accession_id = gsub(
          pattern = ".contig[0-9]*",
          replacement = "",
          x = contig
        )
      ) %>%
      select(accession_id, contig, length, prediction, reliability_score, prophage_contam)
  )
}

jaeger_df_simple <- simplify_df(jaeger_df)
phage_df_simple <- simplify_df(phage_df)

# And write to a CSV file (which can easily be concatenated with a script)
write_delim(
  x = jaeger_df_simple,
  file = snakemake@output[["default"]],
  delim = "\t"
)

write_delim(
  x = phage_df_simple,
  file = snakemake@output[["phages"]],
  delim = "\t"
)
