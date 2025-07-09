#!/usr/bin/env Rscript

# Concatenate Jaeger output files per batch while selecting only the most
# relevant, easy-to-use columns.

library(here)
library(tidyverse)

# Read all output files in a batch (one per genome)
jaeger_files <- Sys.glob(paths = here("data", "tmp", "jaeger", snakemake@params[["batch"]], "*", "*_default_jaeger.tsv"))

# Concatenate them all in one dataframe
jaeger_df <- do.call(
  rbind,
  lapply(X = jaeger_files, FUN = read_delim, show_col_types = FALSE)
)

# Simplify the dataframe by extracting essential info
jaeger_df_simple <- jaeger_df %>%
  mutate(
    contig = gsub(pattern = " .*",
                  replacement = "",
                  x = contig_id),
    accession_id = gsub(pattern = ".contig[0-9]*",
                        replacement = "",
                        x = contig)
  ) %>%
  rename("reliability_score" = "realiability_score") %>%
  select(accession_id, contig, length, prediction, reliability_score, prophage_contam)

# And write to a CSV file (which can easily be concatenated with a script)
write_csv(x = jaeger_df_simple,
          file = snakemake@output[[1]])