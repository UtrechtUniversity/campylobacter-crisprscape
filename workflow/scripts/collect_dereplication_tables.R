#!/usr/bin/env Rscript

sink(
  file = file(snakemake@log[[1]], open = "wt"), type = "message"
)

suppressPackageStartupMessages(
  library(tidyverse)
)

cluster_tables <- snakemake@input[["cluster"]]
winner_tables <- snakemake@input[["winner"]]

read_clusters <- function(clusters) {
  batch <- str_split_i(
    string = clusters,
    pattern = "/",
    i = 3
  )
  cluster_table <- read_csv(
    file = clusters,
    show_col_types = FALSE
  )

  cluster_table <- cluster_table %>%
    mutate(batch = batch)

  return(cluster_table)
}

read_winners <- function(winners) {
  batch <- str_split_i(
    string = winners,
    pattern = "/",
    i = 3
  )

  winner_table <- read_csv(
    file = winners,
    show_col_types = FALSE
  ) %>%
    rename("representative" = "genome")

  winner_table <- winner_table %>%
    mutate(batch = batch)

  return(winner_table)
}

cluster_table <- do.call(
  rbind,
  lapply(X = cluster_tables, FUN = read_clusters)
)

winner_table <- do.call(
  rbind,
  lapply(X = winner_tables, FUN = read_winners)
)

dereplication_df <- left_join(
  x = cluster_table,
  y = winner_table,
  by = c("batch",
    "secondary_cluster" = "cluster"
  )
)

write_delim(
  x = dereplication_df %>%
    # Manually set column order:
    select(
      batch, genome, secondary_cluster, representative, score, cluster_method,
      comparison_algorithm, threshold, primary_cluster
    ),
  file = snakemake@output[[1]],
  delim = "\t"
)
