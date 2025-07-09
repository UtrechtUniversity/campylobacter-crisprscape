#!/usr/bin/env Rscript

# Find output files for geNomad and create one overall prediction report.

library(tidyverse)
library(here)

genomad_scores_files <- snakemake@input[["aggregated_classification"]]
genomad_plasmid_files <- snakemake@input[["plasmid_summary"]]
genomad_virus_files <- snakemake@input[["virus_summary"]]

read_stats <- function(filename, name_position) {
  # Cut the sample name from the file path using the name's position in the
  # absolute file path (counting from the end: 1 = last, 2 = second last, etc.)
  sample_name <- str_split_1(string = filename, pattern = "/") %>%
    tail(name_position) %>%
    head(1)
  
  df <- read_delim(
    file = filename,
    show_col_types = FALSE
  ) %>%
    mutate(batch = sample_name)
  
  return(df)
}

genomad_scores <- do.call(
  rbind,
  lapply(X = genomad_scores_files, FUN = read_stats, name_position = 3)
) %>%
  rename("contig" = "seq_name")

genomad_scores <- genomad_scores %>%
  mutate(genome = gsub(pattern = "\\.contig[0-9]*",
                       replacement = "",
                       x = contig)
  )

plasmid_classifications <- do.call(
  rbind,
  lapply(X = genomad_plasmid_files, FUN = read_stats, name_position = 3)
) %>%
  rename(
    "contig" = "seq_name",
    "plasmid_topology" = "topology",
    "plasmid_genes" = "n_genes"
  )

virus_classifications <- do.call(
  rbind,
  lapply(X = genomad_virus_files, FUN = read_stats, name_position = 3)
) %>%
  rename(
    "contig" = "seq_name",
    "virus_topology" = "topology",
    "virus_genes" = "n_genes",
    "virus_taxonomy" = "taxonomy"
  )

genomad_df <- left_join(
  x = genomad_scores %>%
    select(genome, batch, contig, chromosome_score, plasmid_score, virus_score),
  y = plasmid_classifications %>%
    select(contig, plasmid_topology, plasmid_genes, conjugation_genes, amr_genes),
  by = "contig"
) %>%
  left_join(
    x = .,
    y = virus_classifications %>%
      select(contig, virus_topology, virus_genes, virus_taxonomy),
    by = "contig"
  ) %>%
  mutate(
    genomad_prediction = case_when(
      !is.na(plasmid_topology) ~ "plasmid",
      !is.na(virus_topology) ~ "virus",
      TRUE ~ "chromosome"
      )
  )

write_csv(
  x = genomad_df,
  file = snakemake@output[[1]] #here("data", "processed", "genomad_predictions.csv")
)