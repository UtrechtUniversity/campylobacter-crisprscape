#!/usr/bin/env python3

# Convert spacer table data to visualisation formats for
# CCTK crisprdiff <https://crispr-comparison-toolkit.readthedocs.io/en/latest/content/minced.html#array-ids-txt>
# and SpacerPlacer <https://github.com/fbaumdicker/SpacerPlacer#spacer_fasta-input-format>

import pandas as pd

spacer_table = snakemake.input[0]
table_numbers = snakemake.output["table_numbers"]
table_sequences = snakemake.output["table_sequences"]
fasta_numbers = snakemake.output["fasta_numbers"]
fasta_sequences = snakemake.output["fasta_sequences"]

spacer_df = pd.read_csv(spacer_table, sep="\t")

# Convert cluster ID numbers to string to enable 'joining' them
spacer_df["Cluster"] = spacer_df["Cluster"].apply(str)

# Create the tab-separated table with cluster ID numbers per array
numbers_df = (
    spacer_df[["Spacer_base_ID", "Cluster"]]
    .groupby(["Spacer_base_ID"])
    .agg({"Cluster": " ".join})
)
numbers_df.to_csv(table_numbers, sep="\t", header=False)

# Create the tab-separated table with spacer sequences per array
sequence_df = (
    spacer_df[["Spacer_base_ID", "Sequence"]]
    .groupby(["Spacer_base_ID"])
    .agg({"Sequence": " ".join})
)
sequence_df.to_csv(table_sequences, sep="\t", header=False)

# Create fasta file with cluster ID numbers per array
numbers_fasta_df = (
    spacer_df[["Spacer_base_ID", "Cluster"]]
    .groupby(["Spacer_base_ID"])
    .agg({"Cluster": ", ".join})
)
with open(fasta_numbers, "w") as outfile:
    for name, numbers in numbers_fasta_df.agg("".join, axis=1).items():
        outfile.write(f">{name}\n{numbers}\n")

# Create fasta file with spacer sequences per array
numbers_fasta_df = (
    spacer_df[["Spacer_base_ID", "Sequence"]]
    .groupby(["Spacer_base_ID"])
    .agg({"Sequence": ", ".join})
)
with open(fasta_sequences, "w") as outfile:
    for name, sequences in numbers_fasta_df.agg("".join, axis=1).items():
        outfile.write(f">{name}\n{sequences}\n")

exit(0)
