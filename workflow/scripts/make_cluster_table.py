#!/usr/bin/env python3


################################################################
# Read CD-HIT's output cluster file to generate a table of     #
# clustered sequences.                                         #
# - The program:                                               #
#   1) reads the cluster file to                               #
#      a) create a sequence-based table of:                    #
#         - Seq. length                                        #
#         - Cluster ID                                         #
#         - Sequence                                           #
#         - # Occurrences                                      #
#                                                              #
#      b) calculate per cluster:                               #
#         - the number of members (spacers in the cluster)     #
#         - the shortest, longest and most common sequence     #
#                                                              #
#      c) match cluster ID back with the spacer table          #
#         - that is, 'results/spacers-final_prep.tsv'          #
#                                                              #
#   2) Save (b) as cluster table, and (c) as spacer table      #
#         - both as tab-separated values (TSV)                 #
################################################################

import re
import pandas as pd

cluster_file = snakemake.input["clstr"]
cluster_table = snakemake.output["cluster"]
all_spacer_table = snakemake.input["spacer"]
spacer_table = snakemake.output["spacer"]


def read_clstr_file(inputfile=str):
    """
    Read a CD-HIT generated .clstr file and parse its elements into
     a dictionary of: 1) spacer sequence, 2) sequence length,
                      3) number of occurrences, 4) cluster ID
    """
    spacer_dict = {"Length": [], "Cluster": [], "Sequence": [], "Occurrences": []}

    cluster_regex = r"(>Cluster *)(\d*)"
    sequence_regex = r"^(\d+)\s+(\d+nt), >(\d+)\.(\w+)\.\.\. (.*)$"

    with open(inputfile, "r") as infile:
        for line in infile:
            line = line.strip()  # Remove funny characters

            if line.startswith(">"):
                # Extract the digits from the cluster ID
                cluster = int(re.search(cluster_regex, line).group(2))

            elif len(line) > 1:
                # Use RegEx to extract information
                crispr_info = re.search(sequence_regex, line)

                member_nr = crispr_info.group(1)  # not used
                length = int(crispr_info.group(2).rstrip("nt"))
                prevalence = int(crispr_info.group(3))
                sequence = crispr_info.group(4)
                extra = crispr_info.group(5)  # not used

                spacer_dict["Length"].append(length)
                spacer_dict["Cluster"].append(cluster)
                spacer_dict["Sequence"].append(sequence)
                spacer_dict["Occurrences"].append(prevalence)

            # If the line does not start with '>' or have length > 1, stop
            else:
                break

        return spacer_dict


def make_cluster_overview(df=pd.DataFrame):
    """
    Given a pandas dataframe with spacer clusters, summarise the most common,
    shortest, and longest sequence within each cluster. Also calculate the
    number of spacers present in the clusters.
    """
    longest_per_cluster = (
        df[["Cluster", "Sequence", "Length"]].groupby(["Cluster"]).max()
    )
    shortest_per_cluster = (
        df[["Cluster", "Sequence", "Length"]].groupby("Cluster").min()
    )
    most_common_per_cluster = (
        df[["Cluster", "Sequence", "Length", "Occurrences"]]
        .groupby("Cluster", sort=False)
        .first()
    )
    spacers_per_cluster = df[["Cluster", "Occurrences"]].groupby("Cluster").sum()

    # Rename columns to facilitate joining dataframes
    shortest_per_cluster = shortest_per_cluster.rename(
        columns={"Sequence": "Shortest_sequence", "Length": "Shortest_length"}
    )
    longest_per_cluster = longest_per_cluster.rename(
        columns={"Sequence": "Longest_sequence", "Length": "Longest_length"}
    )
    most_common_per_cluster = most_common_per_cluster.rename(
        columns={
            "Sequence": "Most_common_sequence",
            "Occurrences": "Most_common_number",
            "Length": "Most_common_length",
        }
    )
    spacers_per_cluster = spacers_per_cluster.rename(
        columns={"Occurrences": "Number_of_spacers"}
    )

    # Join cluster IDs with most common, shortest and longest spacer sequence
    spacer_cluster_table = (
        spacers_per_cluster.join(most_common_per_cluster, on="Cluster")
        .join(shortest_per_cluster, on="Cluster")
        .join(longest_per_cluster, on="Cluster")
    )

    return spacer_cluster_table


def update_spacer_table(table=str, df=pd.DataFrame):
    """
    Given the table with spacer IDs, sequences and length, update it to include
    the respective cluster ID and also separate the ID into a 'base' and
    'position' field for parsing the arrays.
    """
    all_spacer_df = pd.read_csv(table, sep="\t")
    all_spacer_df = all_spacer_df.rename(
        columns={"#name": "Spacer_ID", "seq": "Sequence", "length": "Length"}
    )
    all_spacer_df = all_spacer_df.merge(df, on="Sequence")

    # Also separate the spacer ID from the number, to get a 'clean' base ID
    #  and number in the array (1, 2, ..., n)
    # Note: This is different between CCTyper, which uses format '[array]:[spacer#]'
    # and CRISPRidentify, which uses '[array]_spacer_[spacer#]'
    if all_spacer_df["Spacer_ID"].str.contains(":").all():
        # If all values contain a colon, it's the CCTyper format
        array_and_number = all_spacer_df["Spacer_ID"].str.rsplit(":", expand=True)
        all_spacer_df["Spacer_base_ID"] = array_and_number[0]
        all_spacer_df["Spacer_position"] = array_and_number[1]
    else:
        # If there is no colon, it's the CRISPRidentify format
        array_and_number = all_spacer_df["Spacer_ID"].str.rsplit("_", expand=True, n=2)
        all_spacer_df["Spacer_base_ID"] = array_and_number[0]
        all_spacer_df["Spacer_position"] = array_and_number[2]

    return all_spacer_df


def main():
    """
    Main function, running the whole script.
    """
    # Read clstr file, store as dictionary
    spacer_dict = read_clstr_file(inputfile=cluster_file)

    # Convert dictionary to dataframe
    spacer_df = pd.DataFrame(spacer_dict)

    ### 1. Make an overview of clusters ###
    # Find most common, shortest, and longest sequence per cluster,
    # and count total number of spacers per cluster.
    spacer_cluster_table = make_cluster_overview(df=spacer_df)

    # Save as tab-separated text file
    spacer_cluster_table.to_csv(cluster_table, sep="\t", index=True)

    ### 2. Merge cluster IDs into the existing spacer table
    # Take the table of all spacers from all genomes and add in cluster IDs
    all_spacer_df = update_spacer_table(
        table=all_spacer_table, df=spacer_df[["Sequence", "Cluster"]]
    )

    # Save as tab-separated text file
    all_spacer_df.to_csv(spacer_table, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    exit(main())
