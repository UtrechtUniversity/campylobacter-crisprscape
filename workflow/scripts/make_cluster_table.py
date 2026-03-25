#!/usr/bin/env python3

################################################################
# Read CD-HIT's output cluster file, along with the FASTA used #
# to generate it to create a table of clustered sequences.     #
# - Aim of the program:                                        #
#   1) read the cluster file                                   #
#      a) create a cluster-based dictionary (cluster #,        #
#         representative)                                      #
#      b) create a sequence-based dictionary (ID, length,      #
#         representative or % identity and strand)             #
#   2) combine the two dictionaries into a dataframe with      #
#       columns: genome, contig, locus, cluster ID, length,    #
#       longest_sequence (CD-HIT representative), sequence,    #
#       identity (%), strand                                   #
#   3) for each cluster, also mark the shortest and most common#
#       sequence (as alternative representatives)              #
#   4) write the dataframe to a tab-separated table file       #
################################################################

import re
from pyfaidx import Fasta
import pandas as pd

cluster_file = snakemake.input["clstr"]
fasta_file = snakemake.input["fasta"]
table_file = snakemake.output[0]


def read_clstr_file(inputfile=str):
    """
    Read a CD-HIT generated .clstr file and parse its elements into
     two dictionaries: 1) cluster ID + representative (locus ID),
                       2) sequences (locus ID) with length and
                        representative/identity
    """
    cluster_dict = {"Cluster": [], "Longest_sequence": []}
    locus_dict = {
        "Genome": [],
        "Contig": [],
        "Locus": [],
        "Full_locus": [],
        "Length": [],
        "Cluster": [],
        "Strand_to_longest": [],
        "Identity_to_longest": [],
    }

    cluster_regex = r"(>Cluster *)(\d*)"
    locus_regex = r"^(\d+)\s+(\d+nt), >(\w+)([\.-])(contig[\d-]+_.+)\.\.\. (.*)$"

    with open(inputfile, "r") as infile:
        for line in infile:
            line = line.strip()  # Remove funny characters

            if line.startswith(">"):
                # Extract the digits from the cluster ID
                cluster = re.search(cluster_regex, line).group(2)

            elif len(line) > 1:
                # Use RegEx to extract information
                crispr_info = re.search(locus_regex, line)

                member_nr = crispr_info.group(1)  # not used
                length = int(crispr_info.group(2).rstrip("nt"))
                genome = crispr_info.group(3)
                separator = crispr_info.group(4)
                locus = crispr_info.group(5)
                full_locus = f"{genome}{separator}{locus}"
                contig = locus.split("_")[0]
                extra = crispr_info.group(6)

                # Check the final group for representative ('*') or other
                if extra == "*":
                    strand = "NA"
                    identity = "NA"
                    cluster_dict["Cluster"].append(cluster)
                    cluster_dict["Longest_sequence"].append(full_locus)

                else:
                    strand_and_identity = extra.split("/")
                    strand = strand_and_identity[0].replace("at ", "")
                    identity = strand_and_identity[1]

                locus_dict["Genome"].append(genome)
                locus_dict["Contig"].append(contig)
                locus_dict["Locus"].append(locus)
                locus_dict["Full_locus"].append(full_locus)
                locus_dict["Length"].append(length)
                locus_dict["Cluster"].append(cluster)
                locus_dict["Strand_to_longest"].append(strand)
                locus_dict["Identity_to_longest"].append(identity)

            # If the line does not start with '>' or have length > 1, stop
            else:
                break

        return cluster_dict, locus_dict


def generate_sequence_df(fasta=str, ids=list):
    """
    Given a fasta file and list of identifiers, create a dataframe
    of DNA sequences that can be merged with the cluster/locus dataframe.
    """
    sequence_dict = Fasta(fasta, duplicate_action="first")
    sequence_list = []
    for locus in ids:
        sequence_list.append(sequence_dict[locus][:].seq)

    return pd.DataFrame({"Full_locus": ids, "Sequence": sequence_list})


def main():
    """
    Main function, running the whole script.
    """
    # Read clstr file, store as dictionaries
    cluster_dict, locus_dict = read_clstr_file(inputfile=cluster_file)

    # Convert dictionaries to dataframes
    cluster_df = pd.DataFrame(cluster_dict)
    locus_df = pd.DataFrame(locus_dict)

    # Merge dataframes
    combined_df = locus_df.merge(cluster_df, how="inner", on="Cluster")

    # Add sequences
    sequence_df = generate_sequence_df(fasta=fasta_file, ids=locus_df["Full_locus"])

    combined_with_sequences = combined_df.merge(
        sequence_df, how="inner", on="Full_locus"
    )

    ## Not yet implemented:
    # Find shortest sequence per cluster
    # Find most common sequence per cluster

    # Save as tab-separated text file
    combined_with_sequences.to_csv(table_file, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    exit(main())
