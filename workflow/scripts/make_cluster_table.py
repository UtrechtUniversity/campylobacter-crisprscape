#!/usr/bin/env python3

import re
from pyfaidx import Fasta

cluster_file = snakemake.input["clstr"]
fasta_file = snakemake.input["fasta"]
table_file = snakemake.output[0]


def read_clusters(clstr, table, fasta):
    """
    Read a CD-HIT generated file of clusters and the fasta file that was
    used by CD-HIT to generate a tab-separated report of the clusters.

        clstr: the CD-HIT output file with .clstr extensin
        table: output file to write tab-separated output file to
        fasta: fasta file with input sequences
    """
    sequence_dict = find_sequences(fasta)

    HEADER = "Genome\tContig\tLocus\tCluster\tLength\tCluster_representative\tSequence\tIdentity\tStrand\n"
    cluster_regex = r"(>Cluster *)(\d*)"
    locus_regex = r"^(\d+)\s+(\d+nt), >(\w+).(contig[\d-]+_\d+:\d+)... (.*)$"

    with open(table, "w") as outfile:
        outfile.write(HEADER)

        with open(clstr, "r") as infile:
            for line in infile:
                line = line.strip()

                if line.startswith(">"):
                    # Extract the digits from the cluster ID
                    cluster = re.search(cluster_regex, line).group(2)

                elif len(line) > 1:
                    # Use RegEx to extract information
                    crispr_info = re.search(locus_regex, line)

                    member_nr = crispr_info.group(1)  # not used
                    length = crispr_info.group(2)
                    genome = crispr_info.group(3)
                    locus = crispr_info.group(4)
                    full_locus = "%s.%s" % (genome, locus)
                    contig = locus.split("_")[0]
                    extra = crispr_info.group(5)

                    sequence = sequence_dict[full_locus]

                    # Check the final group for representative ('*') or other
                    if extra == "*":
                        representative = full_locus
                        strand = "NA"
                        identity = "NA"
                    else:
                        strand_and_identity = extra.split("/")
                        strand = strand_and_identity[0].replace("at ", "")
                        identity = strand_and_identity[1]

                    # Write the information to the output file
                    crispr_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                        genome,
                        contig,
                        full_locus,
                        cluster,
                        length,
                        representative,
                        sequence,
                        identity,
                        strand,
                    )
                    outfile.write(crispr_line)

                # If the line does not start with '>' or have length > 1, stop
                else:
                    break

    return 0


def find_sequences(fasta_file):
    """
    Look up the DNA sequence in a fasta file and return as dictionary.
    """
    sequence_dict = Fasta(fasta_file, duplicate_action="first")
    return sequence_dict


if __name__ == "__main__":
    exit(read_clusters(clstr=cluster_file, table=table_file, fasta=fasta_file))
