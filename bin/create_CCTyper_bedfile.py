#!/usr/bin/env python3

# Read the tabular output file 'CRISPR_Cas.tab' from CCTyper
# and convert it to a BED file. This way, seqkit can be easily
# used to extract the CRISPR-Cas array and its flanking regions.

import argparse

def parse_arguments():
    """
    Parse argument from the command line:
     -i/--input = input file (CRISPR_Cas.tab)
     -o/--output = output file (BED file)
     -h/--help = show help
    """
    parser = argparse.ArgumentParser(
                        prog="CCTyper_to_bedfile",
                        description="Convert the CRISPR_Cas.tab output file"
                        "from CCTyper to a BED file.",
    )

    required = parser.add_argument_group("Required arguments")

    required.add_argument("-i", "--input",
    dest = "input", required = True, type = str,
    help = "CCTyper CRISPR_Cas.tab file",
    )
    required.add_argument("-o", "--output",
    dest = "output", required = True, type = str,
    help = "BED file to write output to")

    args = parser.parse_args()

    return args

def convert_tab_to_bed(inputfile, outputfile):
    """
    """
    with open(outputfile, 'w') as outfile:
        with open(inputfile, 'r') as infile:
            # Skip the first line, as it is the header
            infile.readline()
            for line in infile:
                # Split by tab
                elements = line.split("\t")

                # Collect required information
                contig = elements[0]
                locus_tag = elements[1]
                start_and_stop = elements[2].split(",")
                start = start_and_stop[0].strip(" [")
                stop = start_and_stop[1].strip(" ]")

                # Combine in single line
                bed_string = "%s\t%s\t%s\t%s\n" % (contig, start, stop, locus_tag)

                # Write to the output file
                outfile.write(bed_string)
    return 0

def main():
    arguments = parse_arguments()

    message = (
        "\n"
        "These are the files you have provided:\n"
        "  INPUT:\n"
        "{0}\n"
        "  OUTPUT:\n"
        "{1}".format(
            arguments.input,
            arguments.output
        )
    )

    print(message)

    exit(
        convert_tab_to_bed(inputfile = arguments.input,
                           outputfile = arguments.output)
    )

if __name__=="__main__":
    exit(main())