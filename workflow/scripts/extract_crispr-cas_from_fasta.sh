#!/usr/bin/env bash

exec > "${snakemake_log[0]}" 2>&1 # send stdout+strerr to a log file

fasta_file="${snakemake_input[genomes]}"
bed_file="${snakemake_input[bed]}"
out_dir="${snakemake_output[bed_dir]}"

echo "Extracting CRISPR-Cas from ${fasta_file}..."

seqkit subseq --bed "${bed_file}" -U "${fasta_file}" |\
# Remove the original ID with positions, up to the space,
#  leaving only the CRISPR ID. Also replace dot by dash.
 seqkit replace -p ".+\s" | seqkit replace -p "\." -r "-"\
 -o "${snakemake_output[crispr_cas]}"

seqkit subseq --bed "${bed_file}" -U "${fasta_file}" -u 5000 -d 5000 |\
 seqkit replace -p ".+\s" | seqkit replace -p "\." -r "-"\
 -o "${snakemake_output[flanks]}"

ls -lh "${out_dir}"
