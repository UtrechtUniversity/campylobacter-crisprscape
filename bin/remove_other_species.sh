#! /usr/bin/env bash
set -euo pipefail
IFS=$'\n'

# Get all accession IDs from all batches
for batch in data/tmp/ATB/batch_*
do
    ls ${batch} | sed 's/.fa//g'
done > data/tmp/ATB/all_samples.txt

# Compare the accession IDs with those from your species of interest -
#  saving only the other ones (not of interest)
comm -23 <(sort data/tmp/ATB/all_samples.txt) <(sort data/ATB/all_samples_of_interest.txt) > data/tmp/other_species.txt

echo -e "Removing $(wc -l data/tmp/other_species.txt) fasta files..."

# Use the list of non-Campylobacters to remove unnecessary fasta files
while read sample;
do
    rm data/tmp/ATB/batch_*/${sample}.fa
done < data/tmp/other_species

echo "Done!"
