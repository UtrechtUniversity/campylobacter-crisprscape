#! /usr/bin/env bash
set -euo pipefail
IFS=$'\n'

downloaded_samples="data/tmp/all_downloaded_samples.txt"
other_species="data/tmp/other_species.txt"

# Get all accession IDs from all batches
for batch in data/tmp/ATB/batch_*
do
    ls ${batch} | sed 's/.fa//g'
done > ${downloaded_samples}

# Compare the accession IDs with those from your species of interest -
#  saving only the other ones (not of interest)
comm -23 <(sort ${downloaded_samples}) <(sort data/ATB/all_samples_of_interest.txt) > ${other_species}

echo -e "Removing $(wc -l ${other_species}) fasta files..."

# Use the list of non-Campylobacters to remove unnecessary fasta files
while read sample;
do
    rm data/tmp/ATB/batch_*/${sample}.fa
done < ${other_species}

echo "Done!"
