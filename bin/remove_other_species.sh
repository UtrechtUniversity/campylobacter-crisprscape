#! /usr/bin/env bash
set -euo pipefail
IFS=$'\n'

downloaded_samples="data/tmp/all_downloaded_samples.txt"
other_species_samples="data/tmp/other_species-samples.txt"

# Get all accession IDs from all batches
for batch in data/tmp/ATB/batch_*
do
    ls ${batch} | sed 's/.fa//g'
done > ${downloaded_samples}

# Compare the accession IDs with those from your species of interest -
#  saving only the other ones (not of interest)
comm -23 <(sort ${downloaded_samples}) <(sort data/ATB/all_samples_of_interest.txt) > ${other_species_samples}

echo -e "Removing $(wc -l ${other_species_samples}) fasta files..."

# Use the list of non-Campylobacters to remove unnecessary fasta files
while read sample;
do
    rm data/tmp/ATB/batch_*/${sample}.fa
done < ${other_species_samples}

# For the curious: make a list of species that were downloaded
zgrep -f ${other_species_samples} -w data/ATB/species_calls.tsv.gz |\
 tee data/tmp/other_species_calls.tsv |\
 cut -f 2 | sort | uniq -c | sort -nr > data/tmp/other_species_numbers.txt

echo "Done!"
