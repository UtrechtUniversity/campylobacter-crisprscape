#! /usr/bin/env bash
set -euo pipefail
IFS=$'\n'

# Get all accession IDs from all batches
for batch in data/tmp/ATB/batch_*
do
    ls ${batch} | sed 's/.fa//g'
done > data/tmp/ATB/all_samples.txt

# Compare the accession IDs with those from Campylobacter -
#  saving only the non-Campylobacter ones
comm -23 <(sort data/tmp/ATB/all_samples.txt) <(sort data/ATB/all_Campylobacter_samples-20240918.txt) > data/tmp/not_Campylobacter.txt

echo -e "Removing $(wc -l data/tmp/not_Campylobacter.txt) fasta files..."

# Use the list of non-Campylobacters to remove unnecessary fasta files
while read sample;
do
    rm data/tmp/ATB/batch_*/${sample}.fa
done < data/tmp/not_Campylobacter.txt

echo "Done!"
