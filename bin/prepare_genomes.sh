#!/usr/bin/env bash

## Prepare genome sequences for CRISPR analysis
##
## Steps:
##  1. Download metadata from AllTheBacteria
##  2. Extract sample names of species of interest
##  3. Filter metadata to species of interest
##  4. Find which batches contain species of interest
##  5. Download genome sequences (in batches) and extract XZ archives
##  6. Remove genomes from species other than species of interest
##
## Note: steps 1, 5 and 6 are separate scripts that are called within this
##  script. The species of interest is read from:
##  `config/species_of_interest.txt`.
## Also note that AllTheBacteria has had two releases as of December 2024.
##  It is possible to download 1) all genomes from both releases combined,
##  2) only genomes from the original release, or 3) genomes from the
##  incremental update (August 2024). This can be provided as input argument
##  for this script.
##
## Usage:
## $ bash prepare_genomes.sh [part]
##  (where [part] = 'all', 'original', 'update' (=default; smallest size))

part=${1:-"update"}

echo "Preparing to download genomes from set '${part}'"

## Step 1: download metadata
echo "Step 1: downloading metadata"
output_dir="data/ATB/"
bash bin/download_ATB_metadata.sh ${output_dir}

echo "Done downloading metadata!"
ls -lh ${output_dir}
echo "----------"

## Step 2: extract sample names of species of interest
#  use `config/species_of_interest.txt` to select for your species,
#  include only high-quality (remove "F" for fail on 'HQ' criterium),
#  and cut column 1 which contains the accession IDs
echo "Step 2: extracting sample accession IDs of species of interest"
species_samples_file="${output_dir}all_samples_of_interest.txt"
zgrep -f config/species_of_interest.txt ${output_dir}species_calls.tsv.gz |\
 grep -v -e "F$" | cut -f 1 > ${species_samples_file}

echo "Done extracting sample names!"
ls -lh ${species_samples_file}
echo "Contains: $(zless ${species_samples_file} | wc -l) entries"
echo "----------"

# Collect the number of genomes included for each species of interest
zgrep -f ${species_samples_file} -w ${output_dir}species_calls.tsv.gz |\
 cut -f 2 | sort | uniq -c | sort -nr > ${output_dir}number_of_genomes_per_species.txt

## Step 3: Filter metadata to the species of interest
echo "Step 3: Filtering metadata for species of interest"
zgrep -w -f ${species_samples_file} ${output_dir}ena_metadata.20240801.tsv.gz |\
 gzip > ${output_dir}ena_metadata.20240801-filtered.tsv.gz

echo "Done filtering!"
ls -lh ${output_dir}ena_metadata.20240801-filtered.tsv.gz
echo "----------"

## Step 4: Find which batches contain the species of interest
#  search for accession IDs and collect unique (deduplicated) file names, URLs and md5 checksums
echo "Step 4: Find relevant genome batches"
zgrep -w -f ${species_samples_file} ${output_dir}file_list.all.20240805.tsv.gz |\
 cut -f 5-7 | sort | uniq > ${output_dir}batches_to_download.tsv

echo "Found all relevant batches!"
ls -lh ${output_dir}batches_to_download.tsv
echo "Top 3 lines:"
echo -e "Filename\tDownload_URL\tmd5sum"
head -3 ${output_dir}batches_to_download.tsv
echo "----------"

## Step 5: Download genome sequences
echo "Step 5: Download genome sequences"
bash bin/download_genomes.sh ${part}
echo "Finished downloading!"
echo "The batches have been written to 'data/tmp/ATB/batch_*'"
echo "----------"

# Step 6: Remove genomes of other species
echo "Step 6: remove genomes of species other than species of interest"
echo "Files before filtering:"
for batch in data/tmp/ATB/batch_*
do
    echo "$(basename ${batch}):"
    ls ${batch} | wc -l
done

bash bin/remove_other_species.sh

echo "Files after filtering:"
for batch in data/tmp/ATB/batch_*
do
    echo "$(basename ${batch}):"
    ls ${batch} | wc -l
done
rmdir --ignore-fail-on-non-empty data/tmp/ATB/batch_*
echo "----------"
echo "Genome files have been prepared!"
