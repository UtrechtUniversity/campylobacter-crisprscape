#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n'

# Above thanks to Aaron Maxwell: http://redsymbol.net/articles/unofficial-bash-strict-mode/

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

echo -e "Species of interest:\n$(cat config/species_of_interest.txt)"

total_genomes=$(zgrep -f config/species_of_interest.txt ${output_dir}species_calls.tsv.gz | wc -l)
hq_genomes=$(zgrep -f config/species_of_interest.txt ${output_dir}species_calls.tsv.gz | grep "T$" | wc -l)
lq_genomes=$(zgrep -f config/species_of_interest.txt ${output_dir}species_calls.tsv.gz | grep "F$" | wc -l)

echo -e "Total_genomes\t${total_genomes}" > data/tmp/total_genomes_of_interest.tsv
echo -e "High-quality_genomes\t${hq_genomes}" >> data/tmp/total_genomes_of_interest.tsv
echo -e "Low-quality_genomes\t${lq_genomes}" >> data/tmp/total_genomes_of_interest.tsv

echo "ATB contains ${total_genomes} genomes of your species of interest."
echo "Of those, ${lq_genomes} are labeled as low-quality, which are not included for further analyses."
echo "That means, ${hq_genomes} are available to work with."

echo
echo "Also see the file data/tmp/total_genomes_of_interest.tsv to see these"
echo "numbers in table form (tab-separated values)."
echo

# Use no further grep options to match anything that contains the species names,
# including subspecies and lineages. Exclude low-quality 'HQ field == F'.
zgrep -f config/species_of_interest.txt ${output_dir}species_calls.tsv.gz |\
 grep -v -e "F$" | cut -f 1 > ${species_samples_file}

echo "Done extracting sample names!"
ls -lh ${species_samples_file}
echo "Contains: $(wc -l ${species_samples_file}) entries"
echo "----------"

# Collect the number of genomes included for each species of interest
zgrep -f ${species_samples_file} -w ${output_dir}species_calls.tsv.gz |\
 cut -f 2 | sort | uniq -c | sort -nr > ${output_dir}number_of_genomes_per_species.txt

# Count the number of genomes of interest in the 'original' and 'update' releases:
original_genomes=$(zgrep -w -f ${species_samples_file}\
 ${output_dir}file_list.all.20240805.tsv.gz | cut -f 5 |\
  grep -c "r0.2")
update_genomes=$(zgrep -w -f ${species_samples_file}\
 ${output_dir}file_list.all.20240805.tsv.gz | cut -f 5 |\
  grep -c "incr_release")

echo -e "Genomes_in_original_release\t${original_genomes}" >> data/tmp/total_genomes_of_interest.tsv
echo -e "Genomes_in_update_release\t${update_genomes}" >> data/tmp/total_genomes_of_interest.tsv

## Step 3: Filter metadata to the species of interest
echo "Step 3: Filtering metadata for species of interest"
# Include the header for simpler processing in e.g. R:
zless ${output_dir}ena_metadata.20240801.tsv.gz | head -1 | gzip >\
 ${output_dir}ena_metadata.20240801-filtered.tsv.gz || true
# N.B. the 'or true' is inserted to avoid the non-zero exit status of head

# Use 'grep -w' to match only whole words; don't match longer accession numbers that
#  contain the sample accession of interest!
zgrep -w -f "${species_samples_file}" "${output_dir}ena_metadata.20240801.tsv.gz" |\
 gzip >> ${output_dir}ena_metadata.20240801-filtered.tsv.gz

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

## Step 4.2: Find out what genomes are in these batches and make a list of
#  genomes to exclude (low-quality and other species).
#  To get there, first collect all samples in the batches to download:
zgrep -f ${output_dir}batches_to_download.tsv ${output_dir}file_list.all.20240805.tsv.gz |\
 cut -f 1 > ${output_dir}samples_in_batches.txt

# Now take out the samples of interest:
grep -v -w -f ${species_samples_file} ${output_dir}samples_in_batches.txt >\
 ${output_dir}samples_not_of_interest.txt

# And for these 'not of interest genomes', give a summary of their
# numbers:
zgrep -w -f ${output_dir}samples_not_of_interest.txt\
 ${output_dir}file_list.all.20240805.tsv.gz | cut -f 2 | sort | uniq -c | sort -nr >\
 data/tmp/other_genomes-numbers.txt

# And get the list of batches + sample accessions, to facilitate removal:
zgrep -w -f ${output_dir}samples_not_of_interest.txt\
 ${output_dir}file_list.all.20240805.tsv.gz |\
 cut -f 4 > data/tmp/samples_to_remove.txt

echo "In the batches to download are $(wc -l ${output_dir}samples_not_of_interest.txt)"
echo "samples that have low-quality genomes or species other than the species of interest."
echo "These are summarised in:"
ls -lh ${output_dir}samples_not_of_interest.txt data/tmp/other_genomes-numbers.txt data/tmp/samples_to_remove.txt

## Step 5: Download genome sequences
echo "Step 5: Download genome sequences"
bash bin/download_genomes.sh ${part}
echo "Finished downloading!"
echo "The batch archives have been written to 'data/tmp/ATB/'"
echo "and their contents extracted to 'data/tmp/assemblies/'"
echo "----------"

## Step 6: Remove genomes of other species
echo "Step 6: remove genomes of species other than species of interest"
for fasta in $(cat data/tmp/samples_to_remove.txt)
do
# Verbose remove: tell what is being removed
    rm -fv "data/tmp/assemblies/${fasta}"
done
echo "----------"

## Step 7: Download functional annotations (Bakta)
echo "Step 7: download functional (gene) annotations"
bash bin/download_bakta_annotations.sh ${part}
echo "Finished downloading!"
echo "The batches have been downloaded to 'data/tmp/ATB/' and extracted in 'data/tmp/annotations'"
echo "----------"

## Step 7.1: Also remove annotation files for non-of-interest samples
# 1: adjust the file basename from 'assembly' to 'bakta'
echo "Step 7.1: remove genome annotations of other/low-quality samples"
for file in $(sed 's|atb.assembly.|atb.bakta.|g' data/tmp/samples_to_remove.txt)
do
    # 2: add 'annotations' subdirectory
    bakta_dir="data/tmp/annotations/"
    # 3: exchange '.fa' extension to 'bakta.json'
    json="${bakta_dir}${file/.fa/.bakta.json}"
    rm -fv ${json}
done
