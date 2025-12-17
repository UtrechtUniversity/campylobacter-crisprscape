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
## $ bash prepare_genomes.sh -p/--part [part] -d/--directory [directory] -h/--help
##  (where [part] = 'all', 'original', 'update' (=default; smallest size))

# Set global (default) variables
part="update"
output_dir="resources/ATB/"
species_of_interest="config/species_of_interest.txt"
genomes_of_interest="${output_dir}stats/total_genomes_of_interest.tsv"
bakta_dir="${output_dir}annotations/"

print_help() {
  # Display help for this script
  echo "Prepare genomes from AllTheBacteria for use with CRISPRscape"
  echo
  echo "Syntax: prepare_genomes.sh -p [part] -d [directory] [-h]"
  echo "Options:"
  echo "-p/--part      Select which part of AllTheBacteria to download"
  echo "               for the selected species ('all', 'original', or"
  echo "               'update', default=update)"
  echo "-d/--directory Directory in which to download the files"
  echo "               (default=resources/ATB/)"
  echo "-h/--help      Print this help message"
  echo
  exit 0
}

while [[ $# -gt 0 ]]
do
  case "$1" in
    -p|--part ) part="$2"
    shift
    shift
    ;;
    -d|--directory ) output_dir="$2"
    shift
    shift
    ;;
    -h|--help ) print_help
    ;;
    * )
    echo "Unknown option: $1"
    print_help
  esac
done

echo "Preparing to download genomes from set '${part}'"
echo "And storing them in the directory: ${output_dir}."
echo "----------"

## Step 1: download metadata
echo "Step 1: downloading metadata"
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

echo -e "Species of interest:\n$(cat ${species_of_interest})"

total_genomes=$(zgrep -f ${species_of_interest} ${output_dir}species_calls.tsv.gz | wc -l)
hq_genomes=$(zgrep -f ${species_of_interest} ${output_dir}species_calls.tsv.gz | grep "T$" | wc -l)
lq_genomes=$(zgrep -f ${species_of_interest} ${output_dir}species_calls.tsv.gz | grep "F$" | wc -l)

mkdir -p "${output_dir}/stats"

echo -e "Total_genomes\t${total_genomes}" > "${genomes_of_interest}"
echo -e "High-quality_genomes\t${hq_genomes}" >> "${genomes_of_interest}"
echo -e "Low-quality_genomes\t${lq_genomes}" >> "${genomes_of_interest}"

echo "ATB contains ${total_genomes} genomes of your species of interest."
echo "Of those, ${lq_genomes} are labeled as low-quality, which are not included for further analyses."
echo "That means, ${hq_genomes} are available to work with."

echo
echo "Also see the file ${genomes_of_interest} to see these"
echo "numbers in table form (tab-separated values)."
echo

# Use no further grep options to match anything that contains the species names,
# including subspecies and lineages. Exclude low-quality 'HQ field == F'.
zgrep -f ${species_of_interest} ${output_dir}species_calls.tsv.gz |\
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

echo -e "Genomes_in_original_release\t${original_genomes}" >> "${genomes_of_interest}"
echo -e "Genomes_in_update_release\t${update_genomes}" >> "${genomes_of_interest}"

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
 ${output_dir}stats/other_genomes-numbers.txt

# And get the list of batches + sample accessions, to facilitate removal:
zgrep -w -f ${output_dir}samples_not_of_interest.txt\
 ${output_dir}file_list.all.20240805.tsv.gz |\
 cut -f 4 > ${output_dir}samples_to_remove.txt

echo "In the batches to download are $(wc -l ${output_dir}samples_not_of_interest.txt)"
echo "samples that have low-quality genomes or species other than the species of interest."
echo "These are summarised in:"
ls -lh ${output_dir}samples_not_of_interest.txt ${output_dir}stats/other_genomes-numbers.txt ${output_dir}samples_to_remove.txt

## Step 5: Download genome sequences
echo "Step 5: Download genome sequences"
bash bin/download_genomes.sh ${part} ${output_dir}
echo "Finished downloading!"
echo "The batch archives have been written to '${output_dir}'"
echo "and their contents extracted to '${output_dir}assemblies/'"
echo "----------"

## Step 6: Remove genomes of other species
echo "Step 6: remove genomes of species other than species of interest"
for fasta in $(cat ${output_dir}samples_to_remove.txt)
do
# Verbose remove: tell what is being removed
    rm -fv "${output_dir}assemblies/${fasta}"
done
echo "----------"

## Step 7: Download functional annotations (Bakta)
echo "Step 7: download functional (gene) annotations"
bash bin/download_bakta_annotations.sh ${part} ${output_dir}
echo "Finished downloading!"
echo "The batches have been downloaded to '${output_dir}' and extracted in '${output_dir}annotations'"
echo "----------"

## Step 7.1: Also remove annotation files for non-of-interest samples
echo "Step 7.1: remove genome annotations of other/low-quality samples"
# adjust the file basename from 'assembly' to 'bakta'
for file in $(sed 's|atb.assembly.|atb.bakta.|g' ${output_dir}samples_to_remove.txt)
do
    # add 'annotations' subdirectory and exchange '.fa' extension to 'bakta.json'
    json="${bakta_dir}${file/.fa/.bakta.json}"
    rm -fv ${json}
done
