#! /usr/bin/env bash

# Read the output directory from the command-line, with
#  'data/ATB/' set as default.
output_dir="${1:-"data/ATB"}/"

mkdir -p ${output_dir}

# Download metadata files for AllTheBacteria
# (https://allthebacteria.readthedocs.io/en/latest/)
# only when the file is not already there:
download_if_not_present() {
    output_file=$1
    url=$2
    if [ ! -s ${output_file} ]
    then
        echo "Downloading ${output_file}..."
        wget -O ${output_file} ${url}
    else
        echo "${output_file} downloaded before:"
        ls -lh ${output_file}
    fi
}

# 1. ENA metadata
download_if_not_present ${output_dir}ena_metadata.0.2.20240606.tsv.gz https://osf.io/download/j47ug/
download_if_not_present ${output_dir}ena_metadata.20240801.tsv.gz https://osf.io/download/ks7yt/

# 2. Sample lists
download_if_not_present ${output_dir}sample_list.txt.gz https://osf.io/download/ev3sj/

# 3. Sylph (species abundances per sample)
download_if_not_present ${output_dir}sylph.tsv.gz  https://osf.io/download/nu5a6/

# 4. Assembly statistics
download_if_not_present ${output_dir}assembly-stats.tsv.gz https://osf.io/download/nbyqv/

# 5. CheckM2 results
download_if_not_present ${output_dir}checkm2.tsv.gz https://osf.io/download/289f5/

# 6. Species calls
download_if_not_present ${output_dir}species_calls.tsv.gz https://osf.io/download/7t9qd/

# Get a list of all files:
download_if_not_present ${output_dir}file_list.all.20240805.tsv.gz https://osf.io/download/dw3h7
