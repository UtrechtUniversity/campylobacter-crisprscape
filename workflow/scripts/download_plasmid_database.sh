#!/usr/bin/env bash
## This script downloads and prepares the plasmid database: PLSDB
## (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) for use with SpacePHARER.
set -euo pipefail

exec > "${snakemake_log[out]}" # send all stdout to a log file
exec 2> "${snakemake_log[err]}" # send stderr to separate log file

. "${snakemake[scriptdir]}/utils.sh"

plasmid_dir=${snakemake_output[plasmid_dir]}
plasmid_archive="${plasmid_dir}/download_meta.tar.gz"

message "Downloading PLSDB plasmid database"
echo "-----"

wget -P "${plasmid_dir}" https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz

message "Extracting PLSDB"
tar -xzf "${plasmid_archive}" -C "${plasmid_dir}"

message "Extracting sequences"
bzip2 -d "${plasmid_dir}/sequences.fasta.bz2"

message "Adjusting metadata delimiters"
sed -i -E ':a;s/"([^"]*),([^"]*)"/"\1\2"/g;ta' "${plasmid_dir}/nuccore.csv"

echo -e "-----\nDone!\n-----"
