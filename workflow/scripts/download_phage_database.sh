#!/usr/bin/env bash
## This script downloads and prepares the phage database: Phagescope
## (https://phagescope.deepomics.org/) for use with SpacePHARER.
set -euo pipefail

exec > "${snakemake_log[out]}" # send stdout to a log file
exec 2> "${snakemake_log[err]}" # also send stderr to a log file

. "${snakemake[scriptdir]}/utils.sh"

threads=${snakemake[threads]}
phage_dir="${snakemake_output[phage_dir]}"
databases=(${snakemake_params[databases]})

message "Downloading PhageScope databases (archived fasta files)"
message "This comprises ${#databases[@]} parts, as specified in 'config/parameters.yaml'"
echo "-----"

for DB in "${databases[@]}"
do
    message "Downloading ${DB}"
    wget -O "${phage_dir}/${DB}.tar.gz"\
     "https://phageapi.deepomics.org/download/phage/fasta/?datasource=${DB}"
done

echo "-----"
message "Extracting databases\n"

# Extracts the .tar.gz archives to subdirectories with separate fasta files.

parallel --jobs "${threads}" 'DB={1}; phage_dir={2};\
 [ -f "${phage_dir}/${DB}.fasta" ] || [ -d "${phage_dir}/${DB}" ] && \
 message "${DB} already extracted, skipping..." || \
 ( message "Extracting ${DB}"; tar -xzf "${phage_dir}/${DB}.tar.gz" -C "${phage_dir}/" )'\
 ::: "${databases[@]}" ::: ${phage_dir}

echo "-----"
message "Concatenating PhageScope sequences\n"

# Concatenate separate fasta files in one long file,
# and then remove the directory with the separate files.

parallel --jobs "${threads}" 'DB={1}; phage_dir={2};\
 ( message "Concatenating ${DB}"; genomes=$(find "${phage_dir}/${DB}" -type f -name "*.fasta");\
  > "${phage_dir}/${DB}.fasta";\
 for files in ${genomes}; do cat ${files} >> "${phage_dir}/${DB}.fasta"; done;
 rm -r "${phage_dir}/${DB}")' \
 ::: Genbank RefSeq DDBJ EMBL PhagesDB GPD GVD MGV TemPhD ::: ${phage_dir}

echo "-----"
message "Downloading Phagescope metadata\n"

# Download metadata as tab-separated values (TSV) text file.

for DB in "${databases[@]}"
do
    metadata_file="${phage_dir}/${DB}_phage_metadata.tsv"
    message "Downloading ${DB}"
    wget -O "${metadata_file}" "https://phageapi.deepomics.org/files/Download/Phage_meta_data/${DB,,}_phage_meta_data.tsv"
    # With ${DB,,} the string variable is converted to lowercase
done

echo "-----"
message "Concatenating PhageScope metadata files\n"

# Make one long file with all metadata from the different databases.
combined_metadata="${snakemake_output[combined_meta]}"

# Rename the first metadata file into the concatenated file
mv "${phage_dir}/${databases[0]}_phage_metadata.tsv" ${combined_metadata}

# Then, for all the other metadata files
for DB in "${databases[@]:1}"
do
    metadata_file="${phage_dir}/${DB}_phage_metadata.tsv"
    # Take all but the first line (header)
    tail -n +2 "${metadata_file}" >> "${combined_metadata}"
    # and remove (clean-up, avoid keeping duplicated data)
    rm ${metadata_file}
done

message "--- Done! ---"
