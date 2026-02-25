#!/usr/bin/env bash

# Reconstruct the GFF and TSV file from Bakta using the tool 'bakta_io'
# to read the information from (ATB's) JSON files.

set -euo pipefail

exec > "${snakemake_log[out]}" # send stdout to a log file
exec 2> "${snakemake_log[err]}" # also send stderr to a log file

. "${snakemake[scriptdir]}/utils.sh"

threads=${snakemake[threads]}
batch_dir="${snakemake_input[batch_dir]}"

mkdir -p "${batch_dir}/gff"
mkdir -p "${batch_dir}/tsv"

message "Converting JSON files from ${batch_dir} to all Bakta file formats,"
message "and moving GFF3 and TSV so subfolders ${batch_dir}/gff and ${batch_dir}/tsv"
echo

find -L ${batch_dir} -mindepth 1 -maxdepth 1 -type f -name "*.bakta.json" -print0 |\
 parallel -0 --jobs ${threads} --retry-failed --halt='now,fail=1'\
 'sample=$(basename -s .bakta {/.});\
 bakta_io --output "{//}/${sample}" --prefix ${sample} {};\
 mv "{//}/${sample}/${sample}.gff3" "{//}/gff/";\
 mv "{//}/${sample}/${sample}.tsv" "{//}/tsv/";\
 rm -r "{//}/${sample}"'

echo
message "--- Done! ---"
