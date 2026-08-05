#!/usr/bin/env bash

set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1 # send stdout+strerr to a log file
. "${snakemake[scriptdir]}/utils.sh"

fasta_file="${snakemake_input[genomes]}"
flank_length=${snakemake_params[flank]}

###

create_fasta() {
    bed_file=$1
    output_fasta=$2
    flank=${3:-0}

    seqkit subseq --bed ${bed_file} -U ${fasta_file} -u ${flank} -d ${flank} |
    # Remove the original ID with positions, up to the space,
    #  leaving only the CRISPR ID.
    # Also replace dot by dash. (CRISPRidentify likes that)
     seqkit replace -p ".+\s" | seqkit replace -p "\." -r "-"\
     -o ${output_fasta}
}

message "Extracting CRISPR-Cas from ${fasta_file}..."

create_fasta "${snakemake_input[crisprcas]}" "${snakemake_output[crispr_cas]}"
create_fasta "${snakemake_input[crisprcas]}" "${snakemake_output[cc_flanks]}" ${flank_length}

###

message "Extracting all CRISPR arrays from ${fasta_file}... (Including from CRISPR-Cas and orphans!)"

create_fasta "${snakemake_input[crispr]}" "${snakemake_output[crispr]}"
create_fasta "${snakemake_input[crispr]}" "${snakemake_output[cr_flanks]}" ${flank_length}

###

message "Extracting orphan CRISPR arrays from ${fasta_file}..."

create_fasta "${snakemake_input[orphan]}" "${snakemake_output[orphan]}"
create_fasta "${snakemake_input[orphan]}" "${snakemake_output[or_flanks]}" ${flank_length}

###

message "Extracting separate cas operons from ${fasta_file}..."

create_fasta "${snakemake_input[cas]}" "${snakemake_output[cas]}"
create_fasta "${snakemake_input[cas]}" "${snakemake_output[ca_flanks]}" ${flank_length}

###

message "--- Done! ---"
