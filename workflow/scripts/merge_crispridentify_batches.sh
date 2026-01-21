#!/usr/bin/env bash
set -euo pipefail

exec 2> "${snakemake_log[0]}" # send all stderr to the log file

spacers_crispr=(${snakemake_input[spacers_crispr]})
summary_crispr=(${snakemake_input[summary_crispr]})

> ${snakemake_output[spacers_crispr]}
for spacers in "${spacers_crispr[@]}"
do
    cat "${spacers}" >> ${snakemake_output[spacers_crispr]}
done

for summary in "${summary_crispr[@]}"
do
    header=$(head -n 1 "${summary}")
    if [ "${header}" == "No arrays found" ]
    then
        continue
    else
        echo "${header}" > ${snakemake_output[summary_crispr]}
        break
    fi
done

for summary in "${summary_crispr[@]}"
do
    tail -n +2 "${summary}" >> ${snakemake_output[summary_crispr]}
done
