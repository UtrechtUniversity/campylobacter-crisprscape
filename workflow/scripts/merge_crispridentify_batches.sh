#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

cat ${snakemake_params[spacers_crispr]} > ${snakemake_output[spacers_crispr]}

for summary in ${snakemake_params[summary_crispr]}
do
    header=$(head -n 1 "${summary}")
    if [ "${header}" == "No arrays found" ]
    then
        continue
    else
        echo "${header}" | tee ${snakemake_output[summary_crispr]}
        break
    fi
done

for summary in ${snakemake_params[summary_crispr]}
do
    tail -n +2 "${summary}" >> ${snakemake_output[summary_crispr]}
done