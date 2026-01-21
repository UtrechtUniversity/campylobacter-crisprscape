#!/usr/bin/env bash

head -1 ${snakemake_input[0]} > ${snakemake_output[0]}

for csv in ${snakemake_input[*]}
do
    sed 1d ${csv} >> ${snakemake_output[0]}
# `sed 1d` deletes the first line and prints the rest!
# (https://stackoverflow.com/a/34504648)
done
