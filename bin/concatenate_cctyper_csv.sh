#!/usr/bin/env bash

batch=$1
batch_name=$(basename ${batch})
output_dir=${2:-"none"}
filename="CRISPR-Cas.csv"

echo ${batch_name}

concatenate_files () {
    if [ ${output_dir} == "none" ]
    then
        output_file="${batch}/${filename/.csv/-${batch_name}.csv}"
    else
        mkdir -p ${output_dir}
        output_file="${output_dir}/${filename/.csv/-${batch_name}.csv}"
    fi

    if [ ! -e ${output_file} ]
    then
        all_files=( $(find ${batch} -mindepth 2 -maxdepth 2 -name "${filename}") )
        echo "Concatenating all ${#all_files[@]} $(basename ${filename}) files!"

        first_file=${all_files[1]}
        # First, only extract the header from the first file:
        head -n 1 ${first_file} > ${output_file}

        # Then, append the contents (without header!) of all files:
        for tab_file in ${all_files[*]}
        do
            sed 1d ${tab_file} >> ${output_file}
        # `sed 1d` deletes the first line and prints the rest!
        # (https://stackoverflow.com/a/34504648)
        done
        ls -lh ${output_file}
        echo

    else
        echo "$(basename ${output_file}) already concatenated!"
        ls -lh ${output_file}
        echo
    fi
}

# Concatenate all CRISPR-Cas.csv files in a batch
concatenate_files
