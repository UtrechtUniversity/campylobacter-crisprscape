#!/usr/bin/env bash

batch=$1
batch_name=$(basename ${batch})

echo ${batch_name}

concatenate_files () {
    file_basename=$1
    output_file="${batch}/${file_basename/.tab/-${batch_name}.tab}"

    if [ ! -e ${output_file} ]
    then
        all_files=( $(find ${batch} -mindepth 2 -maxdepth 2 -name "${file_basename}") )
        echo "Concatenating all ${#all_files[@]} $(basename ${file_basename}) files!"

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

# Concatenate CRISPR_Cas.tab
concatenate_files "CRISPR_Cas.tab"

# Concatenate crisprs_all.tab
concatenate_files "crisprs_all.tab"

# Concatenate crisprs_near_cas.tab
concatenate_files "crisprs_near_cas.tab"

# Concatenate crisprs_orphan.tab
concatenate_files "crisprs_orphan.tab"
