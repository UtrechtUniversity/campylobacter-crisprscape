#!/usr/bin/env bash

batch=$1
batch_name=$(basename ${batch})
output_dir=${2:-"none"}

echo ${batch_name}

concatenate_files () {
    file_basename=$1
    if [ ${output_dir} == "none" ]
    then
        output_file="${batch}/${file_basename/.tab/-${batch_name}.tab}"
    else
        mkdir -p ${output_dir}
        output_file="${output_dir}/${file_basename/.tab/-${batch_name}.tab}"
    fi

    if [ ! -e ${output_file} ]
    then
        all_files=( $(find ${batch} -mindepth 2 -maxdepth 2 -name "${file_basename}") )
        echo "Concatenating all ${#all_files[@]} $(basename ${file_basename}) files!"

        first_file=${all_files[0]}
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

# Concatenate cas_operons_putative.tab
concatenate_files "cas_operons_putative.tab"
# Already filter out the arrays that CCTyper identified:
#  remove the 'putatives'.
if [ ${output_dir} == "none" ]
then
    grep -f <(cut -f 2 ${batch}/CRISPR_Cas-${batch_name}.tab)\
     ${batch}/cas_operons_putative-${batch_name}.tab\
     > ${batch}/cas_operons-${batch_name}.tab

    echo "Filtered the cas operons:"
    ls -lh ${batch}/cas_operons-${batch_name}.tab
else
    grep -f <(cut -f 2 ${output_dir}/CRISPR_Cas-${batch_name}.tab)\
     ${output_dir}/cas_operons_putative-${batch_name}.tab\
     > ${output_dir}/cas_operons-${batch_name}.tab

    echo "Filtered the cas operons:"
    ls -lh ${output_dir}/cas_operons-${batch_name}.tab
fi

