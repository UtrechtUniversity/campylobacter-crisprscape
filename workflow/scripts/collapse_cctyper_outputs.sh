#!/usr/bin/env bash
set -euo pipefail

# Get the working directory from the command-line
work_dir=$1

# Within it should be subdirectories with partial output
subdirectories=( $(find ${work_dir} -mindepth 1 -maxdepth 1 -type d) )

# Concatenate table files

### THIS CODE IS COPIED FROM THE 'concatenate_cctyper_output.sh' SCRIPT ###
concatenate_files () {
    file_basename=$1

    output_file="${work_dir}/${file_basename}"

    if [ ! -e ${output_file} ]
    then
        all_files=( $(find ${work_dir} -mindepth 2 -maxdepth 2 -type f -name "${file_basename}") )
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

# Concatenate hmmer.tab
concatenate_files "hmmer.tab"

# Already filter out the arrays that CCTyper identified:
#  remove the 'putatives'.
grep -f <(cut -f 2 ${work_dir}/CRISPR_Cas.tab)\
 ${work_dir}/cas_operons_putative.tab > ${work_dir}/cas_operons.tab

echo "Filtered the cas operons:"
ls -lh ${work_dir}/cas_operons.tab

### THE ABOVE CODE IS COPIED FROM THE 'concatenate_cctyper_output.sh' SCRIPT ###

mkdir -p ${work_dir}/spacers

for index in "${!subdirectories[@]}"
do
# Move figures separately, with an index number as suffix
    mv ${subdirectories[${index}]}/plot.svg ${work_dir}/plot_${index}.svg ||\
     echo "No plot generated in ${subdirectories[${index}]}"

# and move the spacer fasta files, too!
    mv ${subdirectories[${index}]}/spacers/*.fa ${work_dir}/spacers/ ||\
     echo "No spacers found in ${subdirectories[${index}]}"
done

# Remove all partial outputs in the end
find ${work_dir} -mindepth 1 -maxdepth 1 -type d -name "part*" -exec rm -rf {} \;
