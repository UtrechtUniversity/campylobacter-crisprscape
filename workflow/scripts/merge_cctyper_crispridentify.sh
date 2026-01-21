#!/usr/bin/env bash
set -euo pipefail

exec 2> "${snakemake_log[0]}" # send all stderr to the log file

cctyper=(${snakemake_input[cctyper]})

first=True
for summary in "${cctyper[@]}"
do
    # If it is the first, copy the whole file (including header)
    if [ "${first}" == True ]
    then
        cat "${summary}" > tmp_file1
        first=False
    # otherwise, only take the 'contents', without header
    else
        tail -n +2 "${summary}" >> tmp_file1
    fi
    # And write it all in one concatenated temporary file: 'tmp_file1'
done

# For CRISPRidentify, collect the header as variable
header=$(head -n 1 ${snakemake_input[identify]} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t")
# and write the contents of the 'complete summary' to another temporary file: 'tmp_file2'
tail -n +2 ${snakemake_input[identify]} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t" > tmp_file2

first=True

# Start reading 'tmp_file1' = CCTyper CRISPRs
while read line
do
    if [ "${first}" == True ]
    # For the first line...
    then
        first=False
        # ...concatenate the header with CRISPRidentify's header
        echo -e "$line\t$header" > ${snakemake_output[table]}
    else

    # For all other lines,
        sample=$(echo -e "${line}" | cut -f 1)
        start_cc=$(echo -e "${line}" | cut -f 3)
        start_id=$(expr "${start_cc}" + 1)
        # see if the sample has a match in CRISPRidentify (tmp_file2)
        match=$(grep "${sample}_${start_id}" tmp_file2 || true)

        if [ -z "${match}" ]
        then
        # If there's *no* match (-z = lenth 0)
            echo -e "${line}" >> ${snakemak_output[table]}
            # Just copy the line from CCTyper
        else
        # but if there is a match,
            while read line2
            do
                if [ "${start_cc}" -lt 5000 ];
                # check if the reported start position is greater than 5000 (the length of the flanking regions),
                # and match the lines from CCTyper (line) and CRISPRidentify (match/line2) accordingly
                then
                    echo -e "${line}\t${match}" >> ${snakemake_output[table]}
                else
                    start=$(echo -e "${line2}" | cut -f 2)
                    start=$(expr "${start}" + "${start_cc}" - 5000)
                    length=$(echo -e "${line2}" | cut -f 4)
                    end=$(expr "${length}" + "${start}" - 1)
                    begin=$(echo -e "${line2}" | cut -f 1)
                    rest=$(echo -e "${line2}" | cut -f 4-9)
                    echo -e "${line}\t${begin}\t${start}\t${end}\t${rest}" >> ${snakemake_output[table]}
                fi
            done <<< "${match}"
        fi
    fi
done < tmp_file1

# Remove the temporary files
rm -f tmp_file1 tmp_file2
