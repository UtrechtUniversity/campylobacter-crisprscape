#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

first=True
for summary in ${snakemake_input[cctyper]}
do
    if [ "${first}" == True ]
    then
        cat "${summary}" > tmp_file1
        first=False
    else
        tail -n +2 "${summary}" >> tmp_file1
    fi
done

header=$(head -n 1 ${snakemake_input[identify]} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t")
tail -n +2 ${snakemake_input[identify]} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t" > tmp_file2
first=True

while read line
do
    if [ "${first}" == True ]
    then
        first=False
        echo -e "$line\t$header" > ${snakemake_output[table]}
    else
        sample=$(echo -e "${line}" | cut -f 1)
        start_cc=$(echo -e "${line}" | cut -f 3)
        start_id=$(expr "${start_cc}" + 1)
        match=$(grep "${sample}_${start_id}" tmp_file2 || true)

        if [ -z "${match}" ]
        then
            echo -e "${line}" >> ${snakemak_output[table]}
        else
            while read line2
            do
                if [ "${start_cc}" -lt 5000 ];
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
rm -f tmp_file1 tmp_file2