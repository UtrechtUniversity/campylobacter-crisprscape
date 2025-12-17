#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n'

# Above thanks to Aaron Maxwell: http://redsymbol.net/articles/unofficial-bash-strict-mode/

part=${1:-"update"}
atb_dir=${2:-"resources/ATB/"}
output_dir="${work_dir}archives/"

if [ "${part}" == "update" ]
then
    download_list=$(grep "incr_release" ${atb_dir}/batches_to_download.tsv)

elif [ "${part}" == "original" ]
then
    download_list=$(grep -v "incr_release" ${atb_dir}/batches_to_download.tsv)

elif [ "${part}" == "all" ]
then
    download_list=$(cat ${atb_dir}/batches_to_download.tsv)

else
    download_list=""
    echo "Unknown argument provided! Please use 'all', 'original', or 'update'."
    echo "Or none to use the default (=update)."
fi

mkdir -p ${output_dir}

for line in ${download_list}
do
    filename=$(echo ${line} | cut -f 1 | sed -e 's/assembly/bakta/')
    url=$(grep ${filename} data/ATB/all_atb_files.tsv | cut -f 4)
    checksum=$(grep ${filename} data/ATB/all_atb_files.tsv | cut -f 5)
    echo -e "Filename: ${filename}\tURL: ${url}\tmd5sum: ${checksum}"

    outputfile="${output_dir}${filename}"

    # If the output file is not a file of size greater than zero
    if [ ! -s ${outputfile} ]
    then
        # Then download it
        wget -O ${outputfile} ${url}
    else
        echo "${outputfile} has been downloaded before!"
    fi

    # Check the md5sum
    if echo ${checksum} ${outputfile} | md5sum -c --status;
    then
        echo "OK: md5sum for ${outputfile} is correct!"
    else
        # If it is wrong, delete the file
        echo "BAD: md5sum for ${outputfile} is incorrect... Deleting file"
        rm ${outputfile}
    fi

    # Extract the batch number from the file name
    batchdir="${atb_dir}/annotations"
    mkdir -p ${batchdir}

    # If the batch directory has not been made yet
    if [ ! -d "${batchdir}${outputfile/.tar.xz/}" ]
    then
        echo "Extracting ${outputfile}!"

        # Decompress the XZ archive and send the output to the specified directory.
        tar -Jxf ${outputfile} -C ${batchdir}

    else
        echo "Bakta files have been extracted previously!"
    fi
done
