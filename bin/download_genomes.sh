#! /usr/bin/env bash
set -euo pipefail
IFS=$'\n'

# Above thanks to Aaron Maxwell: http://redsymbol.net/articles/unofficial-bash-strict-mode/

part=${1:-"update"}

if [ "${part}" == "update" ]
then
    download_list=$(grep "incr_release" data/ATB/batches_to_download.tsv)
elif [ "${part}" == "original" ]
then
    download_list=$(grep -v "incr_release" data/ATB/batches_to_download.tsv)
elif [ "${part}" == "all" ]
then
    download_list=$(cat data/ATB/batches_to_download.tsv)
else
    download_list=""
    echo "Unknown argument provided! Please use 'all', 'original', or 'update'."
    echo "Or none to use the default (=update)."
fi

for line in ${download_list}
do
    filename=$(echo ${line} | cut -f 1)
    url=$(echo ${line} | cut -f 2)
    checksum=$(echo ${line} | cut -f 3)
    echo -e "Filename: ${filename}\tURL: ${url}\tmd5sum: ${checksum}"

    mkdir -p data/tmp/ATB

    outputfile="data/tmp/ATB/${filename}"

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
    batchnumber=$(basename ${outputfile} | cut -f 6 -d '.')
    batchdir="data/tmp/ATB/batch_${batchnumber}"

    # If the batch directory has not been made yet
    if [ ! -d ${batchdir} ]
    then
        echo "Extracting ${outputfile}!"

        mkdir -p ${batchdir}

        # Decompress the XZ archive, send the output to the specified directory,
        # and strip the leading directory from within the XZ archive.
        tar -Jxf ${outputfile} -C ${batchdir} --strip-components 1

    else
        echo "Fasta files have been extracted previously!"
    fi
done
