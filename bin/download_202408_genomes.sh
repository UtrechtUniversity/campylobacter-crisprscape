#! /usr/bin/env bash
set -euo pipefail
IFS=$'\n'

# Above thanks to Aaron Maxwell: http://redsymbol.net/articles/unofficial-bash-strict-mode/

for line in $(grep "incr_release.202408" data/ATB/Campylobacter_batches_to_download-20240918.tsv)
do
    filename=$(echo ${line} | cut -f 1)
    url=$(echo ${line} | cut -f 2)
    checksum=$(echo ${line} | cut -f 3)
    echo -e "Filename: ${filename}\tURL: ${url}\tmd5sum: ${checksum}"

    mkdir -p data/tmp/ATB

    outputfile="data/tmp/ATB/${filename}"

    if [ ! -f ${outputfile} ]
    then
        wget -O ${outputfile} ${url}
    fi

    if echo ${checksum} ${outputfile} | md5sum -c --status;
    then
        echo "OK: md5sum for ${outputfile} is correct!"
    else
        echo "BAD: md5sum for ${outputfile} is incorrect..."
    fi
done
