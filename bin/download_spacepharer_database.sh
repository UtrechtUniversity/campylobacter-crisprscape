#! bin/bash

## This script downloads the phage and plasmid databases: Phagescope (https://phagescope.deepomics.org/) and PLSDB (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/).
## They are then also extracted and specifically for Phagescope merged into one database for use in Spacepharer. 
##
## Note: RUN FROM THE BASE FOLDER AS IT WILL NOT WORK OTHERWISE
## Usage: download_spacepharer_database.sh [threads]
## [threads]: the amount of threads used for extracting and merging phagescope data. can at maximum use 14 threads, default is 1. 
threads="${1:-"1"}"
#PLSDB download
mkdir -p data/raw/PLSDB
if  [ -f data/raw/PLSDB/download_meta.tar.gz ]; then
    echo "Already downloaded PLSDB, skipping..."
else
    echo "Downloading PLSDB"
    wget -P data/raw/PLSDB https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz
fi

echo "Extracting download"
if [ -f data/raw/PLSDB/sequences.fasta.bz2 ]; then
    echo "already extracted, skipping..."
else
    tar -xzf data/raw/PLSDB/download_meta.tar.gz -C data/raw/PLSDB/
fi
echo "Unzipping sequences"
bzip2 -d data/raw/PLSDB/sequences.fasta.bz2

echo "correcting metadata delims"
sed -i -E ':a;s/"([^"]*),([^"]*)"/"\1\2"/g;ta' nuccore.csv

#Phagescope download
mkdir -p data/raw/phagescope

echo "Downloading Phagescope databases"

for DB in "Genbank" "RefSeq" "DDBJ" "EMBL" "PhagesDB" "GPD" "GVD" "MGV" "TemPhD" "CHVD" "IGVD" "IMG_VR" "GOV2" "STV"
    do
    if test -f "data/raw/phagescope/$DB.tar.gz"; then
        echo "Already downloaded $DB, skipping..."
    else
        echo "Downloading $DB"
        wget -O "data/raw/phagescope/$DB.tar.gz"  "https://phageapi.deepomics.org/download/phage/fasta/?datasource=$DB"
    fi
done

echo "Extracting databases"
parallel --jobs "$threads" 'DB={}; path=data/raw/phagescope; [ -f "$path/$DB.fasta" ] || [ -d "$path/$DB" ] && \
echo "$DB already extracted, skipping..." || \
( echo "extracting $DB"; tar -xzf "$path/$DB.tar.gz" -C "$path/" )' ::: Genbank RefSeq DDBJ EMBL PhagesDB GPD GVD MGV TemPhD CHVD IGVD IMG_VR GOV2 STV

echo "merging phagescope sequences"
parallel --jobs "$threads" 'DB={}; path=data/raw/phagescope; [ -f "${path}/${DB}.fasta" ] && \
echo "$DB already merged, skipping..." || \
( echo "merging $DB"; genomes=$(find "$path/$DB" -type f -name "*.fasta"); > $path/$DB.fasta ; for files in $genomes; do cat $files >> $path/$DB.fasta; done)' \
::: Genbank RefSeq DDBJ EMBL PhagesDB GPD GVD MGV TemPhD


echo "Downloading Phagescope metadata"

for DB in "genbank" "refseq" "ddbj" "embl" "phagesdb" "gpd" "gvd" "mgv" "temphd" "chvd" "igvd" "img_vr" "gov2" "stv"
    do
    if [ -f "data/raw/phagescope/${DB}_phage_meta_data.tsv" ]; then
        echo "Already downloaded $DB, skipping..."
    else
        echo "Downloading $DB"
        wget -P data/raw/phagescope/ "https://phageapi.deepomics.org/files/Download/Phage_meta_data/${DB}_phage_meta_data.tsv"
    fi
done

echo "Merging phagescope metadata"
first=true
if test -f "data/raw/phagescope/merged_metadata.tsv"; then
    echo "Already merged metadata, skipping..."
else
    for DB in "genbank" "refseq" "ddbj" "embl" "phagesdb" "gpd" "gvd" "mgv" "temphd" "chvd" "igvd" "img_vr" "gov2" "stv";
        do
        echo "Merging $DB"
        if [ "$first" == true ]; then
            cat "data/raw/phagescope/${DB}_phage_meta_data.tsv" > data/raw/phagescope/merged_metadata.tsv
            first=false
        else
            tail -n +2 "data/raw/phagescope/${DB}_phage_meta_data.tsv" >> data/raw/phagescope/merged_metadata.tsv
        fi
    done 
fi
echo "Done!"
