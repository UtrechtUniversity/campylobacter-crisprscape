#! bin/bash

## This script downloads the phage and plasmid databases: Phagescope (https://phagescope.deepomics.org/) and PLSDB (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/).
## They are then also extracted and specifically for Phagescope merged into one database for use in Spacepharer. 
##
## Note: RUN FROM THE BASE FOLDER AS IT WILL NOT WORK OTHERWISE
## Usage: download_spacepharer_database.sh
## TODO: only uses one thread and can take very long for phagescope (unzipping the very large tar files), need to parallelise.

#PLSDB download
mkdir -p data/raw/PLSDB
if test -f data/raw/PLSDB/download_meta.tar.gz; then
    echo "Already downloaded PLSDB, skipping..."
else
    echo "Downloading PLSDB"
    wget -P data/raw/PLSDB https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz
fi

echo "Extracting download"
tar -xzf data/raw/PLSDB/download_meta.tar.gz -C data/raw/PLSDB/

echo "Unzipping sequences"
bzip2 -d data/raw/PLSDB/sequences.fasta.bz2

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

for DB in "Genbank" "RefSeq" "DDBJ" "EMBL" "PhagesDB" "GPD" "GVD" "MGV" "TemPhD" "CHVD" "IGVD" "IMG_VR" "GOV2" "STV"
    do
    echo "extracting $DB"
    tar -xzf "data/raw/phagescope/$DB.tar.gz" -C data/raw/phagescope/

done

echo "Downloading Phagescope metadata"

for DB in "genbank" "refseq" "ddbj" "embl" "phagesdb" "gpd" "gvd" "mgv" "temphd" "chvd" "igvd" "img_vr" "gov2" "stv"
    do
    if test -f "data/raw/phagescope/$DB_phage_meta_data.tsv"; then
        echo "Already downloaded $DB, skipping..."
    else
        echo "Downloading $DB"
        wget -P data/raw/phagescope/ "https://phageapi.deepomics.org/files/Download/Phage_meta_data/${DB}_phage_meta_data.tsv"
    fi
done

echo "Merging phagescope metadata"
first=true
for DB in "genbank" "refseq" "ddbj" "embl" "phagesdb" "gpd" "gvd" "mgv" "temphd" "chvd" "igvd" "img_vr" "gov2" "stv"
    do
    if test -f "data/raw/phagescope/merged_metadata.tsv"; then
        echo "Already merged metadata, skipping..."
    else
        echo "Merging"
        if [ "$first" == true ]; then
            cat $DB > data/raw/phagescope/merged_metadata.tsv
            first=false
        else
            tail -n +2 $DB >> data/raw/phagescope/merged_metadata.tsv
        fi
done 
echo "Done!"