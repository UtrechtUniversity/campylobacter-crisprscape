#!/bin/bash

##this script takes organised ENA metadata from organise_ATB_metadata.sh
##creating .fasta format files with the spacer sequence as the fasta sequence
##
## command structure:
## generate_fasta_from_spacer [spacer data file] [filepath(default=fasta_files)]
## [spacer data file] = tab separated file created from organise_ATB_metadata.sh
## [filepath] = folder where the fasta files will be created including name
### TODO: add fluidity to check where in the columns the sample name and the spacer sequence is present/

#gather structure of the data
columns=$(head -n 1 $1)
samples=$(cut -f 1 $1 | sed 1d | sort -u | head)
filepath=${2:-fasta_files}
#create fasta files folder for output
mkdir -p $filepath

echo "Creating fasta files in $filepath"
#create fasta files
for sample in $samples
do
    > $filepath/$sample.fasta
    batch=$(grep $sample $1)
    while IFS="" read -r line
        do 
        echo "$line" | tr "\t" ":" | cut -d ":" -f 1,3-4 | sed -e "s/^/>/" >> $filepath/$sample.fasta
        echo "$line" | cut -f 2 >> $filepath/$sample.fasta
        done <<< "$batch"       
done


#create metadata file
echo "Metadata file for fasta files in this folder.
The structure of the fasta identification is:" > $filepath/metadata.txt
echo "$columns" | cut -f 1,3-4 | tr "\t" ":" >> $filepath/metadata.txt
echo "Done!"