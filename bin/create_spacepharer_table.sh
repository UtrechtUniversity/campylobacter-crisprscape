#!/usr/bin/env bash

## this script uses the output from the two spacepharer matches created from the spacepharer rules. This script assumes that Phagescope and PLSDB are the used databases.
## Overall the script takes every match and extracts the ID which is directly taking from the database and then tries to match this back to the metadata and attach this to the match.

#spacepharer does not include its own column names and expected metadata column names are also manually added
echo -e "sample_accession\tphage_accession\tp_best_hit\tspacer_start\tspacer_end\tphage_start\tphage_end\t5_3_PAM\t3_5_PAM\tLength\tGC_content\ttaxonomy\tcompleteness\thost\tlifestyle" > ${snakemake_output[phage]}
while read line; do
    ID=$(echo $line | cut -f 2)
    metadata_match=$(grep -w "$ID" ${snakemake_input[meta_phage]}/merged_metadata.tsv | cut -f 2-7)
    echo -e "$line\t$metadata_match" >> ${snakemake_output[phage]}
done < ${snakemake_input[phage]}

echo "starting plasmid"
echo -e "sample_accession\tphage_accession\tp_best_hit\tspacer_start\tspacer_end\tphage_start\tphage_end\t5_3_PAM\t3_5_PAM\ttaxonomy" > ${snakemake_output[plasmid]}
while read line; do
    ID=$(echo $line | cut -f 2)
    #PLSDB uses a metadata system where there are many different files for differing purposes. taxonomy.csv uses a different ID so the taxonomy ID needs to be collected from nuccore and then matched to taxonomy
    nuccore_match=$(grep -w "$ID" ${snakemake_input[meta_plasmid]}/nuccore.csv | cut -f 13 -d ",")
    taxonomy_match=$(grep -w "^$nuccore_match" ${snakemake_input[meta_plasmid]}/taxonomy.csv | cut -f 3 -d ",")
    echo -e "$line\t$taxonomy_match" >> ${snakemake_output[plasmid]}
done < ${snakemake_input[plasmid]}