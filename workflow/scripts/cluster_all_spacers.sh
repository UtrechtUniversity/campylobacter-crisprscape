#!/usr/bin/env bash

###########################################################
#                                                         #
#              --- CRISPR clustering script ---           #
#                                                         #
# We have detected many CRISPR spacers in Campylobacter   #
# genomes using CCTyper. Sometimes these are one nucleo-  #
# tide longer or shorter than otherwise identical coun-   #
# terparts. Therefore, it would be nice to cluster CRISPR #
# spacers based on their similarity and present one con-  #
# sensus sequence as representative for each CRISPR spa-  #
# cer. To do this, I want to cluster all spacers using    #
# CD-HIT-EST with decreasing similarity scores while      #
# making sure that the shorter spacer is always complete- #
# ly covered in the pairwise alignment. Then, to evaluate #
# the clustering cut-offs, I want to plot the number of   #
# clusters for each cut-off in a point (or bar) plot.     #
# This requires:                                          #
#   1. A script to run CD-HIT-EST in steps from 100%      #
#      identity down to ~60%? in steps of 3 (the spacers  #
#      are ~30nt long, so 3% difference ~1nt).            #
#   2. A secondary script to count/parse the number of    #
#      clusters for each cut-off value.                   #
#                                                         #
###########################################################

### Step 1: Cluster using decreasing similarity cut-offs.

all_spacers=$1
output_dir=$2
log_dir=$3

mkdir -p ${log_dir}

for cutoff in 1 0.96 0.93 0.9 0.87 0.84 0.81
# 0.8 is the lowest that CD-HIT-EST can take
do
    cd-hit-est -AS 0 -r 1 -sf 1 -c ${cutoff} -i ${all_spacers}\
     -o ${output_dir}/all_spacers-clustered-${cutoff}\
      > ${log_dir}/all_spacers-clustered-${cutoff}.txt 2>&1
done
# AS 0: match the entire shorter sequence,
# r 1: include reverse complements (default)
# sf 1: sort fasta by cluster size (large to small)
# c [cutoff]: sequence identity threshold
#   (global sequence identity, as fraction of full length of shorter sequence)


### Step 2: Summarise number of clusters for each cut-off.

summary_header="Cutoff\tClusters\tTotal_sequences\n"
summary_file="${output_dir}/spacer_cluster_summary.tsv"

printf ${summary_header} > ${summary_file}
for logfile in ${log_dir}/all_spacers-clustered-*.txt
do
    cutoff=$(basename -s .txt ${logfile} | sed 's/all_spacers-clustered-//')
    clusters_total=$(tac ${logfile} | grep -m 1 "finished" | awk '{printf ("%s\t%s", $3,$1)}')
    printf "${cutoff}\t${clusters_total}\n" >> ${summary_file}
done
# `tac` reverse reads the file, so grep can find the last occurrence of 'finished'
