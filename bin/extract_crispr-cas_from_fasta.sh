#!/usr/bin/env bash

cctyper_dir=$1

sample_name=$(basename ${cctyper_dir})
batch_name=$(basename $(dirname ${cctyper_dir}))
fasta_file="data/tmp/ATB/${batch_name}/${sample_name}.fa"

bed_files=( $(find $1 -mindepth 1 -maxdepth 1 -name "*bed") )


if [ ${#bed_files[@]} > 0 ]
then
    out_dir="${cctyper_dir}/fasta"
    mkdir -p ${out_dir}
    for bedfile in ${bed_files[*]}
    do
        locus_type=$(basename -s .bed ${bedfile})
        echo "Extracting ${locus_type} from ${fasta_file}..."
        seqkit subseq --bed ${bedfile} -U ${fasta_file} -o ${out_dir}/${locus_type}.fasta

        seqkit subseq --bed ${bedfile} -U ${fasta_file}\
         -u 5000 -d 5000 -o ${out_dir}/${locus_type}-with_flanks.fasta

        seqkit subseq --bed ${bedfile} -U ${fasta_file}\
         -u 5000 --only-flank\
         -o ${out_dir}/${locus_type}-upstream.fasta

        seqkit subseq --bed ${bedfile} -U ${fasta_file}\
         -d 5000 --only-flank\
         -o ${out_dir}/${locus_type}-downstream.fasta

         ls -lh ${out_dir}/${locus_type}*.fasta
        echo
    done
else
    echo "Sample ${sample_name} has no CRISPR nor cas"
fi