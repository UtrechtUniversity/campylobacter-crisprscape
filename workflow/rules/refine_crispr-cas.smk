### Refine CRISPR-Cas identifation


rule crispridentify:
    input:
        WORK_DIR + "cctyper/{batch}/subseq",
    output:
        WORK_DIR + "crispridentify/{batch}/complete",
    params:
        out_dir=WORK_DIR + "crispridentify/{batch}",
        arrays=WORK_DIR + "cctyper/{batch}",
    conda:
        "../envs/crispridentify.yml"
    threads: config["crispridentify"]["threads"]
    log:
        "log/crispridentify/{batch}.txt",
    benchmark:
        "log/benchmark/crispridentify/{batch}.txt"
    shell:
        """
cd bin/CRISPRidentify

find ../../{params.arrays}/*/fasta/CRISPR_arrays-with_flanks.fasta -size +0c -print0 |\
parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
'python CRISPRidentify.py --file {{}}\
 --result_folder "../../{params.out_dir}/{{/.}}"\
 --fasta_report True --strand False' > ../../{log} 2>&1

touch ../../{output}
        """


rule merge_crispridentify_batches:
    input:
        expand(WORK_DIR + "crispridentify/{batch}/complete", batch=BATCHES),
    params:
        spacers_crispr=expand(
            WORK_DIR
            + "crispridentify/{batch}/CRISPR_arrays-with_flanks/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
        summary_crispr=expand(
            WORK_DIR
            + "crispridentify/{batch}/CRISPR_arrays-with_flanks/Complete_summary.csv",
            batch=BATCHES,
        ),
    output:
        spacers_crispr=WORK_DIR + "crispridentify/all_spacers.fa",
        summary_crispr=WORK_DIR + "crispridentify/complete_summary.csv",
    threads: 1
    log:
        "log/merge_crispridentify_batches.txt",
    benchmark:
        "log/benchmark/merge_crispridentify_batches.txt"
    shell:
        """
cat {params.spacers_crispr} > {output.spacers_crispr}

for summary in {params.summary_crispr}
do
    header=$(head -n 1 "$summary")
    if [ "$header" == "No arrays found" ]
    then
        continue
    else
        echo $header | tee {output.summary_crispr}
        break
    fi
done

for summary in {params.summary_crispr}
do
    tail -n +2 "$summary" >> {output.summary_crispr}
done
        """


rule merge_cctyper_identify:
    input:
        identify=WORK_DIR + "crispridentify/complete_summary.csv",
        cctyper=expand(
            WORK_DIR + "cctyper/{batch}/crisprs_all-{batch}.tab", batch=BATCHES
        ),
    output:
        OUTPUT_DIR + "all_CRISPRS_with_identify.tab",
    params:
        tmp1="tmp_file1",
        tmp2="tmp_file2",
    threads: 1
    log:
        "log/merge_cctyper_identify",
    shell:
        """
first=True
for summary in {input.cctyper}
do
    if [ $first == True ]
    then
        cat $summary > {params.tmp1}
        first=False
    else
        tail -n +2 $summary >> {params.tmp1}
    fi
done

header=$(head -n 1 {input.identify} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t")
tail -n +2 {input.identify} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t" > {params.tmp2}
first=True
while read line
do
    if [ $first == True ]
    then
        first=False
        echo -e "$line\t$header" > {output}
    else
        sample=$(echo -e "$line" | cut -f 1)
        start_cc=$(echo -e "$line" | cut -f 3)
        start_id=$(expr "$start_cc" + 1)
        match=$(grep "${{sample}}_$start_id" {params.tmp2} || true)
        if [ -z "$match" ]
        then
            echo -e "$line" >> {output}
        else
            while read line2
            do
                if [ "$start_cc" -lt 5000 ];
                then
                    echo -e "$line\t$match" >> {output}
                else
                    start=$(echo -e "$line2" | cut -f 2)
                    start=$(expr "$start" + "$start_cc" - 5000)
                    length=$(echo -e "$line2" | cut -f 4)
                    end=$(expr "$length" + "$start" - 1)
                    begin=$(echo -e "$line2" | cut -f 1)
                    rest=$(echo -e "$line2" | cut -f 4-9)
                    echo -e "$line\t$begin\t$start\t$end\t$rest" >> {output}
                fi
            done <<< "$match"
        fi
    fi
done < {params.tmp1}
rm -f {params.tmp1} {params.tmp2}
        """


rule cluster_unique_spacers_crispridentify:
    input:
        WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        clusters=WORK_DIR + "crispridentify/all_spacers-clustered.clstr",
        spacers=WORK_DIR + "crispridentify/all_spacers-clustered",
        distribution=WORK_DIR + "crispridentify/all_spacers-clustered-distribution.tsv",
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers_identify.txt",
    benchmark:
        "log/benchmark/cluster_unique_spacers_identify.txt"
    shell:
        """
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -sf 1 -d 0 -T {threads}\
 -i {input} -o {output.spacers} > {log} 2>&1

plot_len1.pl {output.clusters}\
 1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999\
 1-10,11-20,21-25,26-30,31-35,36-40,41-50,51-60,61-70,71-999999\
 > {output.distribution}
        """


rule create_crispr_cluster_table_identify:
    input:
        clstr=WORK_DIR + "crispridentify/all_spacers-clustered.clstr",
        fasta=WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        OUTPUT_DIR + "all_spacers_table_identify.tsv",
    conda:
        "../envs/pyfaidx.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table_identify.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table_identify.txt"
    script:
        "../scripts/make_cluster_table_identify.py"
