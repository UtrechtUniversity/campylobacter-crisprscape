### Refine CRISPR-Cas identifation


rule crispridentify:
    input:
        "data/tmp/cctyper/{batch}/subseq",
    output:
        "data/tmp/crispridentify/{batch}/complete",
    params:
        out_dir=subpath(output[0], parent=True),
        arrays=subpath(input[0], parent=True),
    conda:
        "../envs/crispridentify.yaml"
    threads: config["crispridentify"]["threads"]
    log:
        "log/crispridentify/{batch}.txt",
    benchmark:
        "log/benchmark/crispridentify/{batch}.txt"
    shell:
        """
cd bin/CRISPRidentify

find ../../{params.arrays}/*/fasta/CRISPR_arrays-with_flanks.fasta -size +0c -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail,1'\
 'python CRISPRidentify.py --file {{}}\
 --result_folder "../../{params.out_dir}/{{/.}}"\
 --fasta_report True --strand False' > ../../{log} 2>&1

touch ../../{output}
        """


rule merge_crispridentify_batches:
    input:
        flag=expand("data/tmp/crispridentify/{batch}/complete", batch=BATCHES),
        spacers_crispr=expand(
            "data/tmp/crispridentify/{batch}/CRISPR_arrays-with_flanks/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
        summary_crispr=expand(
            "data/tmp/crispridentify/{batch}/CRISPR_arrays-with_flanks/Complete_summary.csv",
            batch=BATCHES,
        ),
    output:
        spacers_crispr="data/tmp/crispridentify/all_spacers.fa",
        summary_crispr="data/tmp/crispridentify/complete_summary.csv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/merge_crispridentify_batches.txt",
    benchmark:
        "log/benchmark/merge_crispridentify_batches.txt"
    script:
        "../scripts/merge_crispridentify_batches.sh"


rule merge_cctyper_identify:
    input:
        identify="data/tmp/crispridentify/complete_summary.csv",
        cctyper=expand(
            "data/tmp/cctyper/{batch}/crisprs_all-{batch}.tab", batch=BATCHES
        ),
    output:
        table="data/processed/all_CRISPRS_with_identify.tab",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/merge_cctyper_identify.txt",
    benchmark:
        "log/benchmark/merge_cctyper_identify.txt"
    script:
        "../scripts/merge_cctyper_crispridentify.sh"


rule cluster_unique_spacers_crispridentify:
    input:
        "data/tmp/crispridentify/all_spacers.fa",
    output:
        clusters="data/tmp/crispridentify/all_spacers-clustered.clstr",
        spacers="data/tmp/crispridentify/all_spacers-clustered",
        distribution="data/tmp/crispridentify/all_spacers-clustered-distribution.tsv",
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
        clstr="data/tmp/crispridentify/all_spacers-clustered.clstr",
        fasta="data/tmp/crispridentify/all_spacers.fa",
    output:
        "data/processed/all_spacers_table_identify.tsv",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table_identify.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table_identify.txt"
    script:
        "../scripts/make_cluster_table_identify.py"
