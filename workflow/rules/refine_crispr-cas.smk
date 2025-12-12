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
    script:
        "../scripts/merge_identify_batches.sh"


rule merge_cctyper_identify:
    input:
        identify=WORK_DIR + "crispridentify/complete_summary.csv",
        cctyper=expand(
            WORK_DIR + "cctyper/{batch}/crisprs_all-{batch}.tab", batch=BATCHES
        ),
    output:
        table=OUTPUT_DIR + "all_CRISPRS_with_identify.tab",
    threads: 1
    log:
        "log/merge_cctyper_identify",
    script:
        "../scripts/merge_cctyper_crispridentify.sh"


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
