### Refine CRISPR-Cas identifation


rule split_crispridentify_input:
    input:
        "results/crispr_fasta/{batch}/CRISPR-Cas-with_flanks.fasta",
    output:
        directory("results/crispr_fasta/{batch}/split_with_flanks"),
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    log:
        "log/split_crispridentify_input/{batch}.txt",
    benchmark:
        "log/benchmark/split_crispridentify_input/{batch}.txt"
    shell:
        """
seqkit split2 {input} -s 1 -N -O {output} > {log} 2>&1
        """


rule crispridentify:
    input:
        folder="results/crispr_fasta/{batch}/split_with_flanks",
        fasta="results/crispr_fasta/{batch}/CRISPR-Cas-with_flanks.fasta",
    output:
        spacers="results/crispridentify/{batch}/Complete_spacer_dataset.fasta",
        summary="results/crispridentify/{batch}/Complete_summary.csv",
    params:
        out_dir=subpath(output[0], parent=True),
    conda:
        "../envs/crispridentify.yaml"
    threads: 1
    log:
        "log/crispridentify/{batch}.txt",
    benchmark:
        "log/benchmark/crispridentify/{batch}.txt"
    shell:
        """
rm -f "{input.folder}/.snakemake_timestamp"
cd bin/CRISPRidentify

python CRISPRidentify.py --input_folder "../../{input.folder}"\
 --result_folder "../../{params.out_dir}" --fasta_report True\
 --strand False > "../../{log}" 2>&1
        """


rule merge_crispridentify_batches:
    input:
        spacers_crispr=expand(
            "results/crispridentify/{batch}/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
        summary_crispr=expand(
            "results/crispridentify/{batch}/Complete_summary.csv",
            batch=BATCHES,
        ),
    output:
        spacers_crispr="results/spacers-crispridentify.fa",
        summary_crispr="results/crisprs-crispridentify.csv",
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
        identify="results/crispridentify/complete_summary.csv",
        cctyper=expand("results/cctyper/{batch}/crisprs_all-{batch}.tab", batch=BATCHES),
    output:
        table="results/all_CRISPRs_with_identify.tab",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/merge_cctyper_identify.txt",
    benchmark:
        "log/benchmark/merge_cctyper_identify.txt"
    script:
        "../scripts/merge_cctyper_crispridentify.sh"


rule cluster_all_spacers_crispridentify:
    input:
        rules.merge_crispridentify_batches.output.spacers_crispr,
    output:
        clusters=expand(
            "results/crispridentify/all_spacers-clustered-{cutoff}.clstr",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        spacers=expand(
            "results/crispridentify/all_spacers-clustered-{cutoff}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        summary="results/crispridentify/spacer_cluster_summary.tsv",
    params:
        work_dir=subpath(input[0], parent=True),
        log_dir="log/spacer_clustering_crispridentify",
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_all_spacers_crispridentify.txt",
    benchmark:
        "log/benchmark/cluster_all_spacers_crispridentify.txt"
    shell:
        """
bash workflow/scripts/cluster_all_spacers.sh\
    {input}\
    {params.work_dir}\
    {params.log_dir} > {log} 2>&1
        """


rule cluster_unique_spacers_crispridentify:
    input:
        "results/crispridentify/all_spacers.fa",
    output:
        clusters="results/crispridentify/all_spacers-clustered.clstr",
        spacers="results/crispridentify/all_spacers-clustered",
        distribution="results/crispridentify/all_spacers-clustered-distribution.tsv",
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
        clstr="results/crispridentify/all_spacers-clustered.clstr",
        fasta="results/crispridentify/all_spacers.fa",
    output:
        "results/all_spacers_table_identify.tsv",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table_identify.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table_identify.txt"
    script:
        "../scripts/make_cluster_table_identify.py"
