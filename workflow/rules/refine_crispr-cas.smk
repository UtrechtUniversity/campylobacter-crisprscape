### Refine CRISPR-Cas identifation


rule split_crispridentify_input:
    input:
        crisprcas="results/crispr_fasta/{batch}/CRISPR-Cas-with_flanks.fasta",
        orphancrispr="results/crispr_fasta/{batch}/Orphan_CRISPRs-with_flanks.fasta"
    output:
        directory("results/crispr_fasta/{batch}/split_with_flanks"),
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/split_crispridentify_input/{batch}.txt",
    benchmark:
        "log/benchmark/split_crispridentify_input/{batch}.txt"
    shell:
        r"""
for fasta in {input}
do
    seqkit split2 ${{fasta}} -s 1 -N -O {output} > {log} 2>&1
done
        """


rule crispridentify:
    input:
        folder="results/crispr_fasta/{batch}/split_with_flanks",
    output:
        spacers="results/crispridentify/{batch}/Complete_spacer_dataset.fasta",
        summary="results/crispridentify/{batch}/Complete_summary.csv",
    params:
        out_dir=subpath(output[0], parent=True),
    conda:
        "../envs/crispridentify.yaml"
    threads: config["crispridentify"]["threads"]
    resources:
        mem_mb=int(config["crispridentify"]["memory"]),
        walltime=int(config["crispridentify"]["time"]),
        runtime=int(config["crispridentify"]["time"]),
    log:
        "log/crispridentify/{batch}.txt",
    benchmark:
        "log/benchmark/crispridentify/{batch}.txt"
    shell:
        r"""
rm -f "{input.folder}/.snakemake_timestamp"
rm -rf "{params.out_dir}"
cd bin/CRISPRidentify

python CRISPRidentify.py --input_folder "../../{input.folder}"\
 --result_folder "../../{params.out_dir}" --fasta_report True\
 --strand False --min_repeats 2 > "../../{log}" 2>&1
        """


rule create_crispridentify_crispr_table:
    input:
        table=expand(
            "results/crispridentify/{batch}/Complete_summary.csv",
            batch=BATCHES,
        ),
    output:
        "results/crisprs-final.tsv",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/create_crispridentify_crispr_table.txt",
    benchmark:
        "log/benchmark/create_crispridentify_crispr_table.txt"
    script:
        "../scripts/summarise_crisprs.R"


rule concatenate_crispridentify_spacers:
    input:
        expand(
            "results/crispridentify/{batch}/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
    output:
        "results/spacers-final.fasta",
    conda:
        "../envs/bash.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/concatenate_crispridentify_spacers.txt",
    benchmark:
        "log/benchmark/concatenate_crispridentify_spacers.txt"
    shell:
        r"""
cat {input} > {output} 2> {log}
        """


rule prepare_crispridentify_spacer_table:
    input:
        rules.concatenate_crispridentify_spacers.output[0],
    output:
        temp("results/spacers-final_prep.tsv"),
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/prepare_crispridentify_spacer_table.txt",
    benchmark:
        "log/benchmark/prepare_crispridentify_spacer_table.txt"
    shell:
        r"""
seqkit fx2tab -l -HQ {input} > {output} 2> {log}
        """


rule deduplicate_crispridentify_spacers:
    input:
        rules.concatenate_crispridentify_spacers.output[0],
    output:
        "results/spacers-final-deduplicated_and_counted.fasta",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/deduplicate_crispridentify_spacers.txt",
    benchmark:
        "log/benchmark/deduplicate_crispridentify_spacers.txt"
    shell:
        r"""
seqkit seq {input} -s -w 0 |\
 sort | uniq -c | sort -rn |\
 while read count seq;\
 do printf ">${{count}}.${{seq}}\n${{seq}}\n";\
 done > {output} 2> {log}
        """


rule cluster_crispridentify_spacers:
    input:
        rules.concatenate_crispridentify_spacers.output[0],
    output:
        expand(
            "results/cluster-crispridentify/spacers-clustered-{cutoff}{ext}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
            ext=[".clstr", "", ".log"],
        ),
        summary="results/cluster-crispridentify/spacer_cluster_summary.tsv",
    params:
        work_dir=subpath(output[0], parent=True),
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/cluster_crispridentify_spacers.txt",
    benchmark:
        "log/benchmark/cluster_crispridentify_spacers.txt"
    shell:
        r"""
bash workflow/scripts/cluster_all_spacers.sh\
 {input}\
 {params.work_dir} > {log} 2>&1
        """


rule cluster_unique_spacers_crispridentify:
    input:
        rules.deduplicate_crispridentify_spacers.output[0],
    output:
        clusters="results/cluster-crispridentify/spacers-clustered.clstr",
        spacers="results/cluster-crispridentify/spacers-clustered",
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/cluster_unique_spacers_crispridentify.txt",
    benchmark:
        "log/benchmark/cluster_unique_spacers_crispridentify.txt"
    shell:
        r"""
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -d 0 -T {threads}\
 -i {input} -o {output.spacers} > {log} 2>&1
        """


rule create_spacer_table_crispridentify:
    input:
        clstr=rules.cluster_unique_spacers_crispridentify.output.clusters,
        spacer=rules.prepare_crispridentify_spacer_table.output[0],
    output:
        spacer="results/spacers-final.tsv",
        cluster="results/spacer_clusters-final.tsv",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/create_spacer_table_crispridentify.txt",
    benchmark:
        "log/benchmark/create_spacer_table_crispridentify.txt"
    script:
        "../scripts/make_cluster_table.py"


rule prepare_spacer_cluster_fasta_crispridentify:
    input:
        rules.create_spacer_table_crispridentify.output.cluster,
    output:
        common="results/spacers/most_common-final.fasta",
        short="results/spacers/shortest-final.fasta",
        long="results/spacers/longest-final.fasta",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/prepare_spacer_cluster_fasta_crispridentify.txt",
    benchmark:
        "log/benchmark/prepare_spacer_cluster_fasta_crispridentify.txt"
    shell:
        """
grep -v "Cluster" {input} | cut -f 1,3 | seqkit tab2fx -o {output.common} > {log} 2>&1
grep -v "Cluster" {input} | cut -f 1,6 | seqkit tab2fx -o {output.short} > {log} 2>&1
grep -v "Cluster" {input} | cut -f 1,8 | seqkit tab2fx -o {output.long} > {log} 2>&1
        """


rule convert_spacer_formats_crispridentify:
    input:
        "results/spacers-final.tsv",
    output:
        table_numbers="results/arrays/array_IDs-final.txt",
        table_sequences="results/arrays/array_seqs-final.txt",
        fasta_numbers="results/arrays/array_IDs-final.fasta",
        fasta_sequences="results/arrays/array_seqs-final.fasta",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/convert_spacer_formats_crispridentify.txt",
    benchmark:
        "log/benchmark/convert_spacer_formats_crispridentify.txt"
    script:
        "../scripts/convert_spacer_arrays.py"
