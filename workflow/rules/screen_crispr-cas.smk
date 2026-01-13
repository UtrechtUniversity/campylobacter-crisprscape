### Screen for CRISPR-Cas
## Run CCTyper and parse/merge its output


rule crisprcastyper:
    input:
        batch="resources/ATB/assemblies/{batch}/",
    output:
        "results/cctyper/{batch}/complete",
    params:
        out_dir=subpath(output[0], parent=True),
    conda:
        "../envs/cctyper.yaml"
    threads: config["cctyper"]["threads"]
    log:
        "log/cctyper/{batch}.txt",
    benchmark:
        "log/benchmark/cctyper/{batch}.txt"
    shell:
        """
find -L {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 'rm -rf "{params.out_dir}{{/.}}" &&\
 cctyper -t 1 {{}} "{params.out_dir}/{{/.}}"' > {log} 2>&1

touch {output}
        """


rule parse_cctyper:
    input:
        "results/cctyper/{batch}/complete",
    output:
        "results/cctyper/{batch}/parsed",
    conda:
        "../envs/pandas.yaml"
    threads: config["parse_cctyper"]["threads"]
    log:
        "log/parse_cctyper/{batch}.txt",
    benchmark:
        "log/benchmark/parse_cctyper/{batch}.txt"
    shell:
        """
find $(dirname {input}) -mindepth 1 -maxdepth 1 -type d -print0 |\
parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    python workflow/scripts/cctyper_extender.py -d {{.}} > {log} 2>&1

touch {output}
        """


rule extract_sequences:
    input:
        flag="results/cctyper/{batch}/parsed",
        genomes="resources/ATB/assemblies/{batch}",
    output:
        "results/cctyper/{batch}/subseq",
    conda:
        "../envs/seqkit.yaml"
    threads: config["extract_sequences"]["threads"]
    log:
        "log/extract_sequences/{batch}.txt",
    benchmark:
        "log/benchmark/extract_sequences/{batch}.txt"
    shell:
        """
find $(dirname {input.flag}) -mindepth 1 -maxdepth 1 -type d -print0 |\
parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    bash workflow/scripts/extract_crispr-cas_from_fasta.sh {{}} {input.genomes} > {log} 2>&1

touch {output}
        """


rule collect_cctyper:
    input:
        cctyper="results/cctyper/{batch}/complete",
        parser="results/cctyper/{batch}/parsed",
    output:
        crispr_cas="results/cctyper/{batch}/CRISPR_Cas-{batch}.tab",
        crisprs_all="results/cctyper/{batch}/crisprs_all-{batch}.tab",
        crisprs_near_cas="results/cctyper/{batch}/crisprs_near_cas-{batch}.tab",
        crisprs_orphan="results/cctyper/{batch}/crisprs_orphan-{batch}.tab",
        spacers="results/cctyper/{batch}/all_spacers-{batch}.fa",
        cas_putative="results/cctyper/{batch}/cas_operons_putative-{batch}.tab",
        cas="results/cctyper/{batch}/cas_operons-{batch}.tab",
        csv="results/cctyper/{batch}/CRISPR-Cas-{batch}.csv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/cctyper/collect_{batch}.txt",
    benchmark:
        "log/benchmark/cctyper/collect_{batch}.txt"
    shell:
        """
bash workflow/scripts/concatenate_cctyper_output.sh $(dirname {input.cctyper}) > {log} 2>&1
echo "\n========================" >> {log}
bash workflow/scripts/concatenate_cctyper_csv.sh $(dirname {input.parser}) >> {log} 2>&1

find $(dirname {input.cctyper}) -mindepth 3 -maxdepth 3 -name "*.fa" -exec cat {{}} + > {output.spacers} 2>> {log}
        """


rule concatenate_all_spacers:
    input:
        expand("results/cctyper/{batch}/all_spacers-{batch}.fa", batch=BATCHES),
    output:
        "results/cctyper/all_spacers.fa",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/concatenate_all_spacers.txt",
    benchmark:
        "log/benchmark/concatenate_all_spacers.txt"
    shell:
        """
cat {input} > {output} 2> {log}
        """


rule cluster_all_spacers:
    input:
        "results/cctyper/all_spacers.fa",
    output:
        clusters=expand(
            "results/cctyper/all_spacers-clustered-{cutoff}.clstr",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        spacers=expand(
            "results/cctyper/all_spacers-clustered-{cutoff}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        summary="results/cctyper/spacer_cluster_summary.tsv",
    params:
        work_dir=subpath(input[0], parent=True),
        log_dir="log/spacer_clustering",
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_all_spacers.txt",
    benchmark:
        "log/benchmark/cluster_all_spacers.txt"
    shell:
        """
bash workflow/scripts/cluster_all_spacers.sh\
    {input}\
    {params.work_dir}\
    {params.log_dir} > {log} 2>&1
        """


rule cluster_unique_spacers:
    input:
        "results/cctyper/all_spacers.fa",
    output:
        clusters="results/cctyper/all_spacers-clustered.clstr",
        spacers="results/cctyper/all_spacers-clustered",
        distribution="results/cctyper/all_spacers-clustered-distribution.tsv",
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers.txt",
    benchmark:
        "log/benchmark/cluster_unique_spacers.txt"
    shell:
        """
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -sf 1 -d 0 -T {threads}\
 -i {input} -o {output.spacers} > {log} 2>&1

plot_len1.pl {output.clusters}\
 1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999\
 1-10,11-20,21-25,26-30,31-35,36-40,41-50,51-60,61-70,71-999999\
 > {output.distribution}
        """


rule create_crispr_cluster_table:
    input:
        clstr="results/cctyper/all_spacers-clustered.clstr",
        fasta="results/cctyper/all_spacers.fa",
    output:
        "results/all_spacers_table.tsv",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table.txt"
    script:
        "../scripts/make_cluster_table.py"
