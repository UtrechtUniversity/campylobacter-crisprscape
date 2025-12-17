### Screen for CRISPR-Cas
## Run CCTyper and parse/merge its output


rule crisprcastyper:
    input:
        batch="data/tmp/assemblies/{batch}/",
    output:
        "data/tmp/cctyper/{batch}/complete",
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
 cctyper -t 1 {{}} "{params.out_dir}{{/.}}"' > {log} 2>&1

touch {output}
        """


rule parse_cctyper:
    input:
        "data/tmp/cctyper/{batch}/complete",
    output:
        "data/tmp/cctyper/{batch}/parsed",
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
        "data/tmp/cctyper/{batch}/parsed",
    output:
        "data/tmp/cctyper/{batch}/subseq",
    conda:
        "../envs/seqkit.yaml"
    threads: config["extract_sequences"]["threads"]
    log:
        "log/extract_sequences/{batch}.txt",
    benchmark:
        "log/benchmark/extract_sequences/{batch}.txt"
    shell:
        """
find $(dirname {input}) -mindepth 1 -maxdepth 1 -type d -print0 |\
parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    bash workflow/scripts/extract_crispr-cas_from_fasta.sh {{}} > {log} 2>&1

touch {output}
        """


rule concatenate_all_spacers:
    input:
        expand("data/tmp/cctyper/{batch}/all_spacers-{batch}.fa", batch=BATCHES),
    output:
        "data/tmp/cctyper/all_spacers.fa",
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
        "data/tmp/cctyper/all_spacers.fa",
    output:
        clusters=expand(
            "data/tmp/cctyper/all_spacers-clustered-{cutoff}.clstr",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        spacers=expand(
            "data/tmp/cctyper/all_spacers-clustered-{cutoff}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        summary="data/tmp/cctyper/spacer_cluster_summary.tsv",
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
        "data/tmp/cctyper/all_spacers.fa",
    output:
        clusters="data/tmp/cctyper/all_spacers-clustered.clstr",
        spacers="data/tmp/cctyper/all_spacers-clustered",
        distribution="data/tmp/cctyper/all_spacers-clustered-distribution.tsv",
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
        clstr="data/tmp/cctyper/all_spacers-clustered.clstr",
        fasta="data/tmp/cctyper/all_spacers.fa",
    output:
        "data/processed/all_spacers_table.tsv",
    conda:
        "../envs/pyfaidx.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table.txt"
    script:
        "../scripts/make_cluster_table.py"
