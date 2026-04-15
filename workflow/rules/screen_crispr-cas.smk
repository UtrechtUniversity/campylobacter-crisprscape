### Screen for CRISPR-Cas
## Run CCTyper and parse/merge its output


rule crisprcastyper:
    input:
        "resources/ATB/assemblies-concatenated/{batch}.fasta",
    output:
        crispr_cas="results/cctyper/{batch}/CRISPR_Cas.tab",
        cas="results/cctyper/{batch}/cas_operons_putative.tab",
        hmmer="results/cctyper/{batch}/hmmer.tab",
        crispr="results/cctyper/{batch}/crisprs_all.tab",
        orphan="results/cctyper/{batch}/crisprs_orphan.tab",
        spacers=directory("results/cctyper/{batch}/spacers/"),
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
        r"""
rm -r {params.out_dir} &&\
 cctyper -t {threads} {input} --prodigal meta {params.out_dir}\
 --minRL 21 --maxRL 55 --minSL 18 --maxSL 78\
 --minNR 1 --simplelog > {log} 2>&1
        """


rule parse_cctyper:
    input:
        crispr_cas="results/cctyper/{batch}/CRISPR_Cas.tab",
        cas="results/cctyper/{batch}/cas_operons_putative.tab",
        hmmer="results/cctyper/{batch}/hmmer.tab",
        crispr="results/cctyper/{batch}/crisprs_all.tab",
        orphan="results/cctyper/{batch}/crisprs_orphan.tab",
    output:
        table="results/cctyper-parse/{batch}/CRISPR-Cas.tsv",
        locus_bed="results/cctyper-parse/{batch}/CRISPR-Cas.bed",
        array_bed="results/cctyper-parse/{batch}/CRISPR_arrays.bed",
        operon_bed="results/cctyper-parse/{batch}/Cas_operons.bed",
    params:
        input_dir=subpath(input.crispr_cas, parent=True),
    conda:
        "../envs/pandas.yaml"
    threads: 1
    log:
        "log/parse_cctyper/{batch}.txt",
    benchmark:
        "log/benchmark/parse_cctyper/{batch}.txt"
    script:
        "../scripts/cctyper_extender.py"


rule extract_sequences:
    input:
        bed="results/cctyper-parse/{batch}/CRISPR-Cas.bed",
        genomes="resources/ATB/assemblies-concatenated/{batch}.fasta",
    output:
        bed_dir=directory("results/crispr_fasta/{batch}"),
        crispr_cas="results/crispr_fasta/{batch}/CRISPR-Cas.fasta",
        flanks="results/crispr_fasta/{batch}/CRISPR-Cas-with_flanks.fasta",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    log:
        "log/extract_sequences/{batch}.txt",
    benchmark:
        "log/benchmark/extract_sequences/{batch}.txt"
    script:
        "../scripts/extract_crispr-cas_from_fasta.sh"


rule create_cctyper_crispr_table:
    input:
        table=expand("results/cctyper-parse/{batch}/CRISPR-Cas.tsv", batch=BATCHES),
    output:
        "results/crisprs-cctyper.tsv",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/cctyper/create_cctyper_crispr_table.txt",
    benchmark:
        "log/benchmark/cctyper/create_cctyper_crispr_table.txt"
    script:
        "../scripts/summarise_crisprs.R"


rule concatenate_cctyper_spacers:
    input:
        expand("results/cctyper/{batch}/spacers", batch=BATCHES),
    output:
        "results/spacers-cctyper.fasta",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/concatenate_cctyper_spacers.txt",
    benchmark:
        "log/benchmark/concatenate_cctyper_spacers.txt"
    shell:
        r"""
find {input} -name "*.fa" -exec cat {{}} \; > {output} 2> {log}
        """


rule cluster_cctyper_spacers:
    input:
        rules.concatenate_cctyper_spacers.output[0],
    output:
        expand(
            "results/cluster-cctyper/spacers-clustered-{cutoff}{ext}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
            ext=[".clstr", "", ".log"],
        ),
        summary="results/cluster-cctyper/spacer_cluster_summary.tsv",
    params:
        work_dir=subpath(output[0], parent=True),
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_cctyper_spacers.txt",
    benchmark:
        "log/benchmark/cluster_cctyper_spacers.txt"
    shell:
        r"""
bash workflow/scripts/cluster_all_spacers.sh\
 {input}\
 {params.work_dir} > {log} 2>&1
        """


rule cluster_unique_spacers:
    input:
        rules.concatenate_cctyper_spacers.output[0],
    output:
        clusters="results/cluster-cctyper/spacers-clustered.clstr",
        spacers="results/cluster-cctyper/spacers-clustered",
        distribution="results/cluster-cctyper/spacers-clustered-distribution.tsv",
    conda:
        "../envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers.txt",
    benchmark:
        "log/benchmark/cluster_unique_spacers.txt"
    shell:
        r"""
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -sf 1 -d 0 -T {threads}\
 -i {input} -o {output.spacers} > {log} 2>&1

plot_len1.pl {output.clusters}\
 1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999\
 1-10,11-20,21-25,26-30,31-35,36-40,41-50,51-60,61-70,71-999999\
 > {output.distribution}
        """


rule create_spacer_table_cctyper:
    input:
        clstr="results/cluster-cctyper/spacers-clustered.clstr",
        fasta=rules.concatenate_cctyper_spacers.output[0],
    output:
        "results/spacers-cctyper.tsv",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    log:
        "log/create_spacer_table_cctyper.txt",
    benchmark:
        "log/benchmark/create_spacer_table_cctyper.txt"
    script:
        "../scripts/make_cluster_table.py"
