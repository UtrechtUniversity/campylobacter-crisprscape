"""
Author: Sam Nooij
Organisation: Utrecht University
Department: Clinical Infectiology (KLIF), Infectious Diseases & Immunology,
  Biomolecular Health Sciences, Faculty of Veterinary Medicine
Date: 2024-11-06

Workflow for testing CRISPR analysis options
In contrast to the 'regular' Snakefile workflow, which works
per genome file, this workflow works per batch and runs GNU
parallel to parallelise processing of the genomes within
each batch.


Input: Fasta files of Campylobacter whole-genomes
Output: (various)

Example use:
    $ snakemake --profile config

N.B. Variables are set in the configuration files under `config`.
"""

from pathlib import Path
import functools
import operator

### Step 1: Import configuration file ###


configfile: Path("config/parameters.yaml")


# Use Python functions to automatically detect batches of genomes fasta files
# in the input directory as 'BATCHES'
INPUT_DIR = config["input_directory"]

BATCH_PATHS = list(Path(INPUT_DIR).glob("batch_*"))
for batch in BATCH_PATHS:
    assert Path(batch).is_dir(), f"Batches must be directories, got {batch}"

BATCHES = [batch.stem for batch in BATCH_PATHS]

OUTPUT_DIR = config["output_directory"]


### Step 2: Specify output files ###


rule all:
    input:
        # Concatenated CCTyper output
        expand(OUTPUT_DIR + "cctyper/{batch}/{filename}-{batch}.tab",
               batch = BATCHES,
               filename = [ "CRISPR_Cas", "crisprs_all", "crisprs_near_cas", "crisprs_orphan" ]
        ),
        expand(OUTPUT_DIR + "cctyper/{batch}/all_spacers-{batch}.fa",
               batch = BATCHES),

        # CCTyper CRISPR spacer cluster analysis report
        OUTPUT_DIR + "cctyper/spacer_cluster_summary.tsv",

        # Cluster unique CRISPR spacers
        OUTPUT_DIR + "all_spacers-clustered.clstr",

        # geNomad output
        expand(
            OUTPUT_DIR + "genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
            batch=BATCHES,
        ),

         # Jaeger output
         expand(
             OUTPUT_DIR + "jaeger/{batch}/complete",
             batch=BATCHES,
         ),


### Step 3: Define processing steps that generate the output ###


rule crisprcastyper:
    input:
        batch=INPUT_DIR + "{batch}/",
    output:
        OUTPUT_DIR + "cctyper/{batch}/complete"
    params:
        out_dir=OUTPUT_DIR + "cctyper/{batch}/",
    conda:
        "envs/cctyper.yaml"
    threads: config["cctyper"]["threads"]
    log:
        "log/cctyper/{batch}.txt",
    benchmark:
        "log/benchmark/cctyper/{batch}.txt"
    shell:
        """
parallel --jobs {threads} --retry-failed --halt='now,fail=1'\
 rm -rf "{params.out_dir}{{/.}}" &&\
 cctyper -t 1 {{}} "{params.out_dir}{{/.}}" > {log} 2>&1\
 ::: {input.batch}/*.fa

touch {output}
        """


rule collect_cctyper:
    input:
        OUTPUT_DIR + "cctyper/{batch}/complete"
    output:
        crispr_cas=OUTPUT_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.tab",
        crisprs_all=OUTPUT_DIR + "cctyper/{batch}/crisprs_all-{batch}.tab",
        crisprs_near_cas=OUTPUT_DIR + "cctyper/{batch}/crisprs_near_cas-{batch}.tab",
        crisprs_orphan=OUTPUT_DIR + "cctyper/{batch}/crisprs_orphan-{batch}.tab",
        spacers=OUTPUT_DIR + "cctyper/{batch}/all_spacers-{batch}.fa",
    threads: 1
    log:
        "log/cctyper/collect_{batch}.txt"
    benchmark:
        "log/benchmark/cctyper/collect_{batch}.txt"
    shell:
        """
bash bin/concatenate_cctyper_output.sh $(dirname {input}) > {log} 2>&1

find $(dirname {input}) -mindepth 3 -maxdepth 3 -name "*.fa" -exec cat {{}} + > {output.spacers} 2>> {log}
        """


rule concatenate_all_spacers:
    input:
        expand(OUTPUT_DIR + "cctyper/{batch}/all_spacers-{batch}.fa",
               batch = BATCHES)
    output:
        OUTPUT_DIR + "cctyper/all_spacers.fa"
    threads: 1
    log:
        "log/concatenate_all_spacers.txt"
    benchmark:
        "log/benchmark/concatenate_all_spacers.txt"
    shell:
        """
cat {input} > {output} 2> {log}
        """


rule cluster_all_spacers:
    input:
        OUTPUT_DIR + "cctyper/all_spacers.fa"
    output:
        clusters=expand(OUTPUT_DIR + "cctyper/all_spacers-clustered-{cutoff}.clstr",
                       cutoff = [ 1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81 ]),
        spacers=expand(OUTPUT_DIR + "cctyper/all_spacers-clustered-{cutoff}",
                       cutoff = [ 1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81 ]),
        summary=OUTPUT_DIR + "cctyper/spacer_cluster_summary.tsv",
    params:
        output_dir=OUTPUT_DIR + "cctyper/",
        log_dir="log/spacer_clustering/"
    conda:
        "envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_all_spacers.txt"
    benchmark:
        "log/benchmark/cluster_all_spacers.txt"
    shell:
        """
bash bin/cluster_all_spacers.sh\
    {input}\
    {params.output_dir}\
    {params.log_dir} > {log} 2>&1
        """


rule cluster_unique_spacers:
    input:
        OUTPUT_DIR + "cctyper/all_spacers.fa"
    output:
        clusters=OUTPUT_DIR + "all_spacers-clustered.clstr",
        spacers=OUTPUT_DIR + "all_spacers-clustered",
        distribution=OUTPUT_DIR + "all_spacers-clustered-distribution.tsv",
    conda:
        "envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers.txt"
    benchmark:
        "log/benchmark/cluster_unique_spacers.txt"
    shell:
        """
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -sf 1 -d 0 -T {threads}\
 -i {input} -o {output.spacers}

plot_len1.pl {output.clusters}\
 1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999\
 1-10,11-20,21-25,26-30,31-35,36-40,41-50,51-60,61-70,71-999999\
 > {output.distribution}
        """


rule concatenate_batches:
    input:
        INPUT_DIR + "{batch}"
    output:
        temp(OUTPUT_DIR + "{batch}.fasta")
    threads: 1
    log:
        "log/concatenate_{batch}.txt"
    benchmark:
        "log/benchmark/concatenate_{batch}.txt"
    shell:
        """
cat {input}/*.fa > {output} 2> {log}
        """


rule batched_cctyper:
    input:
        OUTPUT_DIR + "{batch}.fasta"
    output:
        arguments=OUTPUT_DIR + "cctyper/test/{batch}/arguments.tab",
        putative_operons=OUTPUT_DIR + "cctyper/test/{batch}/cas_operons_putative.tab",
        genes=OUTPUT_DIR + "cctyper/test/{batch}/genes.tab",
        hmmer=OUTPUT_DIR + "cctyper/test/{batch}/hmmer.tab",
    params:
        out_dir=OUTPUT_DIR + "cctyper/test/{batch}",
    conda:
        "envs/cctyper.yaml"
    threads: config["cctyper"]["threads"]
    log:
        "log/cctyper/{batch}.txt"
    benchmark:
        "log/benchmark/cctyper/{batch}.txt"
    shell:
        """
rm -rf {params.out_dir}
cctyper -t {threads} {input} {params.out_dir} > {log} 2>&1
        """


rule genomad:
    input:
        fasta=OUTPUT_DIR + "{batch}.fasta",
        db=config["genomad_database"],
    output:
        aggregated_classification=OUTPUT_DIR + "genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
        plasmid_summary=OUTPUT_DIR + "genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
        virus_summary=OUTPUT_DIR + "genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
    params:
        output_dir=OUTPUT_DIR + "genomad/{batch}/",
    conda:
        "envs/genomad.yaml"
    threads: config["genomad"]["threads"]
    log:
        "log/genomad/{batch}.txt",
    benchmark:
        "log/benchmark/genomad/{batch}.txt"
    shell:
        """
genomad end-to-end -t {threads} --cleanup --enable-score-calibration\
 {input.fasta} {params.output_dir} {input.db} > {log} 2>&1
        """


rule jaeger:
    input:
        batch=INPUT_DIR + "{batch}/"
    output:
        OUTPUT_DIR + "jaeger/{batch}/complete"
    params:
        output_dir=OUTPUT_DIR + "jaeger/{batch}/"
    conda:
        "envs/jaeger.yaml"
    threads: config["jaeger"]["threads"]
    log:
        "log/jaeger/{batch}.txt"
    benchmark:
        "log/benchmark/jaeger/{batch}.txt"
    shell:
        """
parallel --jobs {threads} --retry-failed --halt='now,fail=1'\
 Jaeger -p --workers 1 -i {{}} -o "{params.output_dir}{{/.}}" --overwrite\
 > {log} 2>&1 ::: {input.batch}/*.fa

touch {output}
        """
