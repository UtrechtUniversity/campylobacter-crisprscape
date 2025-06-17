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
        # Combined CCTyper output as CSV + BED files
        expand(OUTPUT_DIR + "cctyper/{batch}/parsed", batch=BATCHES),
        # Concatenated CCTyper output
        expand(
            OUTPUT_DIR + "cctyper/{batch}/{filename}-{batch}.tab",
            batch=BATCHES,
            filename=[
                "CRISPR_Cas",
                "crisprs_all",
                "crisprs_near_cas",
                "crisprs_orphan",
                "cas_operons",
            ],
        ),
        expand(OUTPUT_DIR + "cctyper/{batch}/all_spacers-{batch}.fa", batch=BATCHES),
        # CCTyper CRISPR spacer cluster analysis report
        OUTPUT_DIR + "cctyper/spacer_cluster_summary.tsv",
        # Cluster unique CRISPR spacers
        OUTPUT_DIR + "all_spacers-clustered.clstr",
        # Extracted CRISPR arrays (as fasta)
        expand(OUTPUT_DIR + "arrays/{batch}/complete", batch=BATCHES),
        #CRISPRidentify output
        expand(OUTPUT_DIR + "crispridentify/{batch}/complete", batch = BATCHES),
        #concatenated CRISPRidentify output
        OUTPUT_DIR + "crispridentify/complete_summary.csv",
        OUTPUT_DIR + "crispridentify/all_spacers.fa"




### Step 3: Define processing steps that generate the output ###


rule crisprcastyper:
    input:
        batch=INPUT_DIR + "{batch}/",
    output:
        OUTPUT_DIR + "cctyper/{batch}/complete",
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
find {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 'rm -rf "{params.out_dir}{{/.}}" &&\
 cctyper -t 1 {{}} "{params.out_dir}{{/.}}"' > {log} 2>&1

touch {output}
        """


rule parse_cctyper:
    input:
        OUTPUT_DIR + "cctyper/{batch}/complete",
    output:
        OUTPUT_DIR + "cctyper/{batch}/parsed",
    conda:
        "envs/pandas.yaml"
    threads: config["parse_cctyper"]["threads"]
    log:
        "log/parse_cctyper/{batch}.txt",
    benchmark:
        "log/benchmark/parse_cctyper/{batch}.txt"
    shell:
        """
find $(dirname {input}) -mindepth 1 -maxdepth 1 -type d -print0 |\
parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    python bin/cctyper_extender.py -d {{.}} > {log} 2>&1

touch {output}
        """


rule extract_sequences:
    input:
        OUTPUT_DIR + "cctyper/{batch}/parsed",
    output:
        OUTPUT_DIR + "cctyper/{batch}/subseq",
    conda:
        "envs/seqkit.yaml"
    threads: config["extract_sequences"]["threads"]
    log:
        "log/extract_sequences/{batch}.txt"
    benchmark:
        "log/benchmark/extract_sequences/{batch}.txt"
    shell:
        """
find $(dirname {input}) -mindepth 1 -maxdepth 1 -type d -print0 |\
parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    bash bin/extract_crispr-cas_from_fasta.sh {{}} > {log} 2>&1

touch {output}
        """

rule collect_cctyper:
    input:
        cctyper = OUTPUT_DIR + "cctyper/{batch}/complete",
        parser = OUTPUT_DIR + "cctyper/{batch}/parsed"
    output:
        crispr_cas=OUTPUT_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.tab",
        crisprs_all=OUTPUT_DIR + "cctyper/{batch}/crisprs_all-{batch}.tab",
        crisprs_near_cas=OUTPUT_DIR + "cctyper/{batch}/crisprs_near_cas-{batch}.tab",
        crisprs_orphan=OUTPUT_DIR + "cctyper/{batch}/crisprs_orphan-{batch}.tab",
        spacers=OUTPUT_DIR + "cctyper/{batch}/all_spacers-{batch}.fa",
        cas_putative=temp(
            OUTPUT_DIR + "cctyper/{batch}/cas_operons_putative-{batch}.tab"
        ),
        cas=OUTPUT_DIR + "cctyper/{batch}/cas_operons-{batch}.tab",
        csv=OUTPUT_DIR + "cctyper/{batch}/CRISPR-Cas-{batch}.csv",
    threads: 1
    log:
        "log/cctyper/collect_{batch}.txt",
    benchmark:
        "log/benchmark/cctyper/collect_{batch}.txt"
    shell:
        """
bash bin/concatenate_cctyper_output.sh $(dirname {input.cctyper}) > {log} 2>&1
echo "\n========================" >> {log}
bash bin/concatenate_cctyper_csv.sh $(dirname {input.parser}) >> {log} 2>&1

find $(dirname {input}) -mindepth 3 -maxdepth 3 -name "*.fa" -exec cat {{}} + > {output.spacers} 2>> {log}
        """


rule extract_crispr_cas_locations:
    input:
        OUTPUT_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.tab",
    output:
        OUTPUT_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.bed",
    threads: 1
    log:
        "log/extract_crispr_cas_location/{batch}.txt",
    benchmark:
        "log/benchmark/extract_crispr_cas_location/{batch}.txt"
    shell:
        """
python bin/create_CCTyper_bedfile.py -i {input} -o {output} > {log} 2>&1
        """


rule extract_crispr_array:
    input:
        batch=INPUT_DIR + "{batch}/",
        bed=OUTPUT_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.bed",
    output:
        OUTPUT_DIR + "arrays/{batch}/complete",
    params:
        out_dir=OUTPUT_DIR + "arrays/{batch}/",
    conda:
        "envs/seqkit.yaml"
    threads: config["extract_arrays"]["threads"]
    log:
        "log/extract_crispr_array/{batch}.txt",
    benchmark:
        "log/benchmark/extract_crispr_array/{batch}.txt"
        ""
    shell:
        """
cut -f 1 -d '.' {input.bed} | parallel --jobs {threads} --retry-failed\
 --halt='now,fail=1'\
 'if [ -e "{input.batch}/{{}}.fa" ];\
 then seqkit subseq --bed {input.bed} "{input.batch}/{{}}.fa"\
 -o "{params.out_dir}{{}}.fa";\
 fi' > {log} 2>&1

touch {output}
        """


rule concatenate_all_spacers:
    input:
        expand(OUTPUT_DIR + "cctyper/{batch}/all_spacers-{batch}.fa", batch=BATCHES),
    output:
        OUTPUT_DIR + "cctyper/all_spacers.fa",
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
        OUTPUT_DIR + "cctyper/all_spacers.fa",
    output:
        clusters=expand(
            OUTPUT_DIR + "cctyper/all_spacers-clustered-{cutoff}.clstr",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        spacers=expand(
            OUTPUT_DIR + "cctyper/all_spacers-clustered-{cutoff}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        summary=OUTPUT_DIR + "cctyper/spacer_cluster_summary.tsv",
    params:
        output_dir=OUTPUT_DIR + "cctyper/",
        log_dir="log/spacer_clustering/",
    conda:
        "envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_all_spacers.txt",
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
        OUTPUT_DIR + "cctyper/all_spacers.fa",
    output:
        clusters=OUTPUT_DIR + "all_spacers-clustered.clstr",
        spacers=OUTPUT_DIR + "all_spacers-clustered",
        distribution=OUTPUT_DIR + "all_spacers-clustered-distribution.tsv",
    conda:
        "envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers.txt",
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


rule create_crispr_cluster_table:
    input:
        clstr=OUTPUT_DIR + "all_spacers-clustered.clstr",
        fasta=OUTPUT_DIR + "cctyper/all_spacers.fa",
    output:
        OUTPUT_DIR + "cctyper/all_spacers_table.tsv",
    conda:
        "envs/pyfaidx.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table.txt"
    script:
        "bin/make_cluster_table.py"

rule crispridentify:
    input:
        OUTPUT_DIR + "arrays/{batch}/complete"
    output:
        OUTPUT_DIR + "crispridentify/{batch}/complete",
    params:
        out_dir=OUTPUT_DIR + "crispridentify/{batch}",
        spacers=OUTPUT_DIR + "arrays/{batch}"
    conda:
        "envs/crispridentify.yml",
    threads: config["crispridentify"]["threads"]
    log:
        "log/crispridentify/{batch}.txt",
    benchmark:
        "log/benchmark/crispridentify/{batch}.txt",
    shell:
        """
    cd bin/CRISPRidentify
    parallel --jobs {threads} --retry-failed --halt='now,fail=1'\
    python CRISPRidentify.py --file {{}} --result_folder "../../{params.out_dir}/{{/.}}" --fasta_report True > ../../{log} 2>&1 \
    ::: ../../{params.spacers}/*.fa
    cd ../../
    touch {output}
    
        """
rule concatenate_crispridentify_output:
    input:
        OUTPUT_DIR + "crispridentify/{batch}/complete"
    output:
        spacers=OUTPUT_DIR + "crispridentify/{batch}/all_spacers_{batch}.fa",
        summary=OUTPUT_DIR + "crispridentify/{batch}/complete_summary_{batch}.csv",
    params:
        spacers=OUTPUT_DIR + "crispridentify/{batch}/*/Complete_Bona-fide_spacer_dataset.fasta",
        summary=OUTPUT_DIR + "crispridentify/{batch}/*/Complete_summary.csv",
    threads: 1
    log:
        "log/concatenate_crispridentify_output_{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_crispridentify_output_{batch}.txt"
    shell:
        """
    cat {params.spacers} > {output.spacers}
    for summary in {params.summary} ; do header=$(head -n 1 "$summary"); if [ "$header" == "No arrays found" ];
    then
        continue;
    else 
        echo $header > {output.summary};
        break;
    fi
    done
    for summary in {params.summary} ; do tail -n +2 "$summary" >> {output.summary} ; done
        """

rule merge_crispridentify_batches:
    input: 
        spacers=expand(OUTPUT_DIR + "crispridentify/{batch}/all_spacers_{batch}.fa", batch=BATCHES),
        summary=expand(OUTPUT_DIR + "crispridentify/{batch}/complete_summary_{batch}.csv", batch=BATCHES)    
    output:
        spacers=OUTPUT_DIR + "crispridentify/all_spacers.fa",
        summary=OUTPUT_DIR + "crispridentify/complete_summary.csv",
    threads: 1
    log:
        "log/merge_crispridentify_batches.txt"
    benchmark:
        "log/benchmark/merge_crispridentify_batches.txt"
    shell:
        """
    cat {input.spacers} > {output.spacers}
    for summary in {input.summary} ; do header=$(head -n 1 "$summary"); if [ "$header" == "No arrays found" ];
    then
        continue;
    else 
        echo $header > {output.summary};
        break;
    fi
    done
    for summary in {input.summary} ; do tail -n +2 "$summary" >> {output.summary} ; done
        """

rule cluster_unique_spacers_crispridentify:
    input:
        OUTPUT_DIR + "crispridentify/all_spacers.fa",
    output:
        clusters=OUTPUT_DIR + "crispridentify/all_spacers-clustered.clstr",
        spacers=OUTPUT_DIR + "crispridentify/all_spacers-clustered",
        distribution=OUTPUT_DIR + "crispridentify/all_spacers-clustered-distribution.tsv",
    conda:
        "envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers_identify.txt",
    benchmark:
        "log/benchmark/cluster_unique_spacers_identify.txt"
    shell:
        """
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -sf 1 -d 0 -T {threads}\
 -i {input} -o {output.spacers}

plot_len1.pl {output.clusters}\
 1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999\
 1-10,11-20,21-25,26-30,31-35,36-40,41-50,51-60,61-70,71-999999\
 > {output.distribution}
        """

rule create_crispr_cluster_table_identify:
    input:
        clstr=OUTPUT_DIR + "crispridentify/all_spacers-clustered.clstr",
        fasta=OUTPUT_DIR + "crispridentify/all_spacers.fa"
    output:
        OUTPUT_DIR + "crispridentify/all_spacers_table.tsv",
    conda:
        "envs/pyfaidx.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table_identify.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table_identify.txt"
    script:
        "bin/make_cluster_table.py"


rule genomad:
    input:
        fasta=OUTPUT_DIR + "{batch}.fasta",
        db=config["genomad_database"],
    output:
        aggregated_classification=OUTPUT_DIR
        + "genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
        plasmid_summary=OUTPUT_DIR
        + "genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
        virus_summary=OUTPUT_DIR
        + "genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
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
        batch=INPUT_DIR + "{batch}/",
    output:
        OUTPUT_DIR + "jaeger/{batch}/complete",
    params:
        output_dir=OUTPUT_DIR + "jaeger/{batch}/",
    conda:
        "envs/jaeger.yaml"
    threads: config["jaeger"]["threads"]
    log:
        "log/jaeger/{batch}.txt",
    benchmark:
        "log/benchmark/jaeger/{batch}.txt"
    shell:
        """
parallel --jobs {threads} --retry-failed --halt='now,fail=1'\
 Jaeger -p --workers 1 -i {{}} -o "{params.output_dir}{{/.}}" --overwrite\
 > {log} 2>&1 ::: {input.batch}/*.fa

touch {output}
        """
