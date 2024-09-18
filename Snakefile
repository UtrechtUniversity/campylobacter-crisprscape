"""
Author: Sam Nooij
Organisation: Utrecht University
Department: Clinical Infectiology (KLIF), Infectious Diseases & Immunology,
  Biomolecular Health Sciences, Faculty of Veterinary Medicine
Date: 2024-09-12

Workflow for testing CRISPR analysis options

Input: Fasta files of Campylobacter whole-genomes
Output: (various)

Example use:
    $ snakemake --profile config

N.B. Variables are set in the configuration files under `config`.
"""

from pathlib import Path

### Step 1: Import configuration file ###

configfile: Path("config/parameters.yaml")

SAMPLES = config["samples"]

INPUT_DIR = config["input_directory"]
OUTPUT_DIR = config["output_directory"]


### Step 2: Specify output files ###

rule all:
    input:
        # CRISPRCasTyper output (generic file that is always created, as not all genomes have CRISPR spacers!)
        expand(OUTPUT_DIR + "cctyper/{sample}/arguments.tab",
               sample = SAMPLES),


### Step 3: Define processing steps that generate the output ###

rule crisprcastyper:
    input:
        INPUT_DIR + "{sample}.fasta"
    output:
        arguments = OUTPUT_DIR + "cctyper/{sample}/arguments.tab",
        putative_operons = OUTPUT_DIR + "cctyper/{sample}/cas_operons_putative.tab",
        genes = OUTPUT_DIR + "cctyper/{sample}/genes.tab",
        hmmer = OUTPUT_DIR + "cctyper/{sample}/hmmer.tab"
    params:
        out_dir = OUTPUT_DIR + "cctyper/{sample}"
    conda:
        "envs/cctyper.yaml"
    threads:
        config["cctyper"]["threads"]
    log:
        "log/cctyper-{sample}.txt"
    benchmark:
        "log/benchmark/cctyper-{sample}.txt"
    shell:
        """
rm -rf {params.out_dir}
cctyper -t {threads} {input} {params.out_dir}
        """
