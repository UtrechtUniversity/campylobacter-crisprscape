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

TYPES = ["CRISPR_arrays", "Cas_operons", "CRISPR-Cas"]

OUTPUT_DIR = config["output_directory"]


### Step 2: Specify output files ###


rule all:
    input:
        # Multilocus Sequence Types (ST) for Campylobacter
        expand(OUTPUT_DIR + "mlst/{batch}-concatenated.tsv", batch=BATCHES),
        # Virus and plasmid predictions per contig
        "data/processed/genomad_predictions.csv",
        "data/processed/jaeger_predictions.csv",
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
        # Combined CCTyper output as CSV + BED files
        expand(OUTPUT_DIR + "cctyper/{batch}/parsed", batch=BATCHES),
        # CCTyper CRISPR spacer cluster analysis report
        OUTPUT_DIR + "cctyper/spacer_cluster_summary.tsv",
        # Cluster unique CRISPR spacers
        OUTPUT_DIR + "all_spacers-clustered.clstr",
        # Extracted CRISPR arrays (as fasta)
        expand(OUTPUT_DIR + "arrays/{batch}/complete", batch=BATCHES),
        #CRISPRidentify output
        expand(OUTPUT_DIR + "crispridentify/{batch}/complete", batch=BATCHES),
        #concatenated CRISPRidentify output
        expand(OUTPUT_DIR + "crispridentify/complete_{type}_summary.csv", type=TYPES),
        expand(OUTPUT_DIR + "crispridentify/all_{type}_spacers.fa", type=TYPES),
        #spacepharer output
        expand(
            OUTPUT_DIR + "spacepharer/predicted_{type}_phage_matches.tsv", type=TYPES
        ),
        expand(
            OUTPUT_DIR + "spacepharer/predicted_{type}_plasmid_matches.tsv", type=TYPES
        ),
        #KMA output
        expand(OUTPUT_DIR + "kma/output/{type}.aln", type=TYPES),


### Step 3: Define processing steps that generate the output ###


rule download_mlst_database:
    output:
        OUTPUT_DIR + "mlst/campylobacter.db",
    params:
        species=config["mlst"]["species"],
    conda:
        "envs/pymlst.yaml"
    threads: 1
    log:
        "log/download_mlst_database.txt",
    benchmark:
        "log/benchmark/download_mlst_database.txt"
    shell:
        """
claMLST import -r pubmlst --no-prompt {output} {params.species} > {log} 2>&1
        """


rule mlst:
    input:
        batch=INPUT_DIR + "{batch}/",
        db=OUTPUT_DIR + "mlst/campylobacter.db",
    output:
        OUTPUT_DIR + "mlst/{batch}/complete",
    conda:
        "envs/pymlst.yaml"
    threads: config["mlst"]["threads"]
    log:
        "log/mlst/{batch}.txt",
    benchmark:
        "log/benchmark/mlst/{batch}.txt"
    shell:
        """
find {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 claMLST search {input.db} {{}} -o "$(dirname {output})/{{/.}}.txt" > {log} 2>&1

touch {output}
        """


rule concatenate_mlst:
    input:
        OUTPUT_DIR + "mlst/{batch}/complete",
    output:
        OUTPUT_DIR + "mlst/{batch}-concatenated.tsv",
    threads: config["mlst"]["threads"]
    log:
        "log/concatenate_mlst/{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_mlst/{batch}.txt"
    shell:
        """
echo -e "Genome\tST" > {output}
find $(dirname {input}) -mindepth 1 -maxdepth 1 -type f -name "*.txt" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 'tail -n 1 {{}} | cut -f 1-2 >> {output}'
        """


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


rule download_padloc_database:
    output:
        OUTPUT_DIR + "padloc/database",
    conda:
        "envs/padloc.yaml"
    threads: 1
    log:
        "log/download_padloc_database.txt",
    benchmark:
        "log/benchmark/download_padloc_database.txt"
    shell:
        """
padloc --db-install v2.0.0
touch {output}
        """


rule padloc:
    input:
        batch=INPUT_DIR + "{batch}/",
        db=OUTPUT_DIR + "padloc/database",
    output:
        OUTPUT_DIR + "padloc/{batch}/complete",
    conda:
        "envs/padloc.yaml"
    threads: config["padloc"]["threads"]
    log:
        "log/padloc/{batch}.txt",
    benchmark:
        "log/benchmark/padloc/{batch}.txt"
    shell:
        """
find {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 'mkdir -p "$(dirname {output})/{{/.}}" && padloc --cpu 1 --fna {{}} --outdir "$(dirname {output})/{{/.}}"' > {log} 2>&1

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
        "log/extract_sequences/{batch}.txt",
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
        cctyper=OUTPUT_DIR + "cctyper/{batch}/complete",
        parser=OUTPUT_DIR + "cctyper/{batch}/parsed",
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

find $(dirname {input.cctyper}) -mindepth 3 -maxdepth 3 -name "*.fa" -exec cat {{}} + > {output.spacers} 2>> {log}
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
        OUTPUT_DIR + "cctyper/{batch}/subseq",
    output:
        OUTPUT_DIR + "crispridentify/{batch}/complete",
    params:
        out_dir=OUTPUT_DIR + "crispridentify/{batch}",
        arrays=OUTPUT_DIR + "cctyper/{batch}",
    conda:
        "envs/crispridentify.yml"
    threads: config["crispridentify"]["threads"]
    log:
        "log/crispridentify/{batch}.txt",
    benchmark:
        "log/benchmark/crispridentify/{batch}.txt"
    shell:
        """
    
    cd bin/CRISPRidentify
    find ../../{params.arrays}/*/fasta/CRISPR-Cas-with_flanks.fasta -size +0c -print0 | \
    parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    python CRISPRidentify.py --file {{}} --result_folder "../../{params.out_dir}/{{/.}}" --fasta_report True --strand False > ../../{log} 2>&1   
    
    find ../../{params.arrays}/*/fasta/Cas_operons-with_flanks.fasta -size +0c -print0 | \
    parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    python CRISPRidentify.py --file {{}} --result_folder "../../{params.out_dir}/{{/.}}" --fasta_report True --strand False >> ../../{log} 2>&1 \

    find ../../{params.arrays}/*/fasta/CRISPR_arrays-with_flanks.fasta -size +0c -print0 | \
    parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    python CRISPRidentify.py --file {{}} --result_folder "../../{params.out_dir}/{{/.}}" --fasta_report True --strand False >> ../../{log} 2>&1 \

    touch ../../{output}
    
        """


rule merge_crispridentify_batches:
    input:
        expand(OUTPUT_DIR + "crispridentify/{batch}/complete", batch=BATCHES),
    params:
        spacers_array=expand(
            OUTPUT_DIR
            + "crispridentify/{batch}/CRISPR_arrays-with_flanks/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
        summary_array=expand(
            OUTPUT_DIR
            + "crispridentify/{batch}/CRISPR_arrays-with_flanks/Complete_summary.csv",
            batch=BATCHES,
        ),
        spacers_cas=expand(
            OUTPUT_DIR
            + "crispridentify/{batch}/Cas_operons-with_flanks/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
        summary_cas=expand(
            OUTPUT_DIR
            + "crispridentify/{batch}/Cas_operons-with_flanks/Complete_summary.csv",
            batch=BATCHES,
        ),
        spacers_crispr=expand(
            OUTPUT_DIR
            + "crispridentify/{batch}/CRISPR-Cas-with_flanks/Complete_spacer_dataset.fasta",
            batch=BATCHES,
        ),
        summary_crispr=expand(
            OUTPUT_DIR
            + "crispridentify/{batch}/CRISPR-Cas-with_flanks/Complete_summary.csv",
            batch=BATCHES,
        ),
    output:
        spacers_array=OUTPUT_DIR + "crispridentify/all_CRISPR_arrays_spacers.fa",
        summary_array=OUTPUT_DIR + "crispridentify/complete_CRISPR_arrays_summary.csv",
        spacers_cas=OUTPUT_DIR + "crispridentify/all_Cas_operons_spacers.fa",
        summary_cas=OUTPUT_DIR + "crispridentify/complete_Cas_operons_summary.csv",
        spacers_crispr=OUTPUT_DIR + "crispridentify/all_CRISPR-Cas_spacers.fa",
        summary_crispr=OUTPUT_DIR + "crispridentify/complete_CRISPR-Cas_summary.csv",
    threads: 1
    log:
        "log/merge_crispridentify_batches.txt",
    benchmark:
        "log/benchmark/merge_crispridentify_batches.txt"
    shell:
        """
    cat {params.spacers_array} > {output.spacers_array}
    cat {params.spacers_cas} > {output.spacers_cas}
    cat {params.spacers_crispr} > {output.spacers_crispr}
    for summary in {params.summary_array} ; do header=$(head -n 1 "$summary"); if [ "$header" == "No arrays found" ];
    then
        continue;
    else 
        echo $header | tee {output.summary_array} {output.summary_cas} {output.summary_crispr};
        break;
    fi
    done
    for summary in {params.summary_array} ; do tail -n +2 "$summary" >> {output.summary_array} ; done
    for summary in {params.summary_cas} ; do tail -n +2 "$summary" >> {output.summary_cas} ; done
    for summary in {params.summary_crispr} ; do tail -n +2 "$summary" >> {output.summary_crispr} ; done
        """


rule cluster_unique_spacers_crispridentify:
    input:
        OUTPUT_DIR + "crispridentify/all_{type}_spacers.fa",
    output:
        clusters=OUTPUT_DIR + "crispridentify/all_{type}_spacers-clustered.clstr",
        spacers=OUTPUT_DIR + "crispridentify/all_{type}_spacers-clustered",
        distribution=OUTPUT_DIR
        + "crispridentify/all_{type}_spacers-clustered-distribution.tsv",
    conda:
        "envs/cdhit.yaml"
    threads: 1
    log:
        "log/cluster_unique_spacers_{type}_identify.txt",
    benchmark:
        "log/benchmark/cluster_unique_spacers_{type}_identify.txt"
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


rule collect_genomad_predictions:
    input:
        aggregated_classification=expand(OUTPUT_DIR
        + "genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
        batch = BATCHES),
        plasmid_summary=expand(OUTPUT_DIR
        + "genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
        batch = BATCHES),
        virus_summary=expand(OUTPUT_DIR
        + "genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
        batch = BATCHES),
    output:
        "data/processed/genomad_predictions.csv"
    conda:
        "envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_genomad_predictions.txt"
    benchmark:
        "log/benchmark/collect_genomad_predictions.txt"
    script: "bin/collect_genomad_predictions.R"


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


rule collect_jaeger_batch:
    input:
        OUTPUT_DIR + "jaeger/{batch}/complete"
    output:
        OUTPUT_DIR + "jaeger/{batch}/jaeger-{batch}.csv"
    params:
        batch="{batch}"
    conda:
        "envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_jaeger_{batch}.txt"
    benchmark:
        "log/benchmark/collect_jaeger_{batch}.txt"
    script: "bin/collect_jaeger_batch.R"


rule collect_jaeger_predictions:
    input:
        expand(OUTPUT_DIR + "jaeger/{batch}/jaeger-{batch}.csv",
        batch = BATCHES),
    output:
        "data/processed/jaeger_predictions.csv"
    threads: 1
    log:
        "log/collect_jaeger_predictions.txt"
    benchmark:
        "log/benchmark/collect_jaeger_predictions.txt"
    script: "bin/collect_jaeger_predictions.sh"


rule spacepharer_spacer_setup:
    input:
        spacers=OUTPUT_DIR + "crispridentify/all_{type}_spacers.fa",
    output:
        spacer_DB=OUTPUT_DIR + "spacepharer/DB_{type}/querysetDB",
    params:
        tmp_folder=OUTPUT_DIR + "spacepharer/tmpFolder",
    threads: 48
    log:
        "log/spacepharer/spacepharer_spacer_{type}_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_spacer_{type}_setup.txt"
    shell:
        """
        spacer_DB=$(dirname {output.spacer_DB})
        rm -rf $spacer_DB/* > {log} 2>&1 
        spacepharer createsetdb {input.spacers} {output.spacer_DB} {params.tmp_folder} --extractorf-spacer 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_phage_setup:
    input:
        DB=config["spacepharer_phage_database"],
    output:
        phage_DB=OUTPUT_DIR + "spacepharer/phage_DB/targetsetDB",
        phage_control_DB=OUTPUT_DIR + "spacepharer/phage_DB/controlsetDB",
    params:
        tmp_folder=OUTPUT_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_phage_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_setup.txt"
    shell:
        """
        phage_DB=$(dirname {output.phage_DB})
        rm -rf $phage_DB/* > {log} 2>&1
        spacepharer createsetdb {input.DB} {output.phage_DB} {params.tmp_folder} --threads {threads} >> {log} 2>&1
        spacepharer createsetdb {input.DB} {output.phage_control_DB} {params.tmp_folder} --reverse-fragments 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_phage:
    input:
        spacer_DB=OUTPUT_DIR + "spacepharer/DB_{type}/querysetDB",
        phage_DB=OUTPUT_DIR + "spacepharer/phage_DB/targetsetDB",
        phage_control_DB=OUTPUT_DIR + "spacepharer/phage_DB/controlsetDB",
    output:
        result=OUTPUT_DIR + "spacepharer/predicted_{type}_phage_matches.tsv",
    params:
        tmp_folder=OUTPUT_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_{type}_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_{type}_phage.txt"
    shell:
        """
        spacepharer predictmatch {input.spacer_DB} {input.phage_DB} {input.phage_control_DB} {output.result} {params.tmp_folder} --threads {threads} > {log} 2>&1
        rm -r {params.tmp_folder} >> {log} 2>&1
        """


rule spacepharer_plasmid_setup:
    input:
        DB=config["spacepharer_plasmid_database"],
    output:
        DB=OUTPUT_DIR + "spacepharer/plasmid_DB/targetsetDB",
        control_DB=OUTPUT_DIR + "spacepharer/plasmid_DB/controlsetDB",
    params:
        tmp_folder=OUTPUT_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_plasmid_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_plasmid_setup.txt"
    shell:
        """
        plasmid_DB=$(dirname {output.DB})
        rm -f $plasmid_DB/* > {log} 2>&1
        spacepharer createsetdb {input.DB} {output.DB} {params.tmp_folder} --threads {threads} >> {log} 2>&1
        spacepharer createsetdb {input.DB} {output.control_DB} {params.tmp_folder} --reverse-fragments 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_plasmid:
    input:
        phage_DB=OUTPUT_DIR + "spacepharer/plasmid_DB/targetsetDB",
        phage_control_DB=OUTPUT_DIR + "spacepharer/plasmid_DB/controlsetDB",
        spacer_DB=OUTPUT_DIR + "spacepharer/DB_{type}/querysetDB",
    output:
        result=OUTPUT_DIR + "spacepharer/predicted_{type}_plasmid_matches.tsv",
    params:
        tmp_folder=OUTPUT_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_{type}_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_{type}_phage.txt"
    shell:
        """
        spacepharer predictmatch {input.spacer_DB} {input.phage_DB} {input.phage_control_DB} {output.result} {params.tmp_folder} --threads {threads} > {log} 2>&1
        rm -r {params.tmp_folder} >> {log} 2>&1
        """
