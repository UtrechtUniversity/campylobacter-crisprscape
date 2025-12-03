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
WORK_DIR = config["working_directory"]

BATCH_PATHS = list((Path(WORK_DIR) / "assemblies").glob("atb.assembly.*"))
for batch in BATCH_PATHS:
    assert Path(batch).is_dir(), f"Batches must be directories, got {batch}"

BATCHES = [batch.name for batch in BATCH_PATHS]

OUTPUT_DIR = config["output_directory"]


### Step 2: Specify output files ###


rule all:
    input:
        # Multilocus Sequence Types (ST) for Campylobacter
        OUTPUT_DIR + "mlst_table.tsv",
        # Virus and plasmid predictions per contig
        OUTPUT_DIR + "genomad_predictions.csv",
        OUTPUT_DIR + "jaeger_predictions.csv",
        # Phage defence systems
        OUTPUT_DIR + "padloc_table.csv",
        # Concatenated CCTyper output
        expand(
            WORK_DIR + "cctyper/{batch}/{filename}-{batch}.tab",
            batch=BATCHES,
            filename=[
                "CRISPR_Cas",
                "crisprs_all",
                "crisprs_near_cas",
                "crisprs_orphan",
                "cas_operons",
            ],
        ),
        expand(WORK_DIR + "cctyper/{batch}/all_spacers-{batch}.fa", batch=BATCHES),
        # Combined CCTyper output as CSV + BED files
        expand(WORK_DIR + "cctyper/{batch}/parsed", batch=BATCHES),
        # CCTyper CRISPR spacer cluster analysis report
        WORK_DIR + "cctyper/spacer_cluster_summary.tsv",
        # Cluster unique CRISPR spacers
        WORK_DIR + "cctyper/all_spacers-clustered.clstr",
        WORK_DIR + "crispridentify/all_spacers-clustered.clstr",
        OUTPUT_DIR + "all_spacers_table_identify.tsv",
        OUTPUT_DIR + "all_spacers_table.tsv",
        # Extracted CRISPR arrays (as fasta)
        expand(WORK_DIR + "arrays/{batch}/complete", batch=BATCHES),
        #CRISPRidentify output
        expand(WORK_DIR + "crispridentify/{batch}/complete", batch=BATCHES),
        #concatenated CRISPRidentify output
        WORK_DIR + "crispridentify/complete_summary.csv",
        WORK_DIR + "crispridentify/all_spacers.fa",
        #merged CRISPRidentify and CCtyper output
        OUTPUT_DIR + "all_CRISPRS_with_identify.tab",
        #spacepharer output
        OUTPUT_DIR + "phage_matches.tsv",
        OUTPUT_DIR + "plasmid_matches.tsv",
        #KMA output
        WORK_DIR + "kma/output/CRISPR.frag.gz",
        WORK_DIR + "kma/CRISPR_alignment",


### Step 3: Define processing steps that generate the output ###


rule download_mlst_database:
    output:
        WORK_DIR + "mlst/campylobacter.db",
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
        batch=WORK_DIR + "assemblies/{batch}/",
        db=WORK_DIR + "mlst/campylobacter.db",
    output:
        WORK_DIR + "mlst/{batch}/complete",
    conda:
        "envs/pymlst.yaml"
    threads: config["mlst"]["threads"]
    log:
        "log/mlst/{batch}.txt",
    benchmark:
        "log/benchmark/mlst/{batch}.txt"
    shell:
        """
find -L {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 claMLST search {input.db} {{}} -o "$(dirname {output})/{{/.}}.txt" > {log} 2>&1

touch {output}
        """


rule concatenate_mlst_batches:
    input:
        WORK_DIR + "mlst/{batch}/complete",
    output:
        WORK_DIR + "mlst/{batch}-concatenated.tsv",
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
 'tail -n +2 {{}} | cut -f 1-2 >> {output}'
        """


rule concatenate_mlst_all:
    input:
        expand(WORK_DIR + "mlst/{batch}-concatenated.tsv", batch=BATCHES),
    output:
        OUTPUT_DIR + "mlst_table.tsv",
    threads: 1
    log:
        "log/concatenate_mlst_all.txt",
    benchmark:
        "log/benchmark/concatenate_mlst_all.txt"
    shell:
        """
batches=( {input} )
head -1 ${{batches[0]}} > {output}
sed --separate 1d ${{batches[@]}} >> {output}
        """


rule crisprcastyper:
    input:
        batch=WORK_DIR + "assemblies/{batch}/",
    output:
        WORK_DIR + "cctyper/{batch}/complete",
    params:
        out_dir=WORK_DIR + "cctyper/{batch}/",
    conda:
        "envs/cctyper.yaml"
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


rule download_padloc_database:
    output:
        WORK_DIR + "padloc/database",
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
        batch=WORK_DIR + "assemblies/{batch}/",
        db=WORK_DIR + "padloc/database",
    output:
        WORK_DIR + "padloc/{batch}/complete",
    conda:
        "envs/padloc.yaml"
    threads: config["padloc"]["threads"]
    log:
        "log/padloc/{batch}.txt",
    benchmark:
        "log/benchmark/padloc/{batch}.txt"
    shell:
        """
find -L {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 'mkdir -p "$(dirname {output})/{{/.}}" && padloc --cpu 1 --fna {{}} --outdir "$(dirname {output})/{{/.}}"' > {log} 2>&1

touch {output}
        """


rule concatenate_padloc_batches:
    input:
        WORK_DIR + "padloc/{batch}/complete",
    output:
        WORK_DIR + "padloc/{batch}-concatenated.csv",
    threads: config["padloc"]["threads"]
    log:
        "log/concatenate_padloc/{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_padloc/{batch}.txt"
    shell:
        """
file_array=( $(find $(dirname {input}) -mindepth 2 -maxdepth 2 -type f -name "*_padloc.csv") )
head -1 ${{file_array[0]}} > {output}
parallel --jobs {threads} --retry-failed --halt='now,fail=1'\
 'tail -n +2 {{}} >> {output}' ::: ${{file_array[@]}}
        """


rule concatenate_padloc_all:
    input:
        expand(WORK_DIR + "padloc/{batch}-concatenated.csv", batch=BATCHES),
    output:
        OUTPUT_DIR + "padloc_table.csv",
    threads: 1
    log:
        "log/concatenate_padloc_all.txt",
    benchmark:
        "log/benchmark/concatenate_padloc_all.txt"
    shell:
        """
batches=( {input} )
head -1 ${{batches[0]}} > {output}
sed --separate 1d ${{batches[@]}} >> {output}
        """


rule parse_cctyper:
    input:
        WORK_DIR + "cctyper/{batch}/complete",
    output:
        WORK_DIR + "cctyper/{batch}/parsed",
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
        WORK_DIR + "cctyper/{batch}/parsed",
    output:
        WORK_DIR + "cctyper/{batch}/subseq",
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
        cctyper=WORK_DIR + "cctyper/{batch}/complete",
        parser=WORK_DIR + "cctyper/{batch}/parsed",
    output:
        crispr_cas=WORK_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.tab",
        crisprs_all=WORK_DIR + "cctyper/{batch}/crisprs_all-{batch}.tab",
        crisprs_near_cas=WORK_DIR + "cctyper/{batch}/crisprs_near_cas-{batch}.tab",
        crisprs_orphan=WORK_DIR + "cctyper/{batch}/crisprs_orphan-{batch}.tab",
        spacers=WORK_DIR + "cctyper/{batch}/all_spacers-{batch}.fa",
        cas_putative=temp(WORK_DIR + "cctyper/{batch}/cas_operons_putative-{batch}.tab"),
        cas=WORK_DIR + "cctyper/{batch}/cas_operons-{batch}.tab",
        csv=WORK_DIR + "cctyper/{batch}/CRISPR-Cas-{batch}.csv",
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
        WORK_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.tab",
    output:
        WORK_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.bed",
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
        batch=WORK_DIR + "assemblies/{batch}/",
        bed=WORK_DIR + "cctyper/{batch}/CRISPR_Cas-{batch}.bed",
    output:
        WORK_DIR + "arrays/{batch}/complete",
    params:
        out_dir=WORK_DIR + "arrays/{batch}/",
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
        expand(WORK_DIR + "cctyper/{batch}/all_spacers-{batch}.fa", batch=BATCHES),
    output:
        WORK_DIR + "cctyper/all_spacers.fa",
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
        WORK_DIR + "cctyper/all_spacers.fa",
    output:
        clusters=expand(
            WORK_DIR + "cctyper/all_spacers-clustered-{cutoff}.clstr",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        spacers=expand(
            WORK_DIR + "cctyper/all_spacers-clustered-{cutoff}",
            cutoff=[1, 0.96, 0.93, 0.9, 0.87, 0.84, 0.81],
        ),
        summary=WORK_DIR + "cctyper/spacer_cluster_summary.tsv",
    params:
        WORK_DIR=WORK_DIR + "cctyper/",
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
    {params.WORK_DIR}\
    {params.log_dir} > {log} 2>&1
        """


rule cluster_unique_spacers:
    input:
        WORK_DIR + "cctyper/all_spacers.fa",
    output:
        clusters=WORK_DIR + "cctyper/all_spacers-clustered.clstr",
        spacers=WORK_DIR + "cctyper/all_spacers-clustered",
        distribution=WORK_DIR + "cctyper/all_spacers-clustered-distribution.tsv",
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
 -i {input} -o {output.spacers} > {log} 2>&1

plot_len1.pl {output.clusters}\
 1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999\
 1-10,11-20,21-25,26-30,31-35,36-40,41-50,51-60,61-70,71-999999\
 > {output.distribution}
        """


rule create_crispr_cluster_table:
    input:
        clstr=WORK_DIR + "cctyper/all_spacers-clustered.clstr",
        fasta=WORK_DIR + "cctyper/all_spacers.fa",
    output:
        OUTPUT_DIR + "all_spacers_table.tsv",
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
        WORK_DIR + "cctyper/{batch}/subseq",
    output:
        WORK_DIR + "crispridentify/{batch}/complete",
    params:
        out_dir=WORK_DIR + "crispridentify/{batch}",
        arrays=WORK_DIR + "cctyper/{batch}",
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
    find ../../{params.arrays}/*/fasta/CRISPR_arrays-with_flanks.fasta -size +0c -print0 | \
    parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
    python CRISPRidentify.py --file {{}} --result_folder "../../{params.out_dir}/{{/.}}" --fasta_report True --strand False > ../../{log} 2>&1   

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
    shell:
        """
    cat {params.spacers_crispr} > {output.spacers_crispr}
    for summary in {params.summary_crispr} ; do header=$(head -n 1 "$summary"); if [ "$header" == "No arrays found" ];
    then
        continue;
    else 
        echo $header | tee {output.summary_crispr};
        break;
    fi
    done
    for summary in {params.summary_crispr} ; do tail -n +2 "$summary" >> {output.summary_crispr} ; done
        """


rule merge_cctyper_identify:
    input:
        identify=WORK_DIR + "crispridentify/complete_summary.csv",
        cctyper=expand(
            WORK_DIR + "cctyper/{batch}/crisprs_all-{batch}.tab", batch=BATCHES
        ),
    output:
        OUTPUT_DIR + "all_CRISPRS_with_identify.tab",
    params:
        tmp1="tmp_file1",
        tmp2="tmp_file2",
    threads: 1
    log:
        "log/merge_cctyper_identify",
    shell:
        """
    first=True
    for summary in {input.cctyper} ; do
        if [ $first == True ];
        then
            cat $summary > {params.tmp1}
            first=False
        else
            tail -n +2 $summary >> {params.tmp1}
        fi
    done

    header=$(head -n 1 {input.identify} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t")
    tail -n +2 {input.identify} | cut -f 1,5,6,7,8,9,10,11,14 -d "," | tr "," "\t" > {params.tmp2}
    first=True
    while read line; do
        if [ $first == True ];
        then
            first=False
            echo -e "$line\t$header" > {output}
        else
            sample=$(echo -e "$line" | cut -f 1)
            start_cc=$(echo -e "$line" | cut -f 3)
            start_id=$(expr "$start_cc" + 1)
            match=$(grep "${{sample}}_$start_id" {params.tmp2} || true)
            if [ -z "$match" ]; then
                echo -e "$line" >> {output}
            else
                while read line2; do
                    if [ "$start_cc" -lt 5000 ];
                    then
                        echo -e "$line\t$match" >> {output}
                    else
                        start=$(echo -e "$line2" | cut -f 2)
                        start=$(expr "$start" + "$start_cc" - 5000)
                        length=$(echo -e "$line2" | cut -f 4)
                        end=$(expr "$length" + "$start" - 1)
                        begin=$(echo -e "$line2" | cut -f 1)
                        rest=$(echo -e "$line2" | cut -f 4-9)
                        echo -e "$line\t$begin\t$start\t$end\t$rest" >> {output}
                    fi
                done <<< "$match"
            fi
        fi
    done < {params.tmp1}
    rm -f {params.tmp1} {params.tmp2}
        """


rule cluster_unique_spacers_crispridentify:
    input:
        WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        clusters=WORK_DIR + "crispridentify/all_spacers-clustered.clstr",
        spacers=WORK_DIR + "crispridentify/all_spacers-clustered",
        distribution=WORK_DIR + "crispridentify/all_spacers-clustered-distribution.tsv",
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
        "envs/pyfaidx.yaml"
    threads: 1
    log:
        "log/create_crispr_cluster_table_identify.txt",
    benchmark:
        "log/benchmark/create_crispr_cluster_table_identify.txt"
    script:
        "bin/make_cluster_table_identify.py"


rule concatenate_batches:
    input:
        WORK_DIR + "assemblies/{batch}",
    output:
        temp(WORK_DIR + "assemblies/{batch}.fasta"),
    threads: 1
    log:
        "log/concatenate_{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_{batch}.txt"
    shell:
        """
cat {input}/*.fa > {output} 2> {log}
        """


rule genomad:
    input:
        fasta=WORK_DIR + "assemblies/{batch}.fasta",
        db=config["genomad_database"],
    output:
        aggregated_classification=WORK_DIR
        + "genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
        plasmid_summary=WORK_DIR
        + "genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
        virus_summary=WORK_DIR
        + "genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
    params:
        work_dir=WORK_DIR + "genomad/{batch}/",
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
 {input.fasta} {params.work_dir} {input.db} > {log} 2>&1
        """


rule collect_genomad_predictions:
    input:
        aggregated_classification=expand(
            WORK_DIR
            + "genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
            batch=BATCHES,
        ),
        plasmid_summary=expand(
            WORK_DIR + "genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
            batch=BATCHES,
        ),
        virus_summary=expand(
            WORK_DIR + "genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
            batch=BATCHES,
        ),
    output:
        OUTPUT_DIR + "genomad_predictions.csv",
    conda:
        "envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_genomad_predictions.txt",
    benchmark:
        "log/benchmark/collect_genomad_predictions.txt"
    script:
        "bin/collect_genomad_predictions.R"


rule jaeger:
    input:
        batch=WORK_DIR + "assemblies/{batch}/",
    output:
        WORK_DIR + "jaeger/{batch}/complete",
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
 jaeger run -p --workers 1 -i {{}} -o $(dirname {output}) --overwrite\
 > {log} 2>&1 ::: {input.batch}/*.fa

touch {output}
        """


rule collect_jaeger_batch:
    input:
        WORK_DIR + "jaeger/{batch}/complete",
    output:
        WORK_DIR + "jaeger/{batch}/jaeger-{batch}.csv",
    params:
        batch="{batch}",
    conda:
        "envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_jaeger_{batch}.txt",
    benchmark:
        "log/benchmark/collect_jaeger_{batch}.txt"
    script:
        "bin/collect_jaeger_batch.R"


rule collect_jaeger_predictions:
    input:
        expand(WORK_DIR + "jaeger/{batch}/jaeger-{batch}.csv", batch=BATCHES),
    output:
        OUTPUT_DIR + "jaeger_predictions.csv",
    threads: 1
    log:
        "log/collect_jaeger_predictions.txt",
    benchmark:
        "log/benchmark/collect_jaeger_predictions.txt"
    script:
        "bin/collect_jaeger_predictions.sh"


rule spacepharer_spacer_setup:
    input:
        spacers=WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        spacer_DB=WORK_DIR + "spacepharer/DB_CRISPR/querysetDB",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_spacer_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_spacer_setup.txt"
    shell:
        """
        spacer_DB=$(dirname {output.spacer_DB})
        rm -rf $spacer_DB/* > {log} 2>&1 
        spacepharer createsetdb {input.spacers} {output.spacer_DB} {params.tmp_folder} --extractorf-spacer 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_phage_setup:
    output:
        phage_DB=WORK_DIR + "spacepharer/phage_DB/targetsetDB",
        phage_control_DB=WORK_DIR + "spacepharer/phage_DB/controlsetDB",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
        DB=config["spacepharer_phage_database"] + "*.fasta",
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
        spacepharer createsetdb {params.DB} {output.phage_DB} {params.tmp_folder} --threads {threads} >> {log} 2>&1
        spacepharer createsetdb {params.DB} {output.phage_control_DB} {params.tmp_folder} --reverse-fragments 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_phage:
    input:
        spacer_DB=WORK_DIR + "spacepharer/DB_CRISPR/querysetDB",
        phage_DB=WORK_DIR + "spacepharer/phage_DB/targetsetDB",
        phage_control_DB=WORK_DIR + "spacepharer/phage_DB/controlsetDB",
    output:
        result=WORK_DIR + "spacepharer/predicted_phage_matches.tsv",
        result_sanitised=WORK_DIR + "spacepharer/predicted_phage_matches_san.tsv",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_phage.txt"
    shell:
        """
        spacepharer predictmatch {input.spacer_DB} {input.phage_DB} {input.phage_control_DB} {output.result} {params.tmp_folder} --threads {threads} > {log} 2>&1
        grep -v "#" {output.result} > {output.result_sanitised} 
        rm -r {params.tmp_folder} >> {log} 2>&1
        """


rule spacepharer_plasmid_setup:
    input:
        DB=config["spacepharer_plasmid_database"] + "sequences.fasta",
    output:
        DB=WORK_DIR + "spacepharer/plasmid_DB/targetsetDB",
        control_DB=WORK_DIR + "spacepharer/plasmid_DB/controlsetDB",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
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
        phage_DB=WORK_DIR + "spacepharer/plasmid_DB/targetsetDB",
        phage_control_DB=WORK_DIR + "spacepharer/plasmid_DB/controlsetDB",
        spacer_DB=WORK_DIR + "spacepharer/DB_CRISPR/querysetDB",
    output:
        result=WORK_DIR + "spacepharer/predicted_plasmid_matches.tsv",
        result_sanitised=WORK_DIR + "spacepharer/predicted_plasmid_matches_san.tsv",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
    conda:
        "envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_phage.txt"
    shell:
        """
        spacepharer predictmatch {input.spacer_DB} {input.phage_DB} {input.phage_control_DB} {output.result} {params.tmp_folder} --threads {threads} > {log} 2>&1
        grep -v "#" {output.result} > {output.result_sanitised} 
        rm -r {params.tmp_folder} >> {log} 2>&1
        """


rule create_spacepharer_table:
    input:
        phage=WORK_DIR + "spacepharer/predicted_phage_matches_san.tsv",
        meta_phage=config["spacepharer_phage_database"],
        plasmid=WORK_DIR + "spacepharer/predicted_plasmid_matches_san.tsv",
        meta_plasmid=config["spacepharer_plasmid_database"],
    output:
        phage=OUTPUT_DIR + "phage_matches.tsv",
        plasmid=OUTPUT_DIR + "plasmid_matches.tsv",
    threads: 1
    log:
        "log/create_spacepharer_table.txt",
    script:
        "bin/create_spacepharer_table.sh"


rule kma_indexing:
    input:
        spacers=WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        indexed_spacers=WORK_DIR + "kma/spacer_DB/spacers.name",
    params:
        WORK_DIR + "kma/spacer_DB/spacers",
    conda:
        "envs/kma.yaml"
    threads: 12
    log:
        "log/kma/kma_index.txt",
    benchmark:
        "log/benchmark/kma/kma_index.txt"
    shell:
        """
        kma index -i {input.spacers} -o {params} > {log} 2>&1
        """


rule kma:
    input:
        genomes=expand(WORK_DIR + "assemblies/{batch}/", batch=BATCHES),
        indexed_spacers=WORK_DIR + "kma/spacer_DB/spacers.name",
    output:
        WORK_DIR + "kma/output/CRISPR.frag.gz",
    params:
        output=WORK_DIR + "kma/output/CRISPR",
        indexed_spacers=WORK_DIR + "kma/spacer_DB/spacers",
        spacers=WORK_DIR + "crispridentify/all_spacers.fa",
    conda:
        "envs/kma.yaml"
    threads: 24
    log:
        "log/kma/kma.txt",
    benchmark:
        "log/benchmark/kma/kma.txt"
    shell:
        """     
        grep ">" {params.spacers} | cut -f 2 -d ">" | cut -f 1 -d "-" | sort -u > tmp_file
        find -L {input.genomes} -mindepth 1 -maxdepth 1 -type f -name "*.fa" > all_genomes.txt
        genomes=$(grep -x ".*[0-9]\\.fa" all_genomes.txt | grep -v -f tmp_file)
        kma -hmm -i $genomes -o {params.output} -t_db {params.indexed_spacers} > {log} 2>&1
        rm tmp_file all_genomes.txt
        """


rule collect_kma:
    input:
        WORK_DIR + "kma/output/CRISPR.frag.gz",
    output:
        WORK_DIR + "kma/CRISPR_alignment",
    log:
        "log/kma/collect_kma.txt",
    benchmark:
        "log/benchmark/kma/collect_kma.txt"
    shell:
        """
        echo -e "spacer\tgenome" > {output}
        zcat {input} | cut -f 6,7 | cut -f 1 -d " " > tmp_file
        while read line; do
            match=$(echo $line | cut -f 2)
            crispr=$(echo $line | cut -f 1 | cut -f 1,6,7,10,11 -d "_")
            echo -e "$crispr\t$match" >> {output}
        done < tmp_file
        rm tmp_file
        """
