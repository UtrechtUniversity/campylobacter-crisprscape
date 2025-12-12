### Classify genomes
## 1: determine multilocus sequence type (MLST)


rule download_mlst_database:
    output:
        WORK_DIR + "mlst/campylobacter.db",
    params:
        species=config["mlst"]["species"],
    conda:
        "../envs/pymlst.yaml"
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
        "../envs/pymlst.yaml"
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


## 2. identify whether contig derive from a chromosome, plasmid or virus
# Using both geNomad (chromosome/plasmid/virus)


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
        "../envs/genomad.yaml"
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
        "../envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_genomad_predictions.txt",
    benchmark:
        "log/benchmark/collect_genomad_predictions.txt"
    script:
        "../scripts/collect_genomad_predictions.R"


# And Jaeger (virus (phage/prophage) or not)


rule jaeger:
    input:
        batch=WORK_DIR + "assemblies/{batch}/",
    output:
        WORK_DIR + "jaeger/{batch}/complete",
    conda:
        "../envs/jaeger.yaml"
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
        "../envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_jaeger_{batch}.txt",
    benchmark:
        "log/benchmark/collect_jaeger_{batch}.txt"
    script:
        "../scripts/collect_jaeger_batch.R"


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
        "../scripts/collect_jaeger_predictions.sh"
