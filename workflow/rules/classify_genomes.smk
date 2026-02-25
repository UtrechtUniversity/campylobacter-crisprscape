### Classify genomes
## 1: determine multilocus sequence type (MLST)


rule download_mlst_database:
    output:
        "resources/mlst/campylobacter.db",
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
        batch="resources/ATB/assemblies/{batch}/",
        db="resources/mlst/campylobacter.db",
    output:
        "results/mlst/{batch}/complete",
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
        "results/mlst/{batch}/complete",
    output:
        "results/mlst/{batch}-concatenated.tsv",
    conda:
        "../envs/bash.yaml"
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
        expand("results/mlst/{batch}-concatenated.tsv", batch=BATCHES),
    output:
        "results/mlst_table.tsv",
    conda:
        "../envs/bash.yaml"
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


rule download_genomad_database:
    output:
        db=directory("resources/genomad_db"),
    conda:
        "../envs/genomad.yaml"
    threads: 1
    log:
        "log/download_genomad_database.txt",
    benchmark:
        "log/benchmark/download_genomad_database.txt"
    shell:
        """
mkdir -p $(dirname {output.db})
genomad download-database $(dirname {output.db}) > {log} 2>&1
        """


rule genomad:
    input:
        fasta="resources/ATB/assemblies/{batch}.fasta",
        db="resources/genomad_db",
    output:
        aggregated_classification="results/genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
        plasmid_summary="results/genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
        virus_summary="results/genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
    params:
        work_dir=subpath(output.aggregated_classification, ancestor=2),
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
            "results/genomad/{batch}/{batch}_aggregated_classification/{batch}_aggregated_classification.tsv",
            batch=BATCHES,
        ),
        plasmid_summary=expand(
            "results/genomad/{batch}/{batch}_summary/{batch}_plasmid_summary.tsv",
            batch=BATCHES,
        ),
        virus_summary=expand(
            "results/genomad/{batch}/{batch}_summary/{batch}_virus_summary.tsv",
            batch=BATCHES,
        ),
    output:
        "results/genomad_predictions.csv",
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
        batch="resources/ATB/assemblies/{batch}/",
    output:
        "results/jaeger/{batch}/complete",
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
        "results/jaeger/{batch}/complete",
    output:
        "results/jaeger/{batch}/jaeger-{batch}.csv",
    params:
        batch=subpath(output[0], parent=True),
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
        expand("results/jaeger/{batch}/jaeger-{batch}.csv", batch=BATCHES),
    output:
        "results/jaeger_predictions.csv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/collect_jaeger_predictions.txt",
    benchmark:
        "log/benchmark/collect_jaeger_predictions.txt"
    script:
        "../scripts/collect_jaeger_predictions.sh"


## 3: dereplicate genomes (or at least: mark duplicates)
# (dRep: 99.99% identity)


rule simplify_checkm:
    input:
        "resources/ATB/checkm2.tsv.gz",
    output:
        "resources/ATB/ATB-CheckM2.csv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/simplify_checkm.txt",
    benchmark:
        "log/benchmark/simplify_checkm.txt"
    shell:
        """
echo "genome,contamination,completeness" >  {output}
zless {input} | tail -n +2 | cut -f 1,3,4 | tr "\t" "," >> {output}
        """


rule dereplicate_genomes:
    input:
        batch="resources/ATB/assemblies/{batch}/",
        checkm=rules.simplify_checkm.output,
    output:
        genome_list=temp("results/drep/{batch}/genome_list.txt"),
        data=directory("results/drep/{batch}/data"),
        data_tables=expand(
            "results/drep/{{batch}}/data_tables/{prefix}.csv",
            prefix=["Bdb", "Cdb", "genomeInformation", "Mdb", "Sdb", "Wdb", "Widb"],
        ),
        dereplicated_genomes=directory("results/drep/{batch}/dereplicated_genomes"),
        figures=directory("results/drep/{batch}/figures"),
        log=directory("results/drep/{batch}/log"),
    params:
        prefix=subpath(output["log"], parent=True),
    conda:
        "../envs/drep.yaml"
    threads: config["drep"]["threads"]
    log:
        "log/drep/{batch}.txt",
    benchmark:
        "log/benchmark/drep/{batch}.txt"
    shell:
        """
find {input.batch} -name "*.fa" -print > {output.genome_list}

dRep dereplicate {params.prefix} -g {output.genome_list}\
 -sa 0.9999 -p {threads} --ignoreGenomeQuality\
 --genomeInfo {input.checkm}\
 --clusterAlg complete -pa 0.99 -nc 0.5 > {log} 2>&1
        """


rule collect_dereplications:
    input:
        cluster=expand("results/drep/{batch}/data_tables/Cdb.csv", batch=BATCHES),
        winner=expand("results/drep/{batch}/data_tables/Wdb.csv", batch=BATCHES),
    output:
        "results/dereplication_table.tsv",
    conda:
        "../envs/tidy_here.yaml"
    threads: 1
    log:
        "log/collect_dereplications.txt",
    benchmark:
        "log/benchmark/collect_dereplications.txt"
    script:
        "../scripts/collect_dereplication_tables.R"
