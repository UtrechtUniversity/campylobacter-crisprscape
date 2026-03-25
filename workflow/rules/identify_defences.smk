## Identify anti-phage defence systems using PADLOC


rule download_padloc_database:
    output:
        directory=directory("resources/padloc_db"),
        cm=directory("resources/padloc_db/cm"),
        cm_meta="resources/padloc_db/cm_meta.txt",
        hmm=directory("resources/padloc_db/hmm"),
        hmm_meta="resources/padloc_db/hmm_meta.txt",
        sys=directory("resources/padloc_db/sys"),
        sys_meta="resources/padloc_db/sys_meta.txt",
        system_info="resources/padloc_db/system_info.md",
    conda:
        "../envs/padloc.yaml"
    threads: 1
    log:
        "log/download_padloc_database.txt",
    benchmark:
        "log/benchmark/download_padloc_database.txt"
    shell:
        r"""
padloc --data resources/padloc_db --db-install v2.0.0 > {log} 2>&1
        """


rule padloc:
    input:
        batch="resources/ATB/assemblies/{batch}/",
        db="resources/padloc_db",
    output:
        "results/padloc/{batch}/complete",
    conda:
        "../envs/padloc.yaml"
    threads: config["padloc"]["threads"]
    log:
        "log/padloc/{batch}.txt",
    benchmark:
        "log/benchmark/padloc/{batch}.txt"
    shell:
        r"""
find -L {input.batch} -mindepth 1 -maxdepth 1 -type f -name "*.fa" -print0 |\
 parallel -0 --jobs {threads} --retry-failed --halt='now,fail=1'\
 'mkdir -p "$(dirname {output})/{{/.}}" && padloc --data {input.db}\
  --cpu 1 --fna {{}} --outdir "$(dirname {output})/{{/.}}"' > {log} 2>&1

touch {output}
        """


rule concatenate_padloc_batches:
    input:
        "results/padloc/{batch}/complete",
    output:
        "results/padloc/{batch}-concatenated.csv",
    conda:
        "../envs/bash.yaml"
    threads: config["padloc"]["threads"]
    log:
        "log/concatenate_padloc/{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_padloc/{batch}.txt"
    shell:
        r"""
file_array=( $(find $(dirname {input}) -mindepth 2 -maxdepth 2 -type f -name "*_padloc.csv") )
head -1 ${{file_array[0]}} > {output}
parallel --jobs {threads} --retry-failed --halt='now,fail=1'\
 'tail -n +2 {{}} >> {output}' ::: ${{file_array[@]}}
        """


rule concatenate_padloc_all:
    input:
        expand("results/padloc/{batch}-concatenated.csv", batch=BATCHES),
    output:
        "results/padloc_table.csv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/concatenate_padloc_all.txt",
    benchmark:
        "log/benchmark/concatenate_padloc_all.txt"
    shell:
        r"""
batches=( {input} )
head -1 ${{batches[0]}} > {output}
sed --separate 1d ${{batches[@]}} >> {output}
        """
