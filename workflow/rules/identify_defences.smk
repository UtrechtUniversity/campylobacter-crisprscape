## Identify anti-phage defence systems using PADLOC


rule download_padloc_database:
    output:
        WORK_DIR + "padloc/database",
    conda:
        "../envs/padloc.yaml"
    threads: 1
    log:
        "log/download_padloc_database.txt",
    benchmark:
        "log/benchmark/download_padloc_database.txt"
    shell:
        """
padloc --db-install v2.0.0 > {log} 2>&1
touch {output}
        """


rule padloc:
    input:
        batch=WORK_DIR + "assemblies/{batch}/",
        db=WORK_DIR + "padloc/database",
    output:
        WORK_DIR + "padloc/{batch}/complete",
    conda:
        "../envs/padloc.yaml"
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
