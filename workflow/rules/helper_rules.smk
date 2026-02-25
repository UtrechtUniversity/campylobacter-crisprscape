## General helper functions: define input and output
from pathlib import Path

# Use Python functions to automatically detect batches of genomes fasta files
# in the input directory as 'BATCHES'
BATCH_PATHS = list(Path("resources/ATB/assemblies").glob("atb.assembly.*"))

# Make sure there is at least one batch directory:
assert len(BATCH_PATHS) > 0, (
    "-- No input (batch) directories found in resources/ATB/assemblies.\n"
    "Please run the script bin/prepare_genomes.sh to prepare input. --\n"
)

# And make sure they are actually directories
for batch in BATCH_PATHS:
    assert Path(batch).is_dir(), f"-- Batches must be directories, got {batch} --"

BATCHES = [batch.name for batch in BATCH_PATHS]


## Helper rules (not fitting any particular goal)


rule concatenate_batches:
    input:
        "resources/ATB/assemblies/{batch}",
    output:
        temp("resources/ATB/assemblies/{batch}.fasta"),
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/concatenate_{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_{batch}.txt"
    shell:
        """
cat {input}/*.fa > {output} 2> {log}
        """


rule convert_bakta_annotations:
    input:
        batch_dir="resources/ATB/annotations/{batch}",
    output:
        gff=directory("resources/ATB/annotations/{batch}/gff"),
        tsv=directory("resources/ATB/annotations/{batch}/tsv"),
    conda:
        "../envs/bakta.yaml"
    threads: config["bakta_convert"]["threads"]
    log:
        out="log/convert_bakta_annotations/{batch}.out",
        err="log/convert_bakta_annotations/{batch}.err",
    benchmark:
        "log/benchmark/convert_bakta_annotations/{batch}.txt"
    script:
        "../scripts/convert_bakta_json.sh"
