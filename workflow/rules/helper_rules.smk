## General helper functions: define input and output
from pathlib import Path

# Use Python functions to automatically detect batches of genomes fasta files
# in the input directory as 'BATCHES'
BATCH_PATHS = list(Path("data/tmp/assemblies").glob("atb.assembly.*"))
for batch in BATCH_PATHS:
    assert Path(batch).is_dir(), f"Batches must be directories, got {batch}"

BATCHES = [batch.name for batch in BATCH_PATHS]


## Helper rules (not fitting any particular goal)


rule concatenate_batches:
    input:
        "data/tmp/assemblies/{batch}",
    output:
        temp("data/tmp/assemblies/{batch}.fasta"),
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
