## General helper functions: define input and output
from pathlib import Path

WORK_DIR = config["working_directory"]
OUTPUT_DIR = config["output_directory"]

# Use Python functions to automatically detect batches of genomes fasta files
# in the input directory as 'BATCHES'
BATCH_PATHS = list((Path(WORK_DIR) / "assemblies").glob("atb.assembly.*"))
for batch in BATCH_PATHS:
    assert Path(batch).is_dir(), f"Batches must be directories, got {batch}"

BATCHES = [batch.name for batch in BATCH_PATHS]


## Helper rules (not fitting any particular goal)


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
