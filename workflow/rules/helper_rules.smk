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

ANNOTATION_PATHS = list(Path("resources/ATB/annotations/").glob("atb.bakta.*"))
assert len(ANNOTATION_PATHS) > 0, (
    "-- No annotation files found in resources/ATB/annotations.\n"
    "Please run the script bin/prepare_genomes.sh to prepare input. --\n"
)

for annotation in ANNOTATION_PATHS:
    assert Path(
        annotation
    ).is_dir(), f"-- Batches must be directories, got {annotation} --"

BATCHES = [batch.name for batch in BATCH_PATHS]
ANNOTATIONS = [annotation.name for annotation in ANNOTATION_PATHS]

## Helper rules (not fitting any particular goal)


rule concatenate_batches:
    input:
        "resources/ATB/assemblies/{batch}",
    output:
        "resources/ATB/assemblies-concatenated/{batch}.fasta",
    conda:
        "../envs/bash.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/concatenate_{batch}.txt",
    benchmark:
        "log/benchmark/concatenate_{batch}.txt"
    shell:
        r"""
cat {input}/*.fa > {output} 2> {log}
        """


rule convert_bakta_annotations:
    input:
        annotation_dir="resources/ATB/annotations/{annotation}",
    output:
        gff=directory("resources/ATB/annotations/{annotation}/gff"),
        tsv=directory("resources/ATB/annotations/{annotation}/tsv"),
    conda:
        "../envs/bakta.yaml"
    threads: config["bakta_convert"]["threads"]
    resources:
        mem_mb=int(config["bakta_convert"]["memory"]),
        walltime=int(config["bakta_convert"]["time"]),
        runtime=int(config["bakta_convert"]["time"]),
    log:
        out="log/convert_bakta_annotations/{annotation}.out",
        err="log/convert_bakta_annotations/{annotation}.err",
    benchmark:
        "log/benchmark/convert_bakta_annotations/{annotation}.txt"
    script:
        "../scripts/convert_bakta_json.sh"


rule calculate_pangenomes:
    input:
        "resources/ATB/annotations/{annotation}/gff",
    output:
        multiext(
            "results/panaroo/{annotation}/",
            entropy="alignment_entropy.csv",
            state="alignment_resume_state.json",
            cds="combined_DNA_CDS.fasta",
            prot_cdhit="combined_protein_cdhit_out.txt",
            prot_cdhit_clstr="combined_protein_cdhit_out.txt.clstr",
            protein="combined_protein_CDS.fasta",
            embl_filt="core_alignment_filtered_header.embl",
            embl="core_alignment_header.embl",
            align="core_gene_alignment.aln",
            align_filt="core_gene_alignment_filtered.aln",
            graph="final_graph.gml",
            data="gene_data.csv",
            pres_abs="gene_presence_absence.csv",
            roary="gene_presence_absence_roary.csv",
            rtab="gene_presence_absence.Rtab",
            pan_ref="pan_genome_reference.fa",
            pre_graph="pre_filt_graph.gml",
            struct="struct_presence_absence.Rtab",
            summary="summary_statistics.txt",
        ),
        aligned=directory("results/panaroo/{annotation}/aligned_gene_sequences"),
        general=directory("results/panaroo/{annotation}"),
    params:
        out_dir=subpath(output[0], parent=True),
    conda:
        "../envs/panaroo.yaml"
    threads: config["panaroo"]["threads"]
    resources:
        mem_mb=int(config["panaroo"]["memory"]),
        walltime=int(config["panaroo"]["time"]),
        runtime=int(config["panaroo"]["time"]),
    log:
        "log/calculate_pangenomes/{annotation}.txt",
    benchmark:
        "log/benchmark/calculate_pangenomes/{annotation}.txt"
    shell:
        """
panaroo -i {input}/*.gff3 -o {params.out_dir} --clean-mode strict -t {threads}\
 -a core --aligner mafft --core_threshold 0.95 --remove-invalid-genes\
 > {log} 2>&1

touch {output.entropy} {output.align} {output.align_filt} {output.embl} {output.embl_filt}
        """


rule merge_pangenomes:
    input:
        expand("results/panaroo/{annotation}/", annotation=ANNOTATIONS),
    output:
        multiext(
            "results/panaroo/merged/",
            cds="combined_DNA_CDS.fasta",
            graph="final_graph.gml",
            data="gene_data.csv",
            pres_abs="gene_presence_absence.csv",
            roary="gene_presence_absence_roary.csv",
            rtab="gene_presence_absence.Rtab",
            pan_ref="pan_genome_reference.fa",
            struct="struct_presence_absence.Rtab",
            summary="summary_statistics.txt",
        ),
    params:
        out_dir=subpath(output[0], parent=True),
    conda:
        "../envs/panaroo.yaml"
    threads: config["panaroo"]["threads"]
    resources:
        mem_mb=int(config["panaroo"]["memory"]),
        walltime=int(config["panaroo"]["time"]),
        runtime=int(config["panaroo"]["time"]),
    log:
        "log/merge_pangenomes.txt",
    benchmark:
        "log/benchmark/merge_pangenomes.txt"
    shell:
        """
panaroo-merge -d {input} -o {params.out_dir} -t {threads}\
 --core_threshold 0.95 > {log} 2>&1
        """


rule collect_contig_lengths:
    input:
        "resources/ATB/assemblies-concatenated/{batch}.fasta",
    output:
        "results/contig_lengths/{batch}.tsv",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/collect_contig_lengths/{batch}.txt",
    benchmark:
        "log/benchmark/collect_contig_lengths/{batch}.txt"
    shell:
        """
seqkit fx2tab -H -l -i -Q {input} | cut -f 1,3 > {output} 2> {log}
        """
