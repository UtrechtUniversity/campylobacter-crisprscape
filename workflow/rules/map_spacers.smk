### Map spacers to putative targets/protospacers
## 2. By using SpacePHARER


rule spacepharer_spacer_setup:
    input:
        spacers="results/spacers-final.fasta",
    output:
        spacer_db="results/spacepharer/DB_CRISPR/querysetDB",
    params:
        tmp_folder=subpath(output.spacer_db, parent=True),
    conda:
        "../envs/spacepharer.yaml"
    threads: config["spacepharer"]["threads"]
    resources:
        mem_mb=int(config["spacepharer"]["memory"]),
        walltime=int(config["spacepharer"]["time"]),
        runtime=int(config["spacepharer"]["time"]),
    log:
        "log/spacepharer/spacepharer_spacer_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_spacer_setup.txt"
    shell:
        r"""
rm -rf {params.tmp_folder}/* > {log} 2>&1

spacepharer createsetdb {input.spacers} {output.spacer_db}\
 "{params.tmp_folder}/tmpFolder" --extractorf-spacer 1\
 --threads {threads} >> {log} 2>&1
        """


## First to bacteriophages


rule download_phage_database:
    output:
        phage_dir=directory("resources/phagescope"),
        phage_archives=expand(
            "resources/phagescope/{database}.tar.gz",
            database=config["PhageScope_databases"],
        ),
        phage_fasta=expand(
            "resources/phagescope/{database}.fasta",
            database=config["PhageScope_databases"],
        ),
        combined_meta="resources/phagescope/phagescope_metadata.tsv",
    params:
        databases=config["PhageScope_databases"],
    conda:
        "../envs/bash.yaml"
    threads: config["download_spacepharer_databases"]["threads"]
    resources:
        mem_mb=int(config["download_spacepharer_databases"]["memory"]),
        walltime=int(config["download_spacepharer_databases"]["time"]),
        runtime=int(config["download_spacepharer_databases"]["time"]),
    log:
        out="log/download_phage_database.out",
        err="log/download_phage_database.err",
    benchmark:
        "log/benchmark/download_phage_database.txt"
    script:
        "../scripts/download_phage_database.sh"


rule spacepharer_phage_setup:
    input:
        db=collect(
            "resources/phagescope/{database}.fasta",
            database=[
                "DDBJ",
                "EMBL",
                "Genbank",
                "GPD",
                "GVD",
                "MGV",
                "PhagesDB",
                "RefSeq",
                "TemPhD",
            ],
        ),
    output:
        phage_db="results/spacepharer/phage_DB/targetsetDB",
        phage_control_db="results/spacepharer/phage_DB/controlsetDB",
    params:
        tmp_folder=subpath(output.phage_db, ancestor=2),
    conda:
        "../envs/spacepharer.yaml"
    threads: config["spacepharer"]["threads"]
    resources:
        mem_mb=int(config["spacepharer"]["memory"]),
        walltime=int(config["spacepharer"]["time"]),
        runtime=int(config["spacepharer"]["time"]),
    log:
        "log/spacepharer/spacepharer_phage_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_setup.txt"
    shell:
        r"""
phage_db=$(dirname {output.phage_db})
rm -rf ${{phage_db}}/* > {log} 2>&1

spacepharer createsetdb {input.db} {output.phage_db}\
 "{params.tmp_folder}/tmpFolder" --threads {threads} >> {log} 2>&1
spacepharer createsetdb {input.db} {output.phage_control_db}\
 "{params.tmp_folder}/tmpFolder" --reverse-fragments 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_phage:
    input:
        spacer_db="results/spacepharer/DB_CRISPR/querysetDB",
        phage_db="results/spacepharer/phage_DB/targetsetDB",
        phage_control_db="results/spacepharer/phage_DB/controlsetDB",
    output:
        result="results/spacepharer/predicted_phage_matches.tsv",
        result_sanitised="results/spacepharer/predicted_phage_matches_san.tsv",
    params:
        tmp_folder=subpath(output.result, parent=True),
    conda:
        "../envs/spacepharer.yaml"
    threads: config["spacepharer"]["threads"]
    resources:
        mem_mb=int(config["spacepharer"]["memory"]),
        walltime=int(config["spacepharer"]["time"]),
        runtime=int(config["spacepharer"]["time"]),
    log:
        "log/spacepharer/spacepharer_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_phage.txt"
    shell:
        r"""
spacepharer predictmatch {input.spacer_db} {input.phage_db}\
 {input.phage_control_db} {output.result} {params.tmp_folder}\
 --threads {threads} > {log} 2>&1

grep -v "#" {output.result} > {output.result_sanitised}
        """


## Then also to plasmids


rule download_plasmid_database:
    output:
        plasmid_dir=directory("resources/PLSDB"),
        plasmid_fasta="resources/PLSDB/sequences.fasta",
        plasmid_nuccore="resources/PLSDB/nuccore.csv",
        plasmid_taxonomy="resources/PLSDB/taxonomy.csv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        out="log/download_plasmid_database.out",
        err="log/download_plasmid_database.err",
    benchmark:
        "log/benchmark/download_plasmid_database.txt"
    script:
        "../scripts/download_plasmid_database.sh"


rule spacepharer_plasmid_setup:
    input:
        db="resources/PLSDB/sequences.fasta",
    output:
        db="results/spacepharer/plasmid_DB/targetsetDB",
        control_db="results/spacepharer/plasmid_DB/controlsetDB",
    params:
        tmp_folder=subpath(output.db, parent=True),
    conda:
        "../envs/spacepharer.yaml"
    threads: config["spacepharer"]["threads"]
    resources:
        mem_mb=int(config["spacepharer"]["memory"]),
        walltime=int(config["spacepharer"]["time"]),
        runtime=int(config["spacepharer"]["time"]),
    log:
        "log/spacepharer/spacepharer_plasmid_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_plasmid_setup.txt"
    shell:
        r"""
plasmid_db=$(dirname {output.db})
rm -rf ${{plasmid_db}}/* > {log} 2>&1

spacepharer createsetdb {input.db} {output.db} {params.tmp_folder}\
 --threads {threads} >> {log} 2>&1
spacepharer createsetdb {input.db} {output.control_db} {params.tmp_folder}\
 --reverse-fragments 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_plasmid:
    input:
        phage_db="results/spacepharer/plasmid_DB/targetsetDB",
        phage_control_db="results/spacepharer/plasmid_DB/controlsetDB",
        spacer_db="results/spacepharer/DB_CRISPR/querysetDB",
    output:
        result="results/spacepharer/predicted_plasmid_matches.tsv",
        result_sanitised="results/spacepharer/predicted_plasmid_matches_san.tsv",
    params:
        tmp_folder=subpath(output.result, parent=True),
    conda:
        "../envs/spacepharer.yaml"
    threads: config["spacepharer"]["threads"]
    resources:
        mem_mb=int(config["spacepharer"]["memory"]),
        walltime=int(config["spacepharer"]["time"]),
        runtime=int(config["spacepharer"]["time"]),
    log:
        "log/spacepharer/spacepharer_plasmid.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_plasmid.txt"
    shell:
        r"""
spacepharer predictmatch {input.spacer_db} {input.phage_db}\
 {input.phage_control_db} {output.result} {params.tmp_folder}\
  --threads {threads} > {log} 2>&1

grep -v "#" {output.result} > {output.result_sanitised}
        """


rule create_spacepharer_table:
    input:
        phage="results/spacepharer/predicted_phage_matches_san.tsv",
        phage_meta="resources/phagescope/phagescope_metadata.tsv",
        plasmid="results/spacepharer/predicted_plasmid_matches_san.tsv",
        plasmid_nuccore="resources/PLSDB/nuccore.csv",
        plasmid_taxonomy="resources/PLSDB/taxonomy.csv",
    output:
        phage="results/phage_matches.tsv",
        plasmid="results/plasmid_matches.tsv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/create_spacepharer_table.txt",
    script:
        "../scripts/create_spacepharer_table.sh"


## 2. By using KMA to the input genomes


rule kma_indexing:
    input:
        spacers="results/spacers/most_common-final.fasta",
    output:
        indexed_spacers="results/kma/spacer_DB/spacers.name",
    params:
        subpath(output.indexed_spacers, strip_suffix=".name"),
    conda:
        "../envs/kma.yaml"
    threads: config["kma"]["threads"]
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/kma/kma_index.txt",
    benchmark:
        "log/benchmark/kma/kma_index.txt"
    shell:
        r"""
kma index -i {input.spacers} -o {params} > {log} 2>&1
        """


rule mask_crispr_arrays:
    input:
        assembly="resources/ATB/assemblies-concatenated/{batch}.fasta",
        gene_locations="results/cctyper-parse/{batch}/CRISPR-Cas.bed",
    output:
        temp("resources/ATB/assemblies-concatenated/{batch}-masked_crisprs.fasta"),
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/mask_crispr_arrays/{batch}.txt",
    benchmark:
        "log/benchmark/mask_crispr_arrays/{batch}.txt"
    shell:
        """
maskFastaFromBed -fi {input.assembly} -bed {input.gene_locations}\
 -fo {output} > {log} 2>&1
        """


rule kma:
    input:
        genomes=expand(
            "resources/ATB/assemblies-concatenated/{batch}-masked_crisprs.fasta",
            batch=BATCHES,
        ),
        indexed_spacers="results/kma/spacer_DB/spacers.name",
    output:
        frag="results/kma/CRISPR-Cas.frag.gz",
        aln="results/kma/CRISPR-Cas.aln",
        fsa="results/kma/CRISPR-Cas.fsa",
        res="results/kma/CRISPR-Cas.res",
        sam="results/kma/CRISPR-Cas.sam",
    params:
        output=subpath(output[0], strip_suffix=".frag.gz"),
        indexed_spacers=subpath(input.indexed_spacers, strip_suffix=".name"),
    conda:
        "../envs/kma.yaml"
    threads: config["kma"]["threads"]
    resources:
        mem_mb=int(config["kma"]["memory"]),
        walltime=int(config["kma"]["time"]),
        runtime=int(config["kma"]["time"]),
    log:
        "log/kma/kma.txt",
    benchmark:
        "log/benchmark/kma/kma.txt"
    shell:
        r"""
kma -hmm -sam 2308 -i {input.genomes} -o {params.output}\
 -t_db "{params.indexed_spacers}" -t {threads} > {output.sam} 2> {log}
        """


rule collect_kma:
    input:
        "results/kma/CRISPR-Cas.frag.gz",
    output:
        "results/kma/CRISPR_alignment.tsv",
    conda:
        "../envs/bash.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/kma/collect_kma.txt",
    benchmark:
        "log/benchmark/kma/collect_kma.txt"
    shell:
        r"""
echo -e "spacer\ttarget_contig\tstart\tend" > {output}
zless {input} | cut -f 6-9 >> {output}
        """


rule calculate_kma_mismatches:
    input:
        sam=rules.kma.output.sam,
    output:
        table="results/kma/CRISPR_mismatches.tsv",
    conda:
        "../envs/pyfaidx_pandas.yaml"
    threads: 1
    resources:
        mem_mb=int(config["default_job"]["memory"]),
        walltime=int(config["default_job"]["time"]),
        runtime=int(config["default_job"]["time"]),
    log:
        "log/kma/calculate_kma_mismatches.txt",
    benchmark:
        "log/benchmark/kma/calculate_kma_mismatches.txt"
    script:
        "../scripts/calculate_mapping_mismatches.py"
