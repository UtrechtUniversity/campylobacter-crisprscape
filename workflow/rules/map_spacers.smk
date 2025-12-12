### Map spacers to putative targets/protospacers
## 2. By using SpacePHARER

rule spacepharer_spacer_setup:
    input:
        spacers=WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        spacer_DB=WORK_DIR + "spacepharer/DB_CRISPR/querysetDB",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
    conda:
        "../envs/spacepharer.yaml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_spacer_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_spacer_setup.txt"
    shell:
        """
spacer_DB=$(dirname {output.spacer_DB})
rm -rf $spacer_DB/* > {log} 2>&1

spacepharer createsetdb {input.spacers} {output.spacer_DB} {params.tmp_folder}\
 --extractorf-spacer 1 --threads {threads} >> {log} 2>&1
        """


rule spacepharer_phage_setup:
    output:
        phage_DB=WORK_DIR + "spacepharer/phage_DB/targetsetDB",
        phage_control_DB=WORK_DIR + "spacepharer/phage_DB/controlsetDB",
    params:
        tmp_folder=WORK_DIR + "spacepharer/tmpFolder",
        DB=config["spacepharer_phage_database"] + "*.fasta",
    conda:
        "../envs/spacepharer.yaml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_phage_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_setup.txt"
    shell:
        """
phage_DB=$(dirname {output.phage_DB})
rm -rf $phage_DB/* > {log} 2>&1

spacepharer createsetdb {params.DB} {output.phage_DB} {params.tmp_folder}\
 --threads {threads} >> {log} 2>&1
spacepharer createsetdb {params.DB} {output.phage_control_DB}\
 {params.tmp_folder} --reverse-fragments 1 --threads {threads} >> {log} 2>&1
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
        "../envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_phage.txt"
    shell:
        """
spacepharer predictmatch {input.spacer_DB} {input.phage_DB}\
 {input.phage_control_DB} {output.result} {params.tmp_folder}\
 --threads {threads} > {log} 2>&1

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
        "../envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_plasmid_setup.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_plasmid_setup.txt"
    shell:
        """
plasmid_DB=$(dirname {output.DB})
rm -f $plasmid_DB/* > {log} 2>&1

spacepharer createsetdb {input.DB} {output.DB} {params.tmp_folder}\
 --threads {threads} >> {log} 2>&1
spacepharer createsetdb {input.DB} {output.control_DB} {params.tmp_folder}\
 --reverse-fragments 1 --threads {threads} >> {log} 2>&1
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
        "../envs/spacepharer.yml"
    threads: 48
    log:
        "log/spacepharer/spacepharer_phage.txt",
    benchmark:
        "log/benchmark/spacepharer/spacepharer_phage.txt"
    shell:
        """
spacepharer predictmatch {input.spacer_DB} {input.phage_DB}\
 {input.phage_control_DB} {output.result} {params.tmp_folder}\
  --threads {threads} > {log} 2>&1

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
        "../scripts/create_spacepharer_table.sh"


## 2. By using KMA to the input genomes

rule kma_indexing:
    input:
        spacers=WORK_DIR + "crispridentify/all_spacers.fa",
    output:
        indexed_spacers=WORK_DIR + "kma/spacer_DB/spacers.name",
    params:
        WORK_DIR + "kma/spacer_DB/spacers",
    conda:
        "../envs/kma.yaml"
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
        "../envs/kma.yaml"
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
