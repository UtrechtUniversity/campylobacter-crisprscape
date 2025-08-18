# Output files

## 1. CRISPR-Cas screening

### 1.1 CCTyper

The results of CCTyper are by default written to:

``` bash
data/tmp/cctyper/[batch_id]/[sample_id]/
```

For example:

``` bash
$ ls -F data/tmp/cctyper/batch_22/SAMN36083307
arguments.tab    cas_operons_putative.tab  CRISPR-Cas.bed  CRISPR_Cas.tab   crisprs.gff           fasta/     Flank.njs  Flank.ntf  genes.tab  plot.png  spacers/
Cas_operons.bed  CRISPR_arrays.bed         CRISPR-Cas.csv  crisprs_all.tab  crisprs_near_cas.tab  Flank.ndb  Flank.not  Flank.nto  hmmer.tab  plot.svg
```

For the `*.tab` files, `plot.[png|svg]`, `Flank.*` and `spacers/`, please see
CCTyper's
[own documentation](https://github.com/Russel88/CRISPRCasTyper?tab=readme-ov-file#output-).
**This comprises the fundamental output of our workflow.**

For other (derived) files, please read below.

### 1.2 Collecting many outputs from CCTyper

After running CCTyper, there is a script
[`bin/cctyper_extender`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/cctyper_extender.py)
that reads these tabular output files  (`*.tab` files) and converts them to
BED files for each CRISPR array (spacers and repeats), CRISPR-Cas system
(CRISPR array with _cas_ genes), or _cas_ operon.
These files contain the contig ID in which the locus was found, its start and
stop positions, the CRISPR array or _cas_ operon ID as listed by CCTyper,
the number of CRISPR repeats or _cas_ genes found,
and whether the _cas_ genes were on the template or reverse complement strand.
For example:

``` bash
$ head data/tmp/cctyper/batch_22/SAMN36083307/*bed
==> data/tmp/cctyper/batch_22/SAMN36083307/Cas_operons.bed <==
SAMN36083307.contig00015        38536   42801   SAMN36083307.contig00015@1      3       -

==> data/tmp/cctyper/batch_22/SAMN36083307/CRISPR_arrays.bed <==
SAMN36083307.contig00015        37630   38392   SAMN36083307.contig00015_1      12      -

==> data/tmp/cctyper/batch_22/SAMN36083307/CRISPR-Cas.bed <==
SAMN36083307.contig00015        37630   42801   SAMN36083307.contig00015_1      12      -
```

Furthermore, the script summarises all characteristics of the CRISPR array and
_cas_ genes in a single file: `CRISPR-Cas.csv`.

The columns that it includes are:

``` bash
Sample_accession    Contig              System                  CRISPR_ID           CRISPR_start    CRISPR_end
Trusted             Repeat_subtype      Repeat_type_probability Consensus_repeat    N_repeats
Repeat_len          Repeat_identity     Spacer_len_avg          Spacer_identity     Spacer_len_sem
Operon_ID           Distance_to_CRISPR  Cas_start               Cas_end             Best_type
Best_score          Genes               Complete_interference   Complete_adaptation Strand_cas
N_genes             Gene_lengths_aa
```

Note that while these are taken from CCTyper, some names have been adjusted
or added.

### 1.3 Use BED files to extract sequences

Using the BED files generated in the step above, we use
[`seqkit`](https://bioinf.shenwei.me/seqkit/usage/#subseq)
to extract the loci of interest from the original assemblies based on the
contig IDs and positions.
For this, we use a custom script:
[`bin/extract_crispr-cas_from_fasta.sh`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/extract_crispr-cas_from_fasta.sh)
.
This script extracts the regions from start to stop, and also with (up to)
5000bp up- and downstream of the region. These are then written to FASTA files
in the subdirectory `fasta/`

Besides, for each locus it also generates FASTA files with only either the
upstream or downstream region. These can be used to infer the genomic
neighbourhood and compare between genomes.

``` bash
$ ls -sh data/tmp/cctyper/batch_22/SAMN36083307/fasta/
total 100K
8.0K Cas_operons-downstream.fasta   12K Cas_operons-with_flanks.fasta   8.0K CRISPR_arrays-upstream.fasta     8.0K CRISPR-Cas.fasta
8.0K Cas_operons.fasta             8.0K CRISPR_arrays-downstream.fasta   12K CRISPR_arrays-with_flanks.fasta  4.0K CRISPR-Cas-upstream.fasta
4.0K Cas_operons-upstream.fasta    4.0K CRISPR_arrays.fasta             8.0K CRISPR-Cas-downstream.fasta       16K CRISPR-Cas-with_flanks.fasta
```

### 1.4 Concatenate batches

We currently split our input samples in batches and to facilitate overall
analyses, we concatenate the output files within batches using the rule
`collect_cctyper`, which in turn invokes
[`bin/concatenate_cctyper_output.sh`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/concatenate_cctyper_output.sh)
and
[`bin/concatenate_cctyper_csv.sh`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/concatenate_cctyper_csv.sh)
as well as a `cat` command to concatenate the FASTA files of spacer sequences

``` bash
find $(dirname {input.cctyper}) -mindepth 3 -maxdepth 3 -name "*.fa"\
 -exec cat {{}} + > {output.spacers} 2>> {log}
```

This creates a set of files with the batch name as suffix, for example:

``` bash
$ ls -sh data/tmp/cctyper/batch_22

4.0K all_spacers-batch_22.fa   4.0K CRISPR-Cas-batch_22.csv        4.0K crisprs_orphan-batch_22.tab
4.0K cas_operons-batch_22.tab  4.0K CRISPR_Cas-batch_22.tab        4.0K crisprs_all-batch_22.tab
4.0K CRISPR_Cas-batch_22.bed   4.0K crisprs_near_cas-batch_22.tab
```

### 1.5 Concatenate overall

To then go to one final file for the whole dataset, the files listed above
are further concatenated in rule `concatenate_all_spacers` with `cat`.

``` bash
$ ls -sh data/tmp/cctyper/all_spacers.fa
22M data/tmp/cctyper/all_spacers.fa

$ head data/tmp/cctyper/all_spacers.fa
>SAMN37284265.contig00007_1:1
TTATTTTTGTCGCTAATTGCACCTAAAGAC
>SAMN37284265.contig00007_1:2
TTTAGAAGAAGAAATATCTTATCTTCAAAG
>SAMN37284265.contig00007_1:3
ATAAGAAGTTAATTTGCTTACTTTCAAGTA
>SAMN37284265.contig00007_1:4
GTGGGCATAAGCAACTTAATGGCAGGATTT
>SAMN37284265.contig00007_1:5
TAAAGCTGGAAATTTTTCAACTAATTTACC
```

!!! Warning "Warning: missing piece"

    Note that we currently do not have a separate step to concatenate
    the CSV files derived from CCTyper.

## Cluster spacers: CCTyper

All CRISPR spacer sequences identified by CCTyper are clustered with
CD-HIT to determine the number of unique spacers. This happens in the rule
`cluster_unique_spacers`, and generates:

``` bash
ls -sh data/tmp/cctyper/all_spacers-clustered*
948K data/tmp/cctyper/all_spacers-clustered
 20M data/tmp/cctyper/all_spacers-clustered.clstr
4.0K data/tmp/cctyper/all_spacers-clustered-distribution.tsv
```

The file without extension is the FASTA file that CD-HIT outputs,
`all_spacers-clustered.clstr` is the list of clusters as determined
by CD-HIT (in FASTA-like format). And the file
`all_spacers-clustered-distribution.tsv` uses the script `plot_len1.pl` that
is bundled with CD-HIT to create a tab-separated table of the number of
members in clusters of various sizes.

Also see the user guide on
[bioinformatics.org](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide)
for details on CD-HIT's output files and extra tools.

### Spacer cluster table

To facilitate further analyses on the spacer level, we have included the
step `create_crispr_cluster_table`, which calls
[`bin/make_cluster_table.py`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/make_cluster_table.py)
to parse the output of CD-HIT and create a tab-separated table that looks
like:

``` bash
$ head data/processed/all_spacers_table.tsv
Genome  Contig  Locus   Cluster Length  Cluster_representative  Sequence        Identity        Strand
SAMN30649373    contig00012     SAMN30649373.contig00012_1:1    0       75nt    SAMN30649373.contig00012_1:1    TATAGCAGTAAAATAGGCTTTAGAATGTATTTTATAAAAGCGTAAAAAATCAATATTATACCGACATTTTGCGAC     NA      NA
SAMN30649373    contig00012     SAMN30649373.contig00012_1:2    0       55nt    SAMN30649373.contig00012_1:1    TATAGCAGTAAAATAGGCTTTAGAATGTATTTTATAAAAGCGTAAAAAATCAATA 100.00% +
SAMN19302369    contig00024     SAMN19302369.contig00024_2:1    1       57nt    SAMN19302369.contig00024_2:1    TATTATATAATATATACTAGTGCATTTTGATAATAATTATCGACTTCAAGTAAAATT       NA      NA
SAMN36876523    contig00043     SAMN36876523.contig00043_2:1    2       56nt    SAMN36876523.contig00043_2:1    ATGCTTGGCTTTTTTGTAGTATATAGCAATGCTATATATCCTAGTATATACTTAAA        NA      NA
SAMN11571935    contig00011     SAMN11571935.contig00011_2:1    3       56nt    SAMN11571935.contig00011_2:1    ATAATGGTAAAATAGTCTTTAGAATGTATTTTAAAAAAGCATAAAAAGCTAATATA        NA      NA
SAMN32601407    contig00036     SAMN32601407.contig00036_2:1    4       56nt    SAMN32601407.contig00036_2:1    ATAATGGTAAAATAGTCTTTAGAATGTATTTTAAATAAACTTAAAAAGCTAATATA        NA      NA
SAMEA1676287    contig00014     SAMEA1676287.contig00014_1:2    5       56nt    SAMEA1676287.contig00014_1:2    TAGAAAATTTTTTATAAAAAATGTATATTTTTTTACCAAAAAATGAGACAAAAGCA        NA      NA
SAMEA1676287    contig00014     SAMEA1676287.contig00014_1:4    5       56nt    SAMEA1676287.contig00014_1:2    TAGAAAATTTTTTATAAAAAATGTATATTTTTTTACCAAAAAATGAGACAAAAGCA        100.00% +
SAMN15582136    contig00006     SAMN15582136.contig00006_1:9    6       55nt    SAMN15582136.contig00006_1:9    CTCAAACTTAGATTTATGGTAAAATTCAAATTCAAATATGACATCTAATAGTGAC NA      NA
```

!!! info "Final output: `data/processed/all_spacers_table.tsv`"
