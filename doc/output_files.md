# Output files

!!! note

    Files that we consider to be the 'final output' of this workflow are
    labelled with a box like this:

    !!! info "Final output: `results/filename`"

    They are stored in the directory `results/`.

    Other output files are often not needed for end-users and are stored in
    subfolders within `results/`.

## 1. CRISPR-Cas screening

### 1.1 CCTyper

The results of CCTyper are by default written to:

``` bash
results/cctyper/[batch_id]/
```

For example:

``` bash
$ ls -F results/cctyper/atb.assembly.incr_release.202408.batch.38/
arguments.tab             crisprs_all.tab crisprs_putative.tab  Flank.ntf  plot.png
cas_operons_putative.tab  crisprs.gff           Flank.ndb       Flank.nto  plot.svg
CRISPR_Cas_putative.tab   crisprs_near_cas.tab  Flank.njs       genes.tab  spacers/
CRISPR_Cas.tab            crisprs_orphan.tab    Flank.not       hmmer.tab
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
$ head -2 results/cctyper-parse/atb.assembly.incr_release.202408.batch.38/*bed
==> results/cctyper-parse/atb.assembly.incr_release.202408.batch.38/Cas_operons.bed <==
SAMEA115501754.contig00002 112331 116642 SAMEA115501754.contig00002@4 4 +
SAMEA115501782.contig00002 121539 125851 SAMEA115501782.contig00002@4 3 +

==> results/cctyper-parse/atb.assembly.incr_release.202408.batch.38/CRISPR_arrays.bed <==
SAMEA115501754.contig00002 116789 117088 SAMEA115501754.contig00002_56 5 +
SAMEA115501782.contig00002 125998 127282 SAMEA115501782.contig00002_80 20 +

==> results/cctyper-parse/atb.assembly.incr_release.202408.batch.38/CRISPR-Cas.bed <==
SAMEA115501754.contig00002 112331 117088 SAMEA115501754.contig00002_56 5 +
SAMEA115501782.contig00002 121539 127282 SAMEA115501782.contig00002_80 20 +
```

Furthermore, the script summarises all characteristics of the CRISPR array and
_cas_ genes in a single file: `results/cctyper-parse/[batch]/CRISPR-Cas.tsv`.

It has 28 columns, which are:

``` bash
$ head -1 results/cctyper-parse/atb.assembly.incr_release.202408.batch.38/CRISPR-Cas.tsv | tr "\t" "\n" | nl
     1 Sample_accession
     2 Contig
     3 System
     4 CRISPR_ID
     5 CRISPR_start
     6 CRISPR_end
     7 Trusted
     8 Repeat_subtype
     9 Repeat_type_probability
    10 Consensus_repeat
    11 N_repeats
    12 Repeat_len
    13 Repeat_identity
    14 Spacer_len_avg
    15 Spacer_identity
    16 Spacer_len_sem
    17 Operon_ID
    18 Distance_to_CRISPR
    19 Cas_start
    20 Cas_end
    21 Best_type
    22 Best_score
    23 Genes
    24 Complete_interference
    25 Complete_adaptation
    26 Strand_cas
    27 N_genes
    28 Gene_lengths_aa
```

Note that while these are taken from CCTyper, some names have been adjusted
or added. These per-batch tables are then concatenated into one overall
results table: `results/crisprs-primary.tsv`.

### 1.3 Use BED files to extract sequences

Using the BED files generated in the step above, we use
[`seqkit`](https://bioinf.shenwei.me/seqkit/usage/#subseq)
to extract the loci of interest from the original assemblies based on the
contig IDs and positions.
For this, we use a custom script:
[`bin/extract_crispr-cas_from_fasta.sh`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/extract_crispr-cas_from_fasta.sh)
.
This script extracts the regions from start to stop, with (up to)
5000bp up- and downstream of the region. These are then written to FASTA files
in the subdirectory `results/crispr_fasta/`

``` bash
$ ls -sh results/crispr_fasta/atb.assembly.incr_release.202408.batch.38/
total 188K
 48K CRISPR-Cas.fasta  140K CRISPR-Cas-with_flanks.fasta
```

### 1.4 Concatenate spacers

All spacers detected during the primary screening step of CRISPRscape
are concatenated in one file with `cat`.

``` bash
$ ls -sh results/spacers-primary.fasta 
360K results/spacers-primary.fasta

$ head results/spacers-primary.fasta 
>SAMN39693127.contig00003_1:1
TGCGTGGGCGACCTTTTGTATAAACAAGTT
>SAMN39693127.contig00003_1:2
TGTATCCCTTACCACTACTTTGCCATATTG
>SAMN39693127.contig00003_1:3
CTCTTACAAGCTCTATTGCAATAATTTTAA
>SAMN39693127.contig00003_1:4
TCTAAAAAGGAAATTTTTTAAAGCTGACTT
>SAMN39693127.contig00003_1:5
CAATGTTTTGTCAAGTTTCAAGCGAAGGCG
```

## Cluster spacers: CCTyper

All CRISPR spacer sequences identified by CCTyper are clustered with
CD-HIT to determine the number of unique spacers. This happens in the rule
`cluster_unique_spacers`, and generates:

``` bash
ls -sh results/cluster-cctyper/spacers-clustered*
948K results/cluster-cctyper/spacers-clustered
 20M results/cluster-cctyper/spacers-clustered.clstr
```

The file without extension is the FASTA file that CD-HIT outputs,
`spacers-clustered.clstr` is the list of clusters as determined
by CD-HIT (in FASTA-like format).

Also see the user guide on
[bioinformatics.org](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide)
for details on CD-HIT's output files and extra tools.

### Spacer cluster table

To facilitate further analyses on the spacer level, we have included a script
to summarise the spacer cluster information in one table. This scripts,
[`bin/make_cluster_table.py`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/bin/make_cluster_table.py),
parses the output of CD-HIT to create a tab-separated table
`results/spacers-primary.tsv` that looks like:

``` bash
$ head results/spacers-primary.tsv 
Spacer_ID Sequence Length Cluster Spacer_base_ID Spacer_position
SAMN39693127.contig00003_1:1 TGCGTGGGCGACCTTTTGTATAAACAAGTT 30 1129 SAMN39693127.contig00003_1 1
SAMN39693127.contig00003_1:2 TGTATCCCTTACCACTACTTTGCCATATTG 30 1026 SAMN39693127.contig00003_1 2
SAMN39693127.contig00003_1:3 CTCTTACAAGCTCTATTGCAATAATTTTAA 30 983 SAMN39693127.contig00003_1 3
SAMN39693127.contig00003_1:4 TCTAAAAAGGAAATTTTTTAAAGCTGACTT 30 1156 SAMN39693127.contig00003_1 4
SAMN39693127.contig00003_1:5 CAATGTTTTGTCAAGTTTCAAGCGAAGGCG 30 904 SAMN39693127.contig00003_1 5
SAMN39693127.contig00003_1:6 AAACCTAAATCGTTATTCTCACCGCTCTTT 30 1051 SAMN39693127.contig00003_1 6
SAMN39693127.contig00003_1:7 AATAATGGCTAAATATTTCATGAGAATGGA 30 840 SAMN39693127.contig00003_1 7
SAMN39693127.contig00007_2:1 CGTAGAATACGATGAAAG 18 1635 SAMN39693127.contig00007_2 1
SAMN39693127.contig00010_3:1 TTAAAGTCCCTAAAAATAGGCATTTTTTGCTTTAGGTTAGGTAC 44 781 SAMN39693127.contig00010_3 1
```

The results per cluster are summarised in a separate table  that lists
the number of sequences per cluster, the most common spacer sequence,
the shortest and the longest sequence: `results/spacer_clusters-primary.tsv`.

!!! info "Final output:"

    `results/crisprs-primary.tsv`
    `results/spacers-primary.fasta`
    `results/spacers-primary-deduplicated_and_counted.fasta`
    `results/spacers-primary.tsv`
    `results/spacer_clusters-primary.tsv`
    

## 2. CRISPR-Cas refinement

### 2.1 Evaluation by CRISPRidentify

See the original documentation on
[CRISPRidentify's output](https://github.com/BackofenLab/CRISPRidentify#output-files).

Complete CRISPR-Cas arrays detected by CCTyper are passed on to
CRISPRidentify for further evaluation and refinement. This generates
per array output files, and summaries per batch:

```bash
$ ls -F results/crispridentify/atb.assembly.incr_release.202408.batch.38/SAMEA112866769-contig00010_2/
Alternative_Candidates.txt  Bona-Fide_Candidates.txt  gff_result/               Possible_Candidates.txt            Repeats.fasta  Summary.csv
Arrays.fasta                combined.gff              Low_Score_Candidates.txt  Possible_Discarded_Candidates.txt  Spacers.fasta

$ ls -F results/crispridentify/atb.assembly.incr_release.202408.batch.38/
Complete_array_dataset.fasta   Complete_repeat_dataset.fasta  SAMEA112866769-contig00010_2/
Complete_Cassette_summary.csv  Complete_spacer_dataset.fasta  SAMEA112867010-contig00031_10-11/
Complete_Cas_summary.csv       Complete_summary.csv
```

The 'Complete_summary.csv' files are combined for all batches into one
large table: `results/crisprs-final.tsv`.

```bash
$ head -3 results/crisprs-final.tsv 
Name    Global ID       ID      Region index    Start   End     Length  Consensus repeat        Repeat Length   Average Spacer Length   Number of spacers       Strand  Category        Score   batch
SAMN39693127-contig00003_1      1       1       1       9409    9906    498     ATTTTACCATAAAGAAATTTAAAAAGGGACTAAAAC    36      30      7       Forward (Orientation was not computed)  Bona-fide       0.9462372113536278  atb.assembly.incr_release.202408.batch.23
SAMN39693157-contig00006_6      2       1       1       9410    10039   630     ATTTTACCATAAAGAAATTTAAAAAGGGACTAAAAC    36      30      9       Forward (Orientation was not computed)  Bona-fide       0.9462372113536278  atb.assembly.incr_release.202408.batch.23
```

Also, spacers are collected in one overall file and clustered into unique
sequences, generating: `results/spacers-final.fasta`,
`results/spacers-final-deduplicated_and_counted.fasta`,
`results/spacers-final.tsv`, and `results/spacer_clusters-final.tsv`.

```bash
$ ls -1sh results/spacer*final*
 76K results/spacer_clusters-final.tsv
 52K results/spacers-final-deduplicated_and_counted.fasta
108K results/spacers-final.fasta
176K results/spacers-final.tsv

$ head -2 results/spacer*final*
==> results/spacer_clusters-final.tsv <==
Cluster Number_of_spacers       Most_common_sequence    Most_common_length      Most_common_number      Shortest_sequence       Shortest_length Longest_sequence        Longest_length
0       3       TCTAGCTTGCCTAGTATAGTTACACTTCCAC 31      3       TCTAGCTTGCCTAGTATAGTTACACTTCCAC 31      TCTAGCTTGCCTAGTATAGTTACACTTCCAC 31

==> results/spacers-final-deduplicated_and_counted.fasta <==
>15.AATAATGGCTAAATATTTCATGAGAATGGA
AATAATGGCTAAATATTTCATGAGAATGGA

==> results/spacers-final.fasta <==
>SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906_spacer_1
TGCGTGGGCGACCTTTTGTATAAACAAGTT

==> results/spacers-final.tsv <==
Spacer_ID       Sequence        Length  Cluster Spacer_base_ID  Spacer_position
SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906_spacer_1        TGCGTGGGCGACCTTTTGTATAAACAAGTT  30      252     SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906 1
```

!!! info "Final output:"

    `results/crisprs-final.tsv`
    `results/spacers-final.fasta`
    `results/spacers-final-deduplicated_and_counted.fasta`
    `results/spacers-final.tsv`
    `results/spacer_clusters-final.tsv`

## Target prediction

To predict the targets of CRISPR spacers, CRISPRscape maps the unique
spacer sequences back to the input genomes while ignoring the loci from
which the spacers derive. This is done with KMA, generating output like
described in [its documentation](https://github.com/genomicepidemiology/kma#result-explanation).

From that, we collect a list of spacers and the genome contig to which it maps:
`results/kma/CRISPR-alignments.tsv`.

```bash
$ ls -shF results/kma
total 564K
276K CRISPR_alignment.tsv   44K CRISPR.aln  184K CRISPR.frag.gz   20K CRISPR.fsa   36K CRISPR.res  4.0K spacer_DB/

$ head -5 results/kma/CRISPR_alignment.tsv
spacer  genome
SAMN39693127-contig00003_9409_9906      SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906_spacer_1 SAMEA115501831.contig00001
SAMN39693127-contig00003_9409_9906      SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906_spacer_1 SAMEA115501785.contig00006
SAMN39693127-contig00003_9409_9906      SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906_spacer_1 SAMEA115501713.contig00005
SAMN39693127-contig00003_9409_9906      SAMN39693127-contig00003_1_-_CRISPR_1_9409_9906_spacer_1 SAMEA115501708.contig00005
```

This contig information can be matched with the output of
[geNomad and Jaeger](#virus-and-plasmid-predictions) to predict if spacers
match bacteriophage or plasmid sequences.

## Multilocus Sequence Typing

The workflow includes classical MLST with pyMLST and the pubMLST database.
The output is a simple table with two columns containing the genome's name
and its assigned sequence type (if any):

``` bash
$ head results/mlst_table.tsv
Genome          ST
SAMN10081961    48
SAMN09605354    12411
SAMN09343354    8
SAMN09405115    8
SAMN09489478    9111
```

!!! info "Final output: `results/mlst_table.tsv`"

## Antiviral defence systems

For identifying diverse antiviral defence systems in bacteria, we use
[PADLOC](https://padloc.otago.ac.nz/padloc/).
The developers have explained its output in table format as described
[on GitHub](https://github.com/padlocbio/padloc#interpreting-output).

We concatenate output from all genome batches in one file:
`results/padloc_table.csv`. A comma-separated table with 19
columns:

```bash
$ ls -sh results/padloc_table.csv 
9.0M results/padloc_table.csv

$ head -1 results/padloc_table.csv | tr "," "\n" | nl
     1  system.number
     2  seqid
     3  system
     4  target.name
     5  hmm.accession
     6  hmm.name
     7  protein.name
     8  full.seq.E.value
     9  domain.iE.value
    10  target.coverage
    11  hmm.coverage
    12  start
    13  end
    14  strand
    15  target.description
    16  relative.position
    17  contig.end
    18  all.domains
    19  best.hits
```

!!! info "Final output: `data/processed/padloc_table.csv`"

## Virus and plasmid predictions

For each contig, we collect chromosome/plasmid/virus predictions derived
from geNomad and Jaeger.

For geNomad, we combine the prediction scores with available plasmid or
virus information into one tab-separatede dataframe:
`results/genomad_predictions.tsv`.

``` bash
$ ls -sh results/genomad_predictions.tsv 
2.0M results/genomad_predictions.tsv

$ head -5 results/genomad_predictions.tsv 
genome  batch   contig  chromosome_score        plasmid_score   virus_score     plasmid_topology        plasmid_genes   conjugation_genes       amr_genes       virus_topology  virus_genes     virus_taxonomy  genomad_prediction
SAMN39693127    atb.assembly.incr_release.202408.batch.23       SAMN39693127.contig00001        0.6773  0.2585  0.0641  NA      NA      NA      NA      NA      NA      NA      chromosome
SAMN39693127    atb.assembly.incr_release.202408.batch.23       SAMN39693127.contig00002        0.7124  0.23    0.0576  NA      NA      NA      NA      NA      NA      NA      chromosome
SAMN39693127    atb.assembly.incr_release.202408.batch.23       SAMN39693127.contig00003        0.7214  0.2219  0.0567  NA      NA      NA      NA      NA      NA      NA      chromosome
SAMN39693127    atb.assembly.incr_release.202408.batch.23       SAMN39693127.contig00004        0.7716  0.1783  0.0502  NA      NA      NA      NA      NA      NA      NA      chromosome
```

If geNomad did not assign virus or plasmid information, we classify the contig
as 'chromosome'. This is listed in the final column, 'genomad_prediction'.

!!! info "Final output: `results/genomad_predictions.tsv`"

We also run Jaeger as additional virus check. The output is described in its
[own documentation](https://github.com/Yasas1994/Jaeger#what-is-in-the-output).

```bash
$ ls -sh results/jaeger_predictions.tsv 
440K results/jaeger_predictions.tsv

$ head -5 results/jaeger_predictions.tsv 
accession_id    contig  length  prediction      reliability_score       prophage_contam
SAMN39830054    SAMN39830054.contig00012        10686   non-phage       0.892   TRUE
SAMN39706706    SAMN39706706.contig00037        3106    non-phage       0.957   FALSE
SAMN39847329    SAMN39847329.contig00016        41162   non-phage       0.932   TRUE
SAMN39693127    SAMN39693127.contig00017        6561    non-phage       0.882   TRUE
```

!!! info "Final output: `results/jaeger_predictions.tsv`"

## Genome dereplication

As described in the [extra functions](extra_functions.md#genome-dereplication),
CRISPRscape looks for identical genomes and provides a 'dereplication table'
that may be used to discard duplicates for statistical analyses. The output
is written to `results/dereplication_table.tsv`.
This table lists for each genome the corresponding cluster representative,
which may be used to flag and remove duplicates.

```bash
$ head -4 results/dereplication_table.tsv 
batch   genome  secondary_cluster       representative  score   cluster_method  comparison_algorithm    threshold       primary_cluster
atb.assembly.incr_release.202408.batch.23       SAMN39748578.fa 1_1     SAMN39748578.fa 1.6179841263202184      complete        fastANI 9.9999999999989e-5      1
atb.assembly.incr_release.202408.batch.23       SAMN40530429.fa 1_2     SAMN40530429.fa 1.6180244997520763      complete        fastANI 9.9999999999989e-5      1
atb.assembly.incr_release.202408.batch.23       SAMN39746878.fa 1_3     SAMN39746878.fa 1.6179033569300467      complete        fastANI 9.9999999999989e-5      1
```

!!! info "Final output: `results/dereplication_table.tsv`"
