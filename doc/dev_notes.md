# Development notes

During the development of this workflow, we have made a number of decisions
and ran into unexpected situations when using certain tools. This document
describes the rationale of developing the workflow in its current form and
sheds light on _why we use the tools and parameters we use_.

## Downloading databases

(Rationale for including in the `Snakefile` or providing separate scripts
that need to be run before running the workflow with `snakemake`.)

Databases are downloaded for:

- ATB input genomes (script `bin/prepare_genomes.sh`)

- Multilocus sequence typing (Snakemake rule `download_mlst_database`)

- PADLOC anti-phage defence system screening (Snakemake rule `download_padloc_database`)

- geNomad plasmid/phage predictions (not included, undocumented)

- SpacePHARER CRISPR target prediction (script `download_spacepharer_database.sh`)

## Use of GNU parallel over Snakemake

Snakemake is very good at running processes in parallel and making
efficient use of computational resources. However, the file-based nature
combined with the necessity to generate a directed acyclic graph (DAG)
before running the actual workflow may lead to excessive waiting times
when you have tens of thousands of input files. Therefore, we decided
to use the batches provided by
[AllTheBacteria](https://allthebacteria.readthedocs.io/en/latest/)
as input units and use a separate program,
[GNU parallel](https://www.gnu.org/software/parallel/), to handle
the files within these batches. This solution is taken from this online post
by John van Dam: <https://stackoverflow.com/a/55007872>. As described there,
this makes the whole DAG simpler, and the processing faster, at the cost
of losing checks by Snakemake to see if all expected output files are
generated.

## Using other inputs than ATB

This workflow was designed to work with bacterial genomes from
[AllTheBacteria](https://allthebacteria.readthedocs.io/en/latest/) (ATB)
as input and has provided scripts to automatically download genomes
of the bacterial species of interest.
However, the processing steps should work equally well on other genomes in
FASTA format. Currently, the easiest way to 'trick' the workflow into using
custom input is to either modify the input directory line in
`config/parameters.yaml`:

``` yaml
input_directory: "data/tmp/ATB/"
```

to wherever you stored your files, or create this directory and copy
your genomes to there. Furthermore, Snakemake expects the input to be stored
in batches, using a directory that is named 'batch_*', where the asterisk
can match anything. So for example, you can make a directory:
`data/tmp/ATB/batch_custom`, put your FASTA files in there and Snakemake
should pick them up and have them analysed.

Also note that downstream processing, i.e., statistical analyses and
visualisation, make use of the metadata provided with ATB. This would have
to be adjusted to the custom data as well.
