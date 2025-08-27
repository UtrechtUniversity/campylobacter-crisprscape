# Development notes

During the development of this workflow, we have made a number of decisions
and ran into unexpected situations when using certain tools. This document
describes the rationale of developing the workflow in its current form and
sheds light on _why we use the tools and parameters we use_.

## Downloading databases

(Rationale for including in the `Snakefile` or providing separate scripts
that need to be run before running the workflow with `snakemake`.
_This item needs to be updated!_)

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

## Modified CCTyper installation

We found that the default conda-based installation of CRISPRCasTyper would
create a non-functional program. We have reported this on
[GitHub](https://github.com/Russel88/CRISPRCasTyper/issues/55) and provide
a modified software environment configuration file (`envs/cctyper.yaml`)
that we found to work.

## Trying to be economical with disk space use

Downloading thousands of genomes from ATB and unpacking them into plain FASTA
files needs quite a lot of space. We considered gzipping each fasta file
separately to save disk space while still allowing most tools to work with
them. However, we found CCTyper would not take gzipped input files, so we
left the input files uncompressed.

## Optimising runtime on genome batches

When analysing batches of ATB genomes, we found two ways to process genomes
in parallel and speed up the whole process. One is to use GNU parallel to
read multiple input genomes at the same time. The other is to concatenate
all input files and use multiple threads for processing the whole batch
simultaneously. We hypothesised this would generate the exact same output,
but found that running CCTyper in these two different ways yielded slightly
different output files. We found that the file-by-file approach in one
particular batch would return 7 hits that were not found from the concatenated
file, while the concatenated approach yielded 1 hit that was not found in
the separate genomes. We did not pursue this discrepancy further and
decided to stick with the GNU parallel processing method.

## Unexpected behaviour of tools and databases

### Spacepharer

We found that Spacepharer wants to run within the same folder that you
designate the tmpfolder and where the created its databases.
This needs to be taken into account when writing the commands for the workflow.

Also, following the easy-predict workflow is not recommended as created
`.fasta` files are inconsistently recognized as actual fasta files.

Spacepharer databases are best created using the example names
"`querysetDB`" and "`targetsetDB`" as other names such as "`spacersetDB`"
causes weird errors.

### Phagescope

The Phagescope database says that it can filter genomes based on criteria,
but actually downloading these fastas is impossible due to an error.
Additionally, `wget` and `curl` do not properly download the databases in a way that
spacepharer can identify, requiring a manual upload.

### CRISPRidentify

The installation of CRISPRidentify was, unfortunately, not as straightforward
as running

``` bash
mamba env create -f CRISPRidentify.yml
```

or similar. We found that the YAML file provided on the CRISPRidentify GitHub
page cannot be solved using the recommended 'strict channel priority' setting.
We could only get it to work with the flexible or disabled channel priority.
We have adapted to YAML so that it can be solved using strict priority mode,
but this does make strand prediction in CRISPRidentify non-functional.
