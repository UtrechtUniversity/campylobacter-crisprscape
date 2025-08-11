# Development notes

During the development of this workflow, we have made a number of decisions
and ran into unexpected situations when using certain tools. This document
describes the rationale of developing the workflow in its current form and
sheds light on _why we use the tools and parameters we use_.

## Downloading databases

(Rationale for including in the `Snakefile` or providing separate scripts
that need to be run before running the workflow with `snakemake`.)

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
