# Extra functions included in the workflow

Functions that are not directly related to CRISPR-Cas.
That is, other bacterial genomics analyses that serve to better
characterise and understand the input data.

## Multilocus Sequence Typing

To facilitate comparison with other studies and estimate the genetic diversity
within the input dataset, the workflow can assign Multilocus Sequence Types
(MLST) using
[pyMLST](https://pymlst.readthedocs.io/en/latest/index.html),
which can download databases from [pubmlst.org](https://pubmlst.org/).

Each genome is assigned its closest sequence type (ST) when possible.

## Identification of antiviral defence systems

To identify different antiviral defence systems, we have included
[PADLOC](https://padloc.otago.ac.nz/padloc/).
PADLOC uses a custom database to screen genomes for the presence of various
systems. This information is written to a
[CSV file](https://github.com/padlocbio/padloc#output).
We use this information to assess the prevalence of defence systems and
calculate possible correlations with CRISPR-Cas.

## Plasmid and virus prediction of input genomes

We have included two state-of-the-art tools for predicting the origin of
contigs: chromosome, plasmid or virus (phage). These are
[geNomad](https://portal.nersc.gov/genomad/) and
[Jaeger](https://github.com/Yasas1994/Jaeger).
geNomad compares contigs to a marker database and classifies them using a
neural network, aggregates these results and calculates probability scores
for chromosome, plasmid and virus (these three sum up to 1).
Jaeger uses deep-learning to classify sequences as viral or not. This can
be used to identify 'free' phages as well as integrated prophages.

We map the identified CRISPR spacers back to the input genomes and use
these plasmid/phage predictions to estimate the targets of the CRISPRs.

## Mapping spacers back to input genomes

CRISPRscape maps the identified spacers back to the genomes, ignoring
the genomes from which they derive, to predict what their targets
(protospacers) may be. Spacers are mapped back to genomes using the
fast k-mer based tool
[KMA](https://github.com/genomicepidemiology/kma) (version 1.5.0).

By using the KMA flag `-hmm`, the output files (`*.frag.gz`) become easier
to interpret: it add columns with (1) target name, (2) start and (3) stop
positions. With this option enabled, KMA "uses a HMM to assign template".

## Genome dereplication

When calculating the prevalence of CRISPR-Cas for a species of interest, it
is worth considering duplicated genomes. Multiple copies of the same origin
may be present in ATB, which would affect statistics such as prevalence.
(Simply put, if there are 10 genomes, 9 of which are identical and have
a CRISPR-Cas locus, is the prevalence 9/10 or really 1/2? I would argue
the latter is more appropriate.)

To determine which genomes are identical, we use the tool
[dRep](https://drep.readthedocs.io/en/latest/index.html).
dRep incorporates quality parameters as predicted by
[CheckM](http://ecogenomics.github.io/CheckM/) as well as simple assembly
statistics such as length with genome comparisons based on Average Nucleotide
Identity (ANI). Given an identity threshold and a clustering algorithm,
dRep will clusters all genomes and select the best representative per
cluster based on CheckM and assembly statistics. We opted for a 99.99%
threshold to dereplicate genomes.

Below is the dRep command as implemented in Snakemake, followed by an
explanation of the command-line options:

```bash
find {input.batch} -name "*.fa" -print > {output.genome_list}

dRep dereplicate {params.prefix} -g {output.genome_list}\
 -sa 0.9999 -p {threads} --ignoreGenomeQuality\
 --genomeInfo {input.checkm}\
 --clusterAlg complete -pa 0.99 -nc 0.5 > {log} 2>&1
```

First, we use GNU find to search for the fasta files per batch, to
provide a list of input files to dRep.

Then the workflow runs `dRep dereplicate`, writing output to the
given 'prefix'. The freshly generated genome list is given as input (`-g`).
`-sa 0.9999` specifies the identity threshold for **s**econdary
**A**NI clustering, which is set to 99.99%.
`-p {threads}` sets the number of CPU threads to use.
`--ignoreGenomeQuality` really tells dRep not to calculate CheckM scores
again. These are provided using `--genomeInfo {input.checkm}`. This
uses the completeness and contamination scores calculated by CheckM2,
and downloaded from ATB. (_Note: the developers of ATB calculated these_
_scores for us, so we download that rather than running the same program_
_again!_) The clustering algorithm is set to 'complete' using the
`--clusterAlg` option. (For details on how these clustering algorithms
work, I recommend looking at the figures in
[this guide](https://stackabuse.com/hierarchical-clustering-with-python-and-scikit-learn/#linkagemethods).)
The primary clustering threshold is set to 99% using `-pa 0.99`,
which is already a strict cutoff splitting the genome batch in multiple
clusters. This prevents a compute intensive secondary clustering on large
numbers of genomes. Finally, the `-nc 0.5` option sets the minimum coverage
of genome comparisons to 50%. In other words, only genomes that cover one
another for at least half of their length are compared - which should work
fine, given that we include only high-quality assemblies from ATB.
