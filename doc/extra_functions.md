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
