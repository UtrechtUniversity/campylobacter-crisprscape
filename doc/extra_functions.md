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
