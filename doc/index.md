# CRISPRscape documentation

[![Built with MkDocs](https://cdn.jsdelivr.net/npm/@intergrav/devins-badges@3/assets/cozy/built-with/mkdocs_vector.svg)](https://www.mkdocs.org/) + [![Using Material for MkDocs](https://cdn.jsdelivr.net/gh/Andre601/devins-badges@v3.x-mkdocs-material/assets/compact-minimal/built-with/mkdocs-material_vector.svg)](https://squidfunk.github.io/mkdocs-material)

!!! danger "Work in progress"

## How to use

To learn more about how to use the workflow, see the
[user manual](manual.md)

## What it does

The workflow works as follows:

1. [Download genome sequences](prepare_genomes.md)

2. [Screen for presence of CRISPR-Cas systems](CRISPR_screening.md)

3. [Refine identified CRISPR arrays](CRISPR_refinement.md)

4. [Cluster identified CRISPR spacers](clustering_spacers.md)

## What it creates

For an overview of the output files that are created in the process, please
consult [output files](output_files.md)

## Why it works the way it does

We have summarised some of the rationale of why we have implemented this
workflow the way it is in the [developer notes](dev_notes.md)