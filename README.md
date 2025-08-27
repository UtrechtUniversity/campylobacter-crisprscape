# _Campylobacter_ CRISPRscape

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip) [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![pages-build-deployment](https://github.com/UtrechtUniversity/campylobacter-crisprscape/actions/workflows/pages/pages-build-deployment/badge.svg)](https://utrechtuniversity.github.io/campylobacter-crisprscape/)

Automated bioinformatics workflow for the identification and characterisation
of Clustered Regularly Interspaced Short Palindromic Repeat (CRISPR) arrays
from bacterial genomes, with a focus on _Campylobacter coli_ and _C. jejuni_.

For more detailed documentation, please look at the
[project website](https://utrechtuniversity.github.io/campylobacter-crisprscape/index.html).

## Index

1. [Workflow description](#workflow-description)
   - [download genomes](#download-genomes)
   - [download databases](#download-databases)
   - [analysis workflow](#analysis-workflow)
   - [future suggestions](#suggestions-of-programsanalyses-to-test)
2. [Output files](#output-files)
3. [Problems encountered](#problems-encountered)
4. [Project (file) organisation](#project-organisation)
5. [Licence](#licence)
6. [Citation](#citation)

---

## Release roadmap

- Version 0.1 updates:

    - Updates in documentation of implemented functions and output files
(e.g., expaning documentation on the use of CRISPRidentify)

    - Functions to concatenate and combine existing outputs

    - Code clean-up (moving long commands to separate scripts,
applying standardised formatting, remove unnecessary code)

- Versions 0.2 and further:

    - New functionality (e.g., all-vs-all genome comparisons)

## To do

- [x] Add CRISPRidentify to workflow

   - [x] And make sure it works on a clean install

- [x] Combine results from CCTyper with CRISPRidentify

- [ ] Make and/or correct scripts for combining results into 'Output files' (write to `data/processed/`)

   - [ ] Concatenate MLST results

   - [x] Enable spacer table creation script in Snakefile (add to `rule all`)

- [x] Collect and combine results from geNomad and Jaeger

- [x] Map spacers to genomes and phage/plasmid databases

- [x] Add PADLOC for identifying other anti-phage systems

- [ ] Write documentation for output files

- [x] Rewrite 'Problems encountered' into a rationale for our tool selection (as separate document)

- [x] Write detailed and technical step-by-step description of the workflow

    - [ ] While reviewing the workflow, remove unnecessary pieces and clean-up where possible

- [x] Setup MkDocs-powered documentation (at least locally, integrate with GitHub pages later)

(_Note to self: Add to this list when other ideas come to mind!_)

## User manual

Please consult the
[user manual](https://utrechtuniversity.github.io/campylobacter-crisprscape/manual.html)
on our documentation page, which includes a quick start guide as well as a
detailed step-by-step description.

## Workflow description

### Analysis workflow

The analysis is recorded as a
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.
Its dependencies (bioinformatics tools) are handled by Snakemake using the conda
package manager, or rather its successor
[mamba](https://mamba.readthedocs.io/en/latest/).
If you have not yet done so, please install mamba following the instructions
found here:
<https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>.

After installing mamba, snakemake can be installed using their instructions:
<https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation>
(Note: the workflow was tested with Snakemake version 8.20.3 and is expected
to work with any version since 5.)

When Snakemake has been set up, you can test if the workflow is ready to be run
(dry-run) with:

```bash
snakemake --profile config -n
```

If that returns no errors, run the workflow by removing the `-n` (dry-run) option:

```bash
snakemake --profile config
```

Note that the workflow is currently configured to run on the local machine
(not on a high-performance computing (HPC) cluster or grid) and uses a maximum
of 60 CPU threads. The number of threads to use can be configured in:
[`config/config.yaml`](config/config.yaml) (overall workflow) and
[`config/parameters.yaml`](config/parameters.yaml) (per step/tool).

In its current state, the workflow:

 1. Determines MultiLocus Sequence Types for _Campylobacter jejuni/coli_ using the public typing scheme from [pubMLST](https://pubmlst.org/organisms/campylobacter-jejunicoli) and [pyMLST](https://pymlst.readthedocs.io/) (version 2.2.1)

 2. Identifies CRISPR-Cas loci with [CCTyper](https://github.com/Russel88/CRISPRCasTyper) (version 1.8.0) and resulting loci are processed with [CRISPRidentify](https://github.com/BackofenLab/CRISPRidentify) ([forked](https://github.com/Necopy-byte/CRISPRidentify) from version 1.2.1) to reduce false positives.
    - this includes extra scripts to collect CRISRP-Cas information and extract sequences from the genome fasta files

 3. Collects all CRISPR spacers and creates clusters of identical spacers using
[CD-HIT-EST](https://sites.google.com/view/cd-hit) (version 4.8.1)

 4. Predicts whether contigs of the species of interest derive from chromosomal DNA,
plasmids or viruses using both [geNomad](https://portal.nersc.gov/genomad/index.html) (version 1.8.0)
and [Jaeger](https://github.com/Yasas1994/Jaeger) (version 1.1.26).

 5. Predicts the potential targets of spacers and whether they target
chromosomal DNA (of input genomes), plasmid or viruses using
[Spacepharer](https://github.com/soedinglab/spacepharer) (version 5.c2e680a)
and [kma](https://github.com/genomicepidemiology/kma) (version 1.5.0).

Further steps are added to the workflow after testing!

### Preparing input and databases

Before running the workflow, the user needs to prepare input genomes and
databases.
Also see the corresponding documentation pages for details:

- [Prepare genomes](https://utrechtuniversity.github.io/campylobacter-crisprscape/prepare_genomes.html)
- [Download databases](https://utrechtuniversity.github.io/campylobacter-crisprscape/manual.html#downloading-databases)

### Suggestions of programs/analyses to test

1. Mash with CRISPR loci, and whole genomes (compare all-vs-all)

2. Map CRISPR spacers to all downloaded genomes (bowtie, KMA, and Sassy?),
metagenome assemblies, other databases to predict targets

3. SpacerPlacer (see input file format in <https://github.com/fbaumdicker/SpacerPlacer?tab=readme-ov-file#spacer_fasta-input-format>
 (also requires an extra conversion script?)

## Output files

Ticked boxes indicate that
[documentation](https://utrechtuniversity.github.io/campylobacter-crisprscape/output_files.html)
is available.

- [x] ENA metadata
  - Cleaned-up and filtered metadata of included genomes

- [x] Species classifications
  - Taxonomic classification by Sylph, as collected from AllTheBacteria

- [ ] Contig chromosome/plasmid/virus predictions

- [x] CRISPR-Cas overview table
  - Output from CCTyper, collected and combined in one CSV file
  - Combine with CRISPRidentify, create filtered crispr by adding Cas and orientation data onto crispridentify csv

- [x] CRISPR spacer table

- [ ] MLST
  - Sequence Types (ST) of all included genomes

- [ ] List of spacer-putative targets
  - Output from mapping unique spacers to possible targets seperated by plasmid or phage and merged with database metadata

- [ ] List of anti-phage systems per genome
  - Output from PADLOC, combined in single CSV file

- [ ] Genome comparison all-vs-all
  - by whole-genome MLST, average nucleotide identity (ANI) or similar(?)

## Project organisation

```bash
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── Snakefile          <- Python-based workflow description
├── bin                <- Code and programs used in this project/experiment
├── config             <- Configuration of Snakemake workflow
├── data               <- All project data, divided in subfolders
│   ├── processed      <- Final data, used for visualisation (e.g. tables)
│   ├── raw            <- Raw data, original, should not be modified (e.g. fastq files)
│   └── tmp            <- Intermediate data, derived from the raw data, but not yet ready for visualisation
├── doc                <- Project documentation, notes and experiment records
├── envs               <- Conda environments necessary to run the project/experiment
├── log                <- Log files from programs
└── results            <- Figures or reports generated from processed data
```

---

## Licence

This project is licensed under the terms of the [New BSD licence](LICENSE).

---

## Citation

Please cite this project as described in the [citation file](CITATION.cff).
