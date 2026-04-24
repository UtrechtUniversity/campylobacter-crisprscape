# _Campylobacter_ CRISPRscape

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![pages-build-deployment](https://github.com/UtrechtUniversity/campylobacter-crisprscape/actions/workflows/pages/pages-build-deployment/badge.svg)](https://utrechtuniversity.github.io/campylobacter-crisprscape/)  
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/UtrechtUniversity/campylobacter-crisprscape)

Automated bioinformatics workflow for the identification and characterisation
of Clustered Regularly Interspaced Short Palindromic Repeat (CRISPR) arrays
from bacterial genomes, with a focus on _Campylobacter coli_ and _C. jejuni_.

For more detailed documentation, please look at the
[project website](https://utrechtuniversity.github.io/campylobacter-crisprscape/index.html).
Here you also find the
[user manual](https://utrechtuniversity.github.io/campylobacter-crisprscape/manual.html),
which includes a quick start guide as well as a detailed step-by-step description.
Or look at the computer-generated Wiki with integrated chatbot assistant at
[DeepWiki](https://deepwiki.com/UtrechtUniversity/campylobacter-crisprscape)!

## Index

1. [Roadmap](#release-roadmap)
2. [Workflow description](#workflow-description)
3. [Project (file) organisation](#project-organisation)
4. [Licence](#licence)
5. [Citation](#citation)

---

## Release roadmap

- Future additions:
  - Extract orphan CRISPR arrays from CCTyper -> send to CRISPRidentify
  - CRISPR spacer target prediction
    - map to
      - masked ATB genomes (KMA)
      - PLSDB (SpacePHARER)
      - PhageScope (SpacePHARER)
      - VIRE (t.b.d.)
      - MEGAISurv metagenomes (t.b.d)
    - mini-benchmark different mapping algorithms?
      - Sassy
      - KMA
      - SpacePHARER
    - (where feasible) connect spacer hits with functional annotations
      - Bakta annotations from ATB are available!
  - Integrate downstream analyses with Snakemake?
    - run RMarkdown/Quarto notebooks automatically
  - Build a database like [this spacerdb](https://spacers.jgi.doe.gov/database/overview/)?
  - Create a command-line tool with [snaketool](https://github.com/beardymcjohnface/Snaketool/)
  - Separate 'main' workflow steps from 'optional/extra' steps?

## Workflow description

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

### Preparing input genomes

Before running the workflow, the user needs to prepare input genomes.
These may be downloaded from AllTheBacteria using a convenient script
included with CRISPRscape:

- [Prepare genomes](https://utrechtuniversity.github.io/campylobacter-crisprscape/prepare_genomes.html)

## Project organisation

```bash
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── bin                <- Code and programs used in this projec
├── config             <- Configuration of Snakemake workflow
├── doc                <- Project documentation, notes and experiment records
├── log                <- Log files from programs
├── resources          <- Folder for databases and input files
├── results            <- Workflow output, all results
└── workflow           <- Files describing the Snakemake workflow
    ├── Snakefile      <- Main workflow description file (Python-based)
    ├── envs           <- Software environments (conda) for dependencies
    ├── notebooks      <- Analysis notebooks for downstream statistics/visualisation
    ├── rules          <- Workflow modules with processing steps
    └── scripts        <- Scripts used in the workflow
```

---

## Licence

This project is licensed under the terms of the [New BSD licence](LICENSE).

---

## Citation

Please cite this project as described in the [citation file](CITATION.cff).
