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

1. [x] Verify CRISPR screening functionality (CCTyper + CRISPRidentify);
 generate CRISPR-Cas array table and CRISPR spacer table

2. [ ] Complete documentation for 'basic' functions:

    - [ ] User manual

    - [ ] Step-by-step process description

    - [ ] Output file descriptions

3. [ ] Clean-up code

    - [ ] Remove outdated/unnecessary steps

    - [ ] Apply standardised code formatting

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

- [ ] Rewrite 'Problems encountered' into a rationale for our tool selection (as separate document)

- [ ] Write detailed and technical step-by-step description of the workflow

    - [ ] While reviewing the workflow, remove unnecessary pieces and clean-up where possible

- [x] Setup MkDocs-powered documentation (at least locally, integrate with GitHub pages later)

(_Note to self: Add to this list when other ideas come to mind!_)

## Workflow description

### Input files to prepare

- _Campylobacter jujuni_ and _coli_ genomes from [AllTheBacteria](https://allthebacteria.readthedocs.io/en/latest/)

As of 2024-09-19 this includes 129,080 genomes!
(Up from 104,146 before the incremental update.
That means there are 24,934 extra genomes now.)

_Note: AllTheBacteria has its own [quality criteria](https://allthebacteria.readthedocs.io/en/latest/sample_metadata.html#high-quality-dataset) for inclusion._
_This includes:_

- \>=99% species abundance (practically pure)
- \>= 90% completeness (CheckM2)
- \<= 5% contaminated (CheckM2)
- total length between 100 kbp and 15 Mbp
- \<= 2,000 contigs
- \>= N50 2,000

### Download genomes

This repository includes scripts to automatically download genome files and metadata from ATB.
These can be run as follows:

```bash
git clone https://github.com/UtrechtUniversity/campylobacter-crisprscape.git
cd campylobacter-crisprscape
bash bin/prepare_genomes.sh
```

By default this downloads high-quality _Campylobacter jejuni_ and _C. coli_ genomes from the incremental update.
This [`prepare_genomes.sh`](bin/prepare_genomes.sh) script links to other scripts and has to be run from the 'base' folder as shown above.
The script itself contains a general description of how it works and to use it.
In short, it:

 1. Dowloads the metadata from AllTheBacteria (see [this script](bin/download_ATB_metadata.sh)).
By default, it downloads to the `data/ATB/` subdirectory and a different directory
can be provided as command-line argument.

 2. Extract sample accession IDs of the species of interest, as defined in
[`config/species_of_interest.txt`](config/species_of_interest.txt).
Edit this file if you want to run this workflow for different species!

 3. Look up the metadata of the species of interest by filtering the ENA metadata file.

 4. Find the batches in AllTheBacteria that contain the species of interest.

 5. Download the genome sequences of the species of interest (i.e., the batches identified in step 4).

 6. Remove other species from the downloaded batches.
(Batches may contain a mix of different species.)

As of February 2025, AllTheBacteria consists of an original set of genomes and
an incremental update. The `prepare_genomes.sh` script can download either part,
or all of the genomes at once using command-line options 'all', 'original', or
'update' (default: update).

The genomes are downloaded to the `data/tmp/ATB/` subdirectory. This is also the
default input directory for the analysis workflow.

### Download Databases

#### [geNomad](https://portal.nersc.gov/genomad/index.html)

This workflow uses geNomad to predict whether genomic contigs derive from
chromosomal DNA, plasmids or viruses. This tool uses both a neural network
classifier and a marker-based approach to calculate prediction scores.
For the marker-based method, it requires a database which can be downloaded
using the tool geNomad itself. If you have installed
[mamba](https://mamba.readthedocs.io/en/latest/)
this can be done as follows:

``` bash
mamba env create -f envs/genomad.yaml
mamba activate genomad
genomad download-database data/
```

Note that this will create the subdirectory `data/genomad_db/`,
which is the default that is also defined in
[`config/parameters.yaml`](config/parameters.yaml).

The current version of the database, v1.7, uses 1.4GB disk space.

#### [SpacePHARER](https://github.com/soedinglab/spacepharer)

The bin folder also includes scripts to download and extract pre-selected
databases for use in Spacepharer.
These include [Phagescope](https://phagescope.deepomics.org/) for annotated
phage sequences and [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/)
for annotated plasmid sequences which have been chosen for their broad taxonomy.
By running:

```bash
bash bin/download_spacepharer_database.sh
```

Both databases are downloaded, extracted and then merged for use in Spacepharer.
If you wish to use a different database or add to them, see
[`doc/spacepharer.md`](https://utrechtuniversity.github.io/campylobacter-crisprscape/spacepharer.html)
for advice.

### dependency CRISPRidentify
CRISPRscape uses CRISPRidentify as a second pass over cctyper's output. However this program has no conda environment that contains the program itself and as of writing requires a forked version to function properly.
To install CRISPRidentify, run

```bash
git clone https://github.com/Necopy-byte/CRISPRidentify/tree/master bin/
```

### Analysis workflow

The analysis itself is recorded as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.
Its dependencies (bioinformatics tools) are handled by Snakemake using the conda
package manager, or rather its successor [mamba](https://mamba.readthedocs.io/en/latest/).
If you have not yet done so, please install mamba following the instructions
found here: <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>.

After installing mamba, snakemake can be installed using their instructions:
<https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation>
(Note: the workflow was tested with Snakemake version 8.20.3 and is expected
to work with any version since 5.)

When Snakemake has been set up, you can test if the workflow is ready to be run
(dry-run) with:

```bash
snakemake --profile config -n
```

If that returns no errors, run the workflow by removin the `-n` (dry-run) option:

```bash
snakemake --profile config
```

Note that the workflow is currently configured to run on the local machine
(not on a high-performance computing (HPC) cluster or grid) and uses a maximum
of 24 CPU threads. The number of threads to use can be configured in:
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

 5. Predicts the potential targets of spacers and whether they target chromosomal DNA (of input genomes), plasmid or viruses using [Spacepharer](https://github.com/fbaumdicker/SpacerPlacer) (version 1.0.1) and [kma](https://github.com/genomicepidemiology/kma).

Further steps are added to the workflow after testing!

### Suggestions of programs/analyses to test

1. Mash with CRISPR loci, and whole genomes

2. Map CRISPR spacers to all downloaded genomes (bowtie, and KMA?), metagenome assemblies, other databases?

3. Whole-genome MLST

4. SpacerPlacer (see input file format in <https://github.com/fbaumdicker/SpacerPlacer?tab=readme-ov-file#spacer_fasta-input-format>
 (also requires an extra conversion script?)

## Output files

Unticked boxes indicate that documentation has not been written yet.

- [ ] ENA metadata
  - Cleaned-up and filtered metadata of included genomes

- [ ] Species classifications
  - Taxonomic classification by Sylph, as collected from AllTheBacteria

- [ ] Contig chromosome/plasmid/virus predictions

- [ ] CRISPR-Cas overview table
  - Output from CCTyper, collected and combined in one CSV file
  - Combine with CRISPRidentify, create filtered crispr by adding Cas and orientation data onto crispridentify csv 

- [ ] CRISPR spacer table

- [ ] MLST
  - Sequence Types (ST) of all included genomes

- [ ] List of spacer-putative targets
  - Output from mapping unique spacers to possible targets seperated by plasmid or phage and merged with database metadata

- [ ] List of anti-phage systems per genome
  - Output from PADLOC, combined in single CSV file

- [ ] Genome comparison all-vs-all
  - by whole-genome MLST, average nucleotide identity (ANI) or similar(?)

## Problems encountered

2024-09-12:

- CCTyper does not work from conda installation
    (<https://github.com/Russel88/CRISPRCasTyper/issues/55>)

- CCTyper cannot handle gzipped fasta files as input,
    so the input files need to be uncompressed

- RFPlasmid does not work from local install (file not found, permission denied)
     or conda installation (`rfplasmid --initialize` fails to download database files)

- 'hybrid' RFPlasmid install works: activate conda environment and then run the
     shared executable using the absolute path
     (/mnt/DGK_KLIF/data/rfplasmidweb/pip_package/package_files/RFPlasmid/rfplasmid.py)
     but it seems not to work when only one genome is present in the input directory.

2024-12-10:

- CCTyper returns different output when running single genomes as
    compared to running one concatenated fasta file with contigs of
    many genomes. (See `data/tmp/cctyper/batch_22/CRISPR_Cas-batch_22.tab`
    vs `data/tmp/cctyper/test/batch_22/CRISPR_Cas.tab` - difference is 1KB.)
    The separate method returns 7 contigs that the concatenated did not find,
    and the concatenated method found 1 contig that the separate did not find.
    These may be false positives (how do you check?), but for now I'm sticking
    with the separate method.

2025-02-11:

- Spacepharer wants to run locally within the same folder that you designate the tmpfolder and where the created
   databases are located. Following the easy-predict workflow is not recommended as created .fasta files are
   inconsistently recognized as actual fasta files.
- Phagescope database says that it can filter genomes based on criteria, but actually downloading these fastas is
   impossible through an error. Additionally wget and curl do not properly download the databases in a way that
   spacepharer can identify, requiring a manual upload.

2025-02-21:

- Spacepharer databases are best created using the example names "querysetDB" and "targetsetDB" as other names such as "spacersetDB" causes weird errors.

2025-06-27:

- The yml file provided by CRISPRidentify is only solveable using flexible or disabled channel priority in conda. As of now an adjusted yml file is used that is solveable in strict priority mode, but does make strand prediction in CRISPRidentify non-functional.

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
