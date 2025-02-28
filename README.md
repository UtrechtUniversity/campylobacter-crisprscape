# Test CRISPR workflow

_Date: 2024-09-12_  
_Author: Sam Nooij (s.nooij [at] uu.nl)_


Test a CRISPR analysis workflow for the PINNACLE project

---

## Input files to prepare:

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

### Analysis workflow

The analysis itself is recorded as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.
Its dependencies (bioinformatics tools) are handled by Snakemake using the conda
package manager, or rather its successor [mamba](https://mamba.readthedocs.io/en/latest/).
If you have not yet done so, please install mamba following the instructions
found here: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html.

After installing mamba, snakemake can be installed using their instructions:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation
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

 1. Identifies CRISPR-Cas loci with [CCTyper](https://github.com/Russel88/CRISPRCasTyper) (version 1.8.0)

 2. Collects all CRISPR spacers and creates clusters of identical spacers using
[CD-HIT-EST](https://sites.google.com/view/cd-hit) (version 4.8.1)

 3. Predicts whether contigs of the species of interest derive from chromosomal DNA,
plasmids or viruses using both [geNomad](https://portal.nersc.gov/genomad/index.html) (version 1.8.0)
and [Jaeger](https://github.com/Yasas1994/Jaeger) (version 1.1.26).

Further steps are added to the workflow after testing!

## Suggestions of programs/analyses to test

1. Mash with CRISPR loci, and whole genomes

2. Map (KMA?) CRISPR spacers to all downloaded genomes, metagenome assemblies, other databases?

3. Annotate metagenome contig hits with CAT

### Problems encountered:

2024-09-12:
  - CCTyper does not work from conda installation
    (https://github.com/Russel88/CRISPRCasTyper/issues/55)

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


### Files I want to collect and process

1. Complete CRISPR-Cas locus sequences as fasta
 (make a script to extract from contigs using positions listed by CCTyper)
    - also make one with only CRISPR spacers (concatenate from CCTyper(?))
    - and make one with flanking regions to infer and compare genomic location
    - collect all contigs with CRISPR hits to predict chromosome/plasmid/virus origin?
2. CRISPR input file for SpacerPlacer as defined in https://github.com/fbaumdicker/SpacerPlacer?tab=readme-ov-file#spacer_fasta-input-format
 (also requires an extra conversion script?)
3. Separate CRISPR spacer sequences (provided by CCTyper)
4. List of number of CRISPR arrays per genome


## Project organisation

```
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── Snakefile          <- Python-based workflow description
├── bin                <- Code and programs used in this project/experiment
├── config             <- Configuration of Snakemake workflow
├── data               <- All project data, divided in subfolders
│   ├── processed      <- Final data, used for visualisation (e.g. tables)
│   ├── raw            <- Raw data, original, should not be modified (e.g. fastq files)
│   └── tmp            <- Intermediate data, derived from the raw data, but not yet ready for visualisation
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

Please cite this project as described [here](CITATION.cff).
