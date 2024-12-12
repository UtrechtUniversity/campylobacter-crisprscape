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

This repository includes scripts to automatically download genome files and metadata from ATB.
These can be run as follows:

```bash
git clone https://github.com/UtrechtUniversity/campylobacter-crisprscape.git
cd campylobacter-crisprscape
bash bin/prepare_genomes.sh
```

By default this downloads high-quality _Campylobacter jejuni_ and _C. coli_ genomes from the incremental update.
This `prepare_genomes.sh` script links to other scripts and has to be run from the 'base' folder as shown above.

_Note to self: write user manual with instructions how to edit the configuration!_

## List of programs/analyses to test

1. CCTyper on _Campylobacter_ genomes

    - use CRISPR spacer hits to do saturation experiment: how many new spacers are found when adding more genomes?
      (how likely is it that we cover all spacers with the current dataset?)

2. geNomad for identification of plasmids and phages(?)

  - also compare results to RFPlasmid and ProphET(?)

3. Mash with CRISPR loci

4. kma CRISPR spacers to metagenome assemblies?

5. annotate metagenome contig hits with CAT and geNomad

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
├── doc                <- Project documentation for readthedocs.io
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
