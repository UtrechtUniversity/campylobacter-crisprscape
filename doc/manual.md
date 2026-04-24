# CRISPRscape user manual

## Quick start :material-run-fast:

Install dependencies:
[git](https://git-scm.com/downloads/),
[mamba](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)
and
[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation).

Download the repository (including submodules):

``` bash
git clone --recurse-submodules https://github.com/UtrechtUniversity/campylobacter-crisprscape.git
```

Move into the downloaded directory:

``` bash
cd campylobacter-crisprscape
```

Download genomes from AllTheBacteria

``` bash
bash bin/prepare_genomes.sh
```

_Optional:_ Do a dry-run to check if Snakemake can find all the right input
and output files:

``` bash
snakemake --profile config -n
```

Run the actual analysis workflow:

``` bash
snakemake --profile config
```

_For a more detailed explanation and information for adjusting parameters,_
_please see below._

## 1. Before you start :octicons-info-16:

The CRISPRscape workflow relies on two main tools for managing the workflow
and installing software:
[Snakemake](https://snakemake.readthedocs.io/en/stable/)
and
[mamba](https://mamba.readthedocs.io/en/latest/index.html).
Furthermore, since the project is hosted on GitHub, we expect you to use
[git](https://git-scm.com/).

CRISPRscape is designed to work with the
[AllTheBacteria](https://allthebacteria.readthedocs.io/en/latest/)
resource to download all high-quality genomes of a given species.
(It has been tested on _Campylobacter coli_ and _C. jejuni_ together.)

### Estimated disk use :material-harddisk:

!!! warning

    CRISPRscape requires downloading multiple databases.
    Prepare to use hundreds of GBs!

- [AllTheBacteria metadata](allthebacteria.md): ~1.5GB
    - AllTheBacteria genomes: depends on species
    (e.g., 209GB for ~130,000 _Campylobacter_ genomes)
    - AllTheBacteria genome annotations (Bakta): depends on species
    (e.g., 531GB for ~130,000 _Campylobacter_ genomes)

- Databases of [SpacePHARER](#spacepharer):
    - [PLSDB](spacepharer.md#plsdb-2024_05_31_v2): 14GB initial download
    \+ 68GB after indexing = ~80GB

    - [Phagescope](spacepharer.md#phagescope): 41GB initial download
    \+ 155GB after indexing = ~200GB

- [geNomad database](#genomad): 1.4GB

### Download and install software :material-download:

Before you begin, you need to install:
(follow these links to find installation instructions)

1. [git](https://git-scm.com/downloads/)

2. [mamba](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)

3. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation)

We recommend Snakemake is installed via mamba. This is also the default
and linked above.

When you have these tools installed, you can download CRISPRscape:

``` bash
git clone --recurse-submodules https://github.com/UtrechtUniversity/campylobacter-crisprscape.git
```

(Note: the `--recurse-submodules` option is necessary to also automatically
download [CRISPRidentify](https://github.com/necopy-byte/CRISPRidentify),
which is one of the two CRISPR-Cas screening tools included.)

Move your current working directory into this newly downloaded one
to get started!

``` bash
cd campylobacter-crisprscape
```

(You may of course rename this directory if you want to. In that case,
use the new name instead of `campylobacter-crisprscape`.)

### Tunable parameters :material-tune:

CRISPRscape includes some options that can be modified by the user.
These are stored in YAML and TXT files under the `config` directory.

#### Species of interest :fontawesome-solid-binoculars:

The species of interest can be modified to your liking by changing the file
[`config/species_of_interest.txt`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/main/config/species_of_interest.txt)

This file simply lists species names, one per line. The default is:

``` bash
Campylobacter_D jejuni
Campylobacter_D coli
```

By changing the species name, one can adjust the species that can be
automatically downloaded from AllTheBacteria (ATB). This requires the
taxonomic names as defined in [GTDB](https://gtdb.ecogenomic.org/tree).

This also affects multilocus sequence typing (MLST): CRISPRscape includes
automated MLST, which requires downloading the proper marker gene database.
This information is stored under
[`config/parameters.yaml`](https://github.com/UtrechtUniversity/campylobacter-crisprscape/blob/8f66f87d44ef6a8761e372d67074227f0b64c026/config/parameters.yaml#L15)
.
For finding valid species names, please consult
[pyMLST](https://pymlst.readthedocs.io/en/latest/documentation/clamlst/initialise.html#import-from-pubmlst).

#### Technical parameters :material-tools:

Then, there are some technical parameters that you can adjust to fit your
system. These range from the input directory in which your genomes are stored
(default: `resources/ATB/`) and the location of databases to the number of CPU
threads to use.

Please open `config/config.yaml` and `config/parameters.yaml` to review
the default parameters and adjust where needed to make it work for your
system.

The default CPU settings are:

- Use a maximum total of 60 CPU threads

- Use 20 CPU threads for most compute-intensive tasks

### Download input genomes :material-cloud-download:

CRISPRscape includes a convenient script to autmatically download genomes
of interest from ATB. When the desired species name is saved in
`config/species_of_interest.txt`, you can start downloading with:

``` bash
bash bin/prepare_genomes.sh
```

Here, you may optionally add a command-line parameter to tell which part
of ATB to look into: the complete database (`all` = largest), only the
original version (`original`) or the incremental update (`update` = smallest;
default). For example:

``` bash
bash bin/prepare_genomes.sh --part all
```

One may also set the output directory using the `--directory` or `-d` option:

```bash
bash bin/prepare_genomes.sh --part all --directory my/genomes
```

It also has a built-in help function to show the supported options:

```bash
$ bash bin/prepare_genomes.sh --help
Prepare genomes from AllTheBacteria for use with CRISPRscape

Syntax: prepare_genomes.sh -p [part] -d [directory] [-h]
Options:
-p/--part      Select which part of AllTheBacteria to download
               for the selected species ('all', 'original', or
               'update', default=update)
-d/--directory Directory in which to download the files
               (default=resources/ATB/)
-h/--help      Print this help message
```

### Downloading databases :material-cloud-download:

Reference databases required by the different tools are automatically
downloaded with the Snakemake workflow. No user action required.
Below is a summary of the databases that are used.

#### [PADLOC](https://github.com/padlocbio/padloc)

CRISPRscape uses PADLOC: Prokaryotic Antiviral Defence LOCator to screen
genomes for the presence of different antiviral defence systems, including
_Cas_ genes. PADLOC uses a database of genes and HMM classifications to
predict presence of a large number of
[different defence systems](https://github.com/padlocbio/padloc-db/blob/master/system_info.md).

The current version, v2.0.0, uses 954MB disk space.

#### [geNomad](https://portal.nersc.gov/genomad/index.html)

This workflow uses geNomad to predict whether genomic contigs derive from
chromosomal DNA, plasmids or viruses. This tool uses both a neural network
classifier and a marker-based approach to calculate prediction scores.
For the marker-based method, it requires a database which can be downloaded
using the tool geNomad itself.

The current version of the database, v1.7, uses 1.4GB disk space.

#### [SpacePHARER](https://github.com/soedinglab/spacepharer)

For details of SpacePHARER, please see
[the corresponding page](spacepharer.md).

## 2. Running the workflow :material-run

The workflow is fully automated and should complete with one command.
For details on what happens under the hood, see the tab 'Workflow details'.

One can do a 'dry-run' to test if all preparations have been satisfied:

``` bash
snakemake --profile config -n
```

To run the actual workflow:

``` bash
snakemake --profile config
```

## 3. Interpreting results :material-magnify:

After running the workflow, the user is presented with a number of output files.
These are described in detail under the tab '[Output files](output_files.md)'.
