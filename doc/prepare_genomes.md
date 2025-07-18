# Prepare genomes :material-database-arrow-down

A step-by-step explanation of the `bin/prepare_genomes.sh` script that
downloads genomes of interest from [AllTheBacteria](https://allthebacteria.readthedocs.io/en/latest/)
along with the relevant metadata.

## 1. Download ATB metadata

The script starts by calling a separate script: `bin/download_ATB_metadata.sh`.
This script downloads:

1. ENA metadata (standard metadata from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home))
    - for this, it downloads both the '20240801' and '0.2.20240606' files

2. Sample list (accession IDs for all samples included in ATB)

3. Sylph output (taxonomic classification based on the
[Genome Taxonomy Database](https://gtdb.ecogenomic.org/))

4. Assembly statistics (total length, number of contigs, N50, and more)

5. Assembly quality assessment ([CheckM2](https://github.com/chklovski/CheckM2))

6. Species calls (an interpretation of the Sylph and CheckM2 output)
    - three columns: accession ID, species name and high-quality (T/F)

7. A list of all files
    - this is actually a table with accession IDs, species name, which
    ATB batch file they are part of, batch MD5 checksum and file size

### Output files

By default, these are all downloaded to the directory `data/ATB/`.
You then get:

``` bash
ena_metadata.0.2.20240606.tsv.gz
ena_metadata.20240801.tsv.gz
sample_list.txt.gz
sylph.tsv.gz
assembly-stats.tsv.gz
checkm2.tsv.gz
species_calls.tsv.gz
file_list.all.20240805.tsv.gz
```

## 2. Look up accession IDs of species of interest

Using the file `config/species_of_interest.txt`, which contains one species
name per line, the script looks up all matches in the `species_calls.tsv.gz`
file from ATB and filters the high-quality assemblies.
The accession IDs are stored in a separate file:
`data/ATB/all_samples_of_interest.txt`.

As an extra, the script reads the total number of selected genomes,
which is printed to the command-line (stdout).
Also, the number of genomes per species is collected and stored as:
`data/ATB/number_of_genomes_per_species.txt`.

Note: the system works with the [GTDB taxonomy](https://gtdb.ecogenomic.org/tree?r=d__Bacteria)
and GNU `grep`. This means it can only find names as defined by the GTDB
and can use parts of the names too. E.g., _Campylobacter jejuni_ is known
in GTDB as 'Campylobacter_D jejuni' and by using this as search term,
the script also matches 'Campylobacter_D jejuni_A' or any other suffix.

## 3. Filter metadata

The complete metadata file is big (677MiB compressed, 3,112,707 lines).
To make this a bit easier to work with, we're extracting only the lines
referring to the species of interest and store this as a separate file:
`data/ATB/enametadata.20240801-filtered.tsv.gz`.

## 4. Find batches that contain species of interest

ATB has created batches of genomes to use clever compression and
significantly reduce file sizes. To find which batches contain
the species of interest, the script:

1. reads sample accession IDs from `data/ATB/all_samples_of_interest.txt`

2. uses the accession IDs to filter matching lines in `data/ATB/file_list.all.20240805.tsv.gz`

3. extracts the column containing file names, download URLs and MD5 checksums

4. deduplicates them, keeping one copy of each batch containing species of interest

5. and stores this in a separate file: `data/ATB/batches_to_download.tsv`

## 5. Download the genomes

The actual downloading of the genomes happens here! :material-file-download:
The script calls yet another separate script: `bin/download_genomes.sh`

`download_genomes.sh` reads the files that need to be downloaded
from `data/ATB/batches_to_download.tsv` and downloads them one by one
to `data/tmp/ATB/`.
It checks file integrity with the MD5 checksum and deletes corrupted
files. (It does not retry downloading automatically.)

The files (as batches) are downloaded as XZ archive, which are
extracted to subdirectories named after the batch number.
This yields FASTA files in a directory called `data/tmp/ATB/batch_[number]`.

## 6. Remove other species

Some batches may contain more than only the species of interest.
In those cases, FASTA files that contain genomes from other species are
deleted.
To do this, the final script `bin/remove_other_species.sh` is called.
This script compares the accession IDs (samples) of the species of interest
with the accession IDs present in the downloaded batches to
compose a list of samples with other species (`data/tmp/other_species-samples.txt`).
It counts how many samples have other species and then the corresponding
FASTA files are removed.
For curious users, the other species names and numbers are stored as:
`data/tmp/other_species_calls.tsv` and `data/tmp/other_species_numbers.txt`.

## General remarks

- The whole process is set up as Bash script, and uses GNU tools such
as `grep` and `wget`. This should work on most Unix-like systems.

- The script can use a command-line option to donwload genomes from the
complete ATB, only the first version, or the incremental update.
For this, it uses options 'all', 'original' or 'update', respectively.
By default, the script uses 'update', which has the smallest file sizes.
The option can be provided as:

```bash
bash bin/prepare_genomes.sh all
```

This option is most relevant to `bin/download_genomes.sh`.
[(Step 5)](#5-download-the-genomes)

- The script reads whether files exist already. If they are already there,
they will not be downloaded again.
