# AllTheBacteria details

Last updated: September 2024

## Metadata files and sizes

After downloading all metadata using the bundled script, you have:

``` bash
$ ls -sh data/ATB/
total 1.4G
679M ena_metadata.20240801.tsv.gz            74M assembly-stats.tsv.gz
 18M file_list.all.20240805.tsv.gz           63M checkm2.tsv.gz
5.5M sample_list.txt.gz                     371M ena_metadata.0.2.20240606.tsv.gz
17M species_calls.tsv.gz                    104M sylph.tsv.gz
```

- Assembly statistics `assembly-stats.tsv.gz` 74MB
  - Lists per sample accession the total length, number of contigs, N50 and
   more statistics of the assembly.

- CheckM2 results `checkm2.tsv.gz` 63MB
  - Lists per sample accession the results of CheckM2 including assembly
   completeness and contamination in percentages.

- ENA sample metadata `ena_metadata.20240801.tsv.gz` 677MB
  - Lists per sample all the metadata that have been deposited in the European
   Nucleotide Archive - this is a table with over 100 columns!

- File list `file_list.all.20240805.tsv.gz` 17MB
  - This file lists per sample accession the corresponding batch in which it is
   archived, with download URL, md5sum and file size of the batch archive.
    It can be used to identify which batch archives contain species of interest,
    e.g.

```bash
zgrep -f all_Campylobacter_samples.txt file_list.all.20240805.tsv.gz |\
 cut -f 3,4,5 | sort | uniq
```

- Sample list `sample_list.txt.gz` 5.4MB
  - This file simply lists all the sample accessions that are present in the dataset.

- Species calls `species_calls.tsv.gz` 17MB
  - Lists per sample accession the species identified and whether or not it is of
   high-quality (T/F).
    This file can be used to identify which sample accessions contain
     high-quality genomes of the species of interest,
    e.g.

```bash
zless species_calls.tsv.gz | grep "Campylobacter_D jejuni" |\
 grep -v "F" | cut -f 1
```

- Sylph results `sylph.tsv.gz` 103MB
  - Lists per sample accession the output from Sylph, including relative
   abundance (%) Average Nucleotide Identity score (%) and assigned species name.

## Note on file sizes per batch

- a batch of fasta files as .xz archive takes up 12-242MB of disk space (median ~30MB)
- the extracted fasta files take up 2.2GB (in the case of the 12MB .xz file)
- gzipping each fasta file (in fast mode) shrinks that to 716MB (~ 3x reduction)
