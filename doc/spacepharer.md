# Spacepharer documentation

## [Spacepharer](https://github.com/soedinglab/spacepharer)

As described by the github, Spacepharer is a program that uses the search method of MMSEQ2 to query spacer sequences against a database with the intent to find host-phage matches. Normally Spacer sequences have quite short sequences which confound more conventional methods such as BLASTN and therefore require stringent filtering methods. However phage databases are frequently incomplete and phages are generally highly mutable, necessitating that a certain amount of mismatches are allowed. This is what spacepharer overcomes by prioritising open reading frames which bypass the need for highly similar matches.

Spacepharer was therefore chosen as the method in CRISPRscape to identify spacer to phage matches compared to conventional methods such as BLASTN. With similar rationale, spacepharer was also chosen to match spacers with plasmids as they were expected to have similar problems such as virality causing high mutability and consisting mostly of coding sequences.

## Databases

Spacepharer requires the use of databases to match against. For this, CRISPRscape comes with the script `download_spacepharer_database.sh` which download both Phagescope and PLSDB as databases for phage and plasmid matching respectively and prepare them for use in Spacepharer. It is possible to use other databases but this requires some changes in both `Snakefile` and `parameters.yaml` to account for them.

Specifically `Snakefile` would require changes to rules:

`spacepharer_phage_setup`:

params: DB="path/to/database.fasta"

`spacepharer_phage_setup`:

params: DB="path/to/database.fasta"

`create_spacepharer_table`:

It is recommended to comment out this rule and instead merge the metadata yourself, as it assumes a certain metadata structure that can change between databases.

### [Phagescope](https://phagescope.deepomics.org/)

Phagescope is a large database containing a large set of bacteriophages combined from several sources and further annotated with information such as taxonomy. Specifically this database was chosen as it has a large variety of host species, sufficient annotation for CRISPRscape's goal and is relatively large compared to other databases. Furthermore, depending on further research, the extra annotations are extensive and could be used for further research.

As of writing Phagescope contains 873.718 phage sequences from 4.723 host species.

The full download comprises about 41GB, of which `merged_metadata.tsv` (127MB)
and the fasta files are the most important.

- `merged_metadata.tsv` lists for each phage its genome length, GC-content,
taxonomy, completeness (e.g., 'High-quality'), host species, lifestyle
(e.g, 'virulent'), and which cluster and subcluster it belongs to.

From these, a custom database is build, which uses **279GB**.

``` bash
$ ls -sh data/raw/phagescope/
total 41G
1.4G CHVD.fasta                    11M GPD                         2.6G MGV.tar.gz
5.7M chvd_phage_meta_data.tsv     5.1G GPD.fasta                   336K PhagesDB
400M CHVD.tar.gz                   18M gpd_phage_meta_data.tsv     224M PhagesDB.fasta
 28K DDBJ                         1.6G GPD.tar.gz                  524K phagesdb_phage_meta_data.tsv
 12M DDBJ.fasta                   4.0M GVD                          68M PhagesDB.tar.gz
 40K ddbj_phage_meta_data.tsv     560M GVD.fasta                   356K RefSeq
3.6M DDBJ.tar.gz                  4.9M gvd_phage_meta_data.tsv     318M RefSeq.fasta
 12K EMBL                         172M GVD.tar.gz                  572K refseq_phage_meta_data.tsv
8.5M EMBL.fasta                   241M IGVD.fasta                   97M RefSeq.tar.gz
 20K embl_phage_meta_data.tsv     1.2M igvd_phage_meta_data.tsv     98M STV.fasta
2.6M EMBL.tar.gz                   70M IGVD.tar.gz                 536K stv_phage_meta_data.tsv
164K Genbank                      8.2G IMG_VR.fasta                 28M STV.tar.gz
147M Genbank.fasta                 31M img_vr_phage_meta_data.tsv  5.0M TemPhD
256K genbank_phage_meta_data.tsv  2.4G IMG_VR.tar.gz               2.5G TemPhD.fasta
 45M Genbank.tar.gz               127M merged_metadata.tsv         8.6M temphd_phage_meta_data.tsv
4.0G GOV2.fasta                    15M MGV                         770M TemPhD.tar.gz
 35M gov2_phage_meta_data.tsv     8.4G MGV.fasta
1.1G GOV2.tar.gz                   23M mgv_phage_meta_data.tsv
```

### [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) 2024_05_31_v2

PLSDB is an up to date database containing plasmids taken from NCBI and INSDC that have been scrubbed of chromosomal sequences. Furthermore PLSDB has extensive annotation including taxonomy and antimicrobial resistance genes.

This database was chosen for its wide variety of host species, origin location and being of sufficient size. Furthermore, as of writing it is still getting updates.

The version used by CRISPRscape (PLSDB 2024_05_31_v2) contains 72.360 plasmid sequences from 9436 different host species.

The full download uses 14GB of disk space, but the essential files for
CRISPRscape are smaller:

- `sequences.fasta`: 7GB
    - self-explanatory: fasta file with the plasmid sequences

- `nuccore.csv`: 17MB
    - Metadata file with unique ID, accession number, description, creation
     date, predicted completeness, 'Genome' (plasmid), length, and 14 more
     columns with a.o. source information

- `taxonomy.csv`: 2.3MB
    - Taxonomic tree for all included plasmids as comma-separated table,
    including both taxonomic names for each rank and taxID numbers.

These are used to build a custom database of **68GB**.

``` bash
$ ls -sh data/raw/PLSDB/
total 14G
 37M amr.tsv                         1.7M disease_terms.csv      16M plasmidfinder.csv      16K README.md
3.1M assembly.csv                    3.5G download_meta.tar.gz   14M plsdb_mashdb_sim.tsv  7.0G sequences.fasta
 16K biosample_attributes_plsdb.csv  8.0K ecopaths.csv          569M plsdb_sketch.msh      2.3M taxonomy.csv
 12M biosample.csv                   780K nucc_identical.csv    610M proteins.csv           23M typing.csv
2.5M changes.tsv                      17M nuccore.csv           1.7G proteins.fasta         75M typing_markers.csv
```

## Output

Spacepharers output consists of a singular .tsv file showing all matches to a phage_ID within the database provided. The expected output takes the form of:

```
#prok_acc  phage_acc   S_comb      num_hits
>spacer_acc      phage_acc   p_bh    spacer_start      spacer_end  phage_start phage_end   5'_PAM|3'_PAM    5'_PAM|3'_PAM(reverse strand)
```
where the # indicates the prokaryotic match. However due to the way that input is provided, this ends up being redundant. CRISPRscape therefore filters out the # sections retaining only the individual matches.

of the matches the columns consist of:
`spacer accession, phage accession, p besthit, spacer start and end, phage start and end, possible 5’ PAM|3’ PAM, possible 5’ PAM|3’ PAM on the reverse strand`

Lastly the matches are then matched with their respective metadata and added onto the table. Crisprscape will by default take the metadata from Phagescope and PLSDB and add the relevant columns.

For Phagescope these columns are: `Length, GC_content, taxonomy, completeness, host, lifestyle`

For PLSDB the column added is: `Host taxonomy`
