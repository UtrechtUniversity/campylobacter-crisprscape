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

it is recommended to comment out this rule and instead merge to metadata yourself, as it assumes a certain metadata structure that can change between databases.

### [Phagescope](https://phagescope.deepomics.org/)

Phagescope is a large database containing a large set of bacteriophages combined from several sources and further annotated with information such as taxonomy. Specifically this database was chosen as it has a large variety of host species, sufficient annotation for CRISPRscape's goal and is relatively large compared to other databases. Furthermore, depending on further research, the extra annotations are extensive and could be used for further research.

As of writing Phagescope contains 873.718 phage sequences from 4.723 host species.

### [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) 2024_05_31_v2

PLSDB is an up to date database containing plasmids taken from NCBI and INSDC that have been scrubbed of chromosomal sequences. Furthermore PLSDB has extensive annotation including taxonomy and antimicrobial resistance genes.

This database was chosen for its wide variety of host species, origin location and being of sufficient size. Furthermore, as of writing it is still getting updates.

The version used by CRISPRscape (PLSDB 2024_05_31_v2) contains 72.360 plasmid sequences from 9436 different host species.

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
