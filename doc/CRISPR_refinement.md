# CRISPR-Cas refinement workflow (simplified)

## 2. Assessing putative CRISPR-Cas loci :material-shield-check:

### Rationale of two-step approach

It has been shown that CCTyper, the primary method we use for CRISPR-Cas
identification, may have a high false-positive rate. Or rather, that the tool
that it relies on (
[MinCED](https://github.com/ctSkennerton/minced)
, derived from
[CRT](https://www.room220.com/crt/)
) has a high false-positive rate.
Furthermore, we have found that with highly similar CRISPR arrays in
different genomes, the start and stop positions of repeats and spacers may
shift by one or more positions when comparing these arrays.
Therefore, we include a second tool,
[CRISPRidentify](https://github.com/BackofenLab/CRISPRidentify)
to evaluate the loci identified by CCTyper.
For a detailed description of how CRISPRidentify works, we refer you to the
[publication](https://doi.org/10.1093/nar/gkaa1158).
Also note the related publication of
[CRISPRloci](https://doi.org/10.1093/nar/gkab456) for information on the
different modules of which this tool or suite consists.

!!! question "Word of caution regarding CRISPRidentify's superior accuracy"

    An important comment has been posted below the publication of
    CRISPRidentify, criticising their claims of superiority over
    previous tools.
    [Link to comment](https://academic.oup.com/nar/article/49/4/e20/6027817?login=false#usercomments)

    In short, results should be checked by an expert. There is no CRISPR-Cas
    tool that outclasses the other available tools in all aspects.

In brief, CRISPRidentify has two crucial advantages compared to CCTyper:

1. It has a more sophisticated method of identifying the correct start and
stop positions of CRISPR repeats (details in
[Supplementary file 1](https://academic.oup.com/nar/article/49/4/e20/6027817?login=false#supplementary-data)
)

2. It uses a complex confidence scoring algorithm. This machine learning-based
approach should exclude nearly all false-positive CRISPR-Cas loci.

The downside, however, is that running CRISPRidentify on complete genome
assemblies takes a lot more time than CCTyper. Therefore, we came up with
this two-step approach.

## 2.1 Passing CCTyper output to CRISPRidentify

## 2.2 Re-evaluation of the CRISPR arrays

Using the arrays parsed from CCTyper, CRISPRidentify has two seperate steps to
come to its final conclusion: Candidate generation and Candidate evaluation.
In Candidate generation, CRISPRidentify uses Vmatch to find putative repeat pairs,
which by default need to be between 21 and 55 nucleotides long and 18-78
nucleotides apart. This process is relatively sensitive and usually more than
one repeat candidate is generated. All the repeat candidates are then aligned.
This alignment is created from a maximum element and a minimum element. The
maximum element is the largest repeat string generated from the most common
nucleotides in each base of all candidates. The minimum element is generated
from the most common substring of all repeats, this also by definition has 100%
identity as a substring of the maximum element.

Every possible repeat is then generated between the maximum and minimum element
and put alongside the matches found by Vmatch and has duplicates filtered out,
forming the set of repeat candidates. This set of candidates is then even further
extended by omitting up to 3 nucleotides on each side of the repeats, generating
an additional 15 candidates per repeat. In essence, many possible repeats are
considered for analysis.

CRISPR array candidates are created by string searching the different repeats
in the provided sequence and attempts to minimise the number of editing
operations needed for the consensus repeat, while still allowing mutations
in the repeats to be detected.

## 2.3 Calculating CRISPR confidence

After all CRISPR array candidates are generated, they are evaluated by
CRISPRidentify's internal scoring system. This scoring system was created
by considering 13 features that can predict array viability in multiple ways.
The 13 features are listed in
[Supplementary file 1](https://academic.oup.com/nar/article/49/4/e20/6027817?login=false#supplementary-data),
table S2 (page 20).
Performing feature subset selection on all combinations of these 13 features,
three models containing 8, 9 and 10 of the 13 features achieved similar accuracy.
By default, CRISPRidentify uses the average of these three models to score the
candidate arrays.

The scoring is divided into three possible categories.
0-0.4 are low scoring candidates which are unlikely to be CRISPR.
0.4-0.75 are possible candidate CRISPR arrays.
0.75-1.0 are Bona-Fide CRISPR arrays and are very likely to be valid CRISPR arrays.
In cases where CRISPR array candidates are overlapping but both Bona-Fide,
the lower scoring arrays are instead put into the alternative candidate category.

## 2.4 Combining the output with CCTyper's

### Output files generated in the process :file_folder: :material-file-table: :material-file-table:

Each step in the process generates a number of output files, which by default
are written to:

```bash
results/
  crispridentify/
    [batch]/
      Complete_array_dataset.fasta        # Sequences of detected CRISPR arrays
      Complete_Cassette_summary.csv       # Summary of CRISPR-Cas cassettes detected by CRISPRidentify
                                          # (disabled in CRISPRscape)
      Complete_Cas_summary.csv            # Summary of detected Cas genes (disabled in CRISPRscape)
      Complete_repeat_dataset.fasta       # Sequences of detected CRISPR repeats
      Complete_spacer_dataset.fasta       # Sequences of detected CRISPR spacers
      Complete_summary.csv                # Table with CRISPR-Cas array information
      [array_1]/
        Alternative_Candidates.txt        # CRISPR detection report of CRISPR arrays
                                          #  with certainty scores
        Arrays.fasta                      # Sequences of detected CRISPR array
        Bona-Fide_Candidates.txt          # List of high-scoring CRISPR-Cas arrays
        combined.gff                      # Genome annotation (GFF) file with CRISPR-Cas
        gff_result/
        Low_Score_Candidates.txt          # List of low-scoring CRISPR-Cas arrays
        Possible_Candidates.txt           # List of intermediate-scoring CRISPR-Cas arrays
        Possible_Discarded_Candidates.txt # List of discarded CRISPR-Cas candidates
        Repeats.fasta                     # Sequences of detected CRISPR repeats
        Spacers.fasta                     # Sequences of detected CRISPR spacers
        Summary.csv                       # Summary table per array
      [array_2]/
        ...
      [array_n]/
        ...
```

For more details on the output files, see [output](output_files.md).

#### Next steps

&rarr; [Cluster spacers](clustering_spacers.md)
