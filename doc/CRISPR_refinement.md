# CRISPR-Cas refinement workflow (simplified)

## 2. Assessing putative CRISPR-Cas loci :material-shield-check:

### Rationale of two-step approach

It has been shown that CCTyper, the primary method we use for CRISPR-Cas
identification, has a high false-positive rate. Or rather, that the tool
that it relies on (
[MinCED](https://github.com/ctSkennerton/minced)
, derived from 
[CRT](https://www.room220.com/crt/)
) has a high false-positive rate. Therefore, we include a second tool,
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

1. It has a more sophisticated method of identifying CRISPR repeats and
their correct start and stop positions (details in 
[Supplementary file 1](https://academic.oup.com/nar/article/49/4/e20/6027817?login=false#supplementary-data)
)

2. It uses a complex confidence scoring algorithm. This machine learning-based
approach should exclude nearly all false-positive CRISPR-Cas loci.

The downside, however, is that running CRISPRidentify on complete genome
assemblies takes a lot more time than CCTyper. Therefore, we came up with
this two-step approach.

## 2.1 Passing CCTyper output to CRISPRidentify

## 2.2 Re-evaluation of the CRISPR arrays

(Brief description of how CRISPRidentify works)

## 2.3 Calculating CRISPR confidence

(Brief description of what the machine learning thing does)

## 2.4 Combining the output with CCTyper's

### Output files generated in the process :file_folder: :material-file-table: :material-file-table:

Each step in the process generates a number of output files, which by default
are written to:

``` bash
data/
  tmp/
    crispridentify/                 # Here go overall files, such as 'all_spacers.fa'
      batch_[number]/               # Here is only one subfolder
        CRISPR_arrays-with_flanks/  # In here are subfolders for each CRISPR
                                    #  array identified with CCTyper.
          [CRISPR_ID]/              # Here are CRISPRidentify's output files
```

For more details on the output files, see [output](output_files.md).
