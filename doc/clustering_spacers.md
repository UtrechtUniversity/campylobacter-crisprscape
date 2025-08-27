# Clustering spacers

Among the CRISPR spacers identified are likely to be many duplicates.
To estimate the number of unique spacers, we use
[CD-HIT-EST](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide.pdf)
(page 11)
to cluster spacer sequences and establish clusters of identical spacer
sequences.

The command that we use to cluster as written in the `Snakefile`,
rule `cluster_unique_spacers`, is:

``` bash
cd-hit-est -c 1 -n 8 -r 1 -g 1 -AS 0 -sf 1 -d 0 -T {threads}\
 -i {input} -o {output.spacers}
```

The options used do the following:

- `-c 1` sets the cluster threshold (as fraction) to 1, meaning all sequences in
the cluster share 100% identity. (default: 0.9)

- `-n 8` sets the word size to 8 (8-10 is recommended for thresholds 0.9-1).
This is a very technical thing and was found not to have a large effect on
runtime and end result.

- `-r 1` includes reverse complements (0 = no, 1 = yes). This means spacer
sequences from either DNA strand may be included in a cluster, rather than
being identified as separate sequences.

- `-g 1` enables the slow but accurate mode (default: 0). With this number
of sequences, this is still very fast and it is nice to be accurate when
possible.

- `-AS 0` sets the alingment coverage of the shorter sequences as number of
positions that may be missed. 0 means the complete shorter sequence must
match the longer sequence against which it is compared.

- `-sf 1` sorts the output fasta file by decreasing cluster size (default: 0)

- `-d 0` sets the length of sequence names as displayed in the output
`.clstr` file. In fact, 0 disables the limit, allowing the full length
of sequence IDs to be printed in the output. (default: 20)

- `-T {threads}` sets the number of CPU threads to use, where {threads}
is a variable controlled by Snakemake.
(default: 1. Set to 0 to use all available CPUs)

- `-i` and `-o` allow you to specify input and output files, respectively!

CD-HIT-EST is used to cluster all CRISPR spacer sequences as reported by
CCTyper, and later also those reported (or rather, modified) by CRISRPidentify.
This allows the user to tell how many unique spacer sequences were detected.

The output files are briefly described in the
[output files](output_files.md#cluster-spacers-cctyper)
page. For more details on CD-HIT's output, and for information on
alignment coverage control options (figure on page 9), see
[CD-HIT's documentation](http://www.bioinformatics.org/cd-hit/cd-hit-user-guide.pdf)
on bioinformatics.org.
