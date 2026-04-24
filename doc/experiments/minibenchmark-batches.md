# Minibenchmark batched analysis

Run CCTyper on four batches and compare runtimes between 'parallel' method
and with concatenated input files.

## Times based on parallel

- fast = atb.assembly.r0.2.batch.103 (1h40m)
- 3 x slow =
  - atb.assembly.incr_release.202408.batch.24 (21h28m)
  - atb.assembly.r0.2.batch.40 (23h7m)
  - atb.assembly.r0.2.batch.51 (23h50m)

## Times with concatenated input

- fast = atb.assembly.r0.2.batch.103 (35m?)
- 3 x slow =
  - atb.assembly.incr_release.202408.batch.24 (21h26m?? - calculated from output files, not by `time` --> 10h30m if inferred using batch.40's and 103's times!)
  - atb.assembly.r0.2.batch.40 (10h32m)
  - atb.assembly.r0.2.batch.51 (12h36m)

Based on these numbers, concatenated input makes CCTyper **2-3 times** faster.
(2.86 times for smaller batches, 1.89 times for the large batch.)

## Times with concatenated input and --prodigal meta

- atb.assembly.r0.2.batch.103 (39m)
- atb.assembly.incr_release.202408.batch.24 (748m = 12h28m)
- atb.assembly.r0.2.batch.40 (774m = 12h54m)
- atb.assembly.r0.2.batch.51 (727m = 12h07m)

These times are similar to the above, time speedup is now:

2.56 - 1.72 - 1.79 - 1.97

So still roughly **twice as fast**.

### Try and compare the output

General findings

```bash
wc -l results/cctyper/atb.assembly.incr_release.202408.batch.24/CRISPR_Cas-atb.assembly.incr_release.202408.batch.24.tab results/cctyper/atb.assembly.r0.2.batch.103/CRISPR_Cas-atb.assembly.r0.2.batch.103.tab results/cctyper/atb.assembly.r0.2.batch.40/CRISPR_Cas-atb.assembly.r0.2.batch.40.tab results/cctyper/atb.assembly.r0.2.batch.51/CRISPR_Cas-atb.assembly.r0.2.batch.51.tab

  2693 results/cctyper/atb.assembly.incr_release.202408.batch.24/CRISPR_Cas-atb.assembly.incr_release.202408.batch.24.tab
   213 results/cctyper/atb.assembly.r0.2.batch.103/CRISPR_Cas-atb.assembly.r0.2.batch.103.tab
  1228 results/cctyper/atb.assembly.r0.2.batch.40/CRISPR_Cas-atb.assembly.r0.2.batch.40.tab
  2542 results/cctyper/atb.assembly.r0.2.batch.51/CRISPR_Cas-atb.assembly.r0.2.batch.51.tab
```

```bash
wc -l */CRISPR_Cas.tab
   2652 atb.assembly.incr_release.202408.batch.24/CRISPR_Cas.tab
    213 atb.assembly.r0.2.batch.103/CRISPR_Cas.tab
   1223 atb.assembly.r0.2.batch.40/CRISPR_Cas.tab
   2508 atb.assembly.r0.2.batch.51/CRISPR_Cas.tab

   2640 meta-atb.assembly.incr_release.202408.batch.24/CRISPR_Cas.tab
    213 meta-atb.assembly.r0.2.batch.103/CRISPR_Cas.tab
   1225 meta-atb.assembly.r0.2.batch.40/CRISPR_Cas.tab
   2522 meta-atb.assembly.r0.2.batch.51/CRISPR_Cas.tab
```

Using concatenated input files yields the same number or slightly fewer hits
compared to using separate files per genome (e.g. 213 for atb.assembly.r0.2.batch.103
in both cases, against 2542 vs 2508 = 34 (~1%) for atb.assembly.r0.2.batch.51).

Using the option `--prodigal meta` may results in slightly different numbers of hits
compared to using default settings with concatenated output.
(Note 24 has 12 fewer hits with meta, while 51 has 14 hits more.)
Perhaps this is partly random?

Comparing CRISPR IDs with diff is complicated by the fact that they use different
indexes. E.g. the 'original' has SAMEA95883418.contig00004@3 and the 'concatenated input'
has SAMEA95883418.contig00004@2.

```bash
diff <(cut -f 1,3,15 results/cctyper/atb.assembly.r0.2.batch.51/crisprs_all-atb.assembly.r0.2.batch.51.tab | sort) <(cut -f 1,3,15 ../crisprscape_integration_test/results/test_batches/atb.assembly.r0.2.batch.51/crisprs_all.tab | sort)
```

Indicates very subtle differences: 2 nt in start position,
or same position but different prediction confidence.
If the positions differ, it appears to be a consistent start position +2 for the
concatenated input approach, while the end position remains the same.

It appears that often if both approaches found the same CRISPR with the same
start position, the concatenated input method reports a lower prediction confidence,
e.g., 0.996 vs 0.519 for 'SAMEA112259549.contig00006    41248'.
This is, however, not always the case. Confidence may be (near) identical
or higher for the separate genomes approach...

#### Conclusion

I think it is a reasonable compromise to call these results interchangeable -
differences may mainly be within the 'false positives' -
and rely on CRISPRidentify to get rid of less reliable hits.

#### Detail on time consumpion

**Note: the two main time consuming steps within CCTyper are:**

- Prodigal (~20% of time)
- HMMer (~75% of time)

both of which may be significantly faster if Pyrodigal and PyHMMER were used!

- [Pyrodigal benchmark](https://pyrodigal.readthedocs.io/en/stable/guide/benchmarks.html#single-run)
- [pyHMMER benchmark](https://pyhmmer.readthedocs.io/en/stable/guide/benchmarks.html#v0-11-1-2025-05-08)

The remaining time is:

- Subtyping Cas genes (>1%)
- Predicting arrays with minced (~1%)
- Predicting operons with BLAST (~1%)
- Predicting subtype based on repeats (0%, instant)
- Connecting Cas genes to CRISPR arrays (0%, seconds))
- Plotting loci (<=2%)
- File cleanup (0%)

So further optimisation by swapping MinCED for Diced would be much less impactful and only save less than 1% of the total time.

## Also test with PADLOC

Other processing steps also depend on parallel now, but may be faster and
more convenient if used with concatenated fasta files as input.

### Times with parallel

- fast = atb.assembly.r0.2.batch.103 (31m)
- 3 x slow =
  - atb.assembly.incr_release.202408.batch.24 (8h)
  - atb.assembly.r0.2.batch.40 (8h12m)
  - atb.assembly.r0.2.batch.51 (8h6m)

### Times concatenated input

- fast = atb.assembly.r0.2.batch.103 (68m)
- 3 x slow =
  - atb.assembly.incr_release.202408.batch.24 (11h25m)
  - atb.assembly.r0.2.batch.40 ()
  - atb.assembly.r0.2.batch.51 ()

PADLOC does not seem to like running on such big files!
And is faster on separate genomes.

## And with Jaeger

### Jaeger times with parallel

- fast = atb.assembly.r0.2.batch.103 (17m22s)
- 3 x slow =
  - atb.assembly.incr_release.202408.batch.24 (5h24m)
  - atb.assembly.r0.2.batch.40 (7h43m)
  - atb.assembly.r0.2.batch.51 (7h12m)

### Jaeger times concatenated input

Jaeger appears crazy fast with these concatenated input files...
I should have done this way earlier!

- fast = atb.assembly.r0.2.batch.103 (6m45)
- 3 x slow =
  - atb.assembly.incr_release.202408.batch.24 (1h37m)
  - atb.assembly.r0.2.batch.40 (1h40m)
  - atb.assembly.r0.2.batch.51 (1h53m)

Jaeger works **2.5-4.5 times faster** using concatenated input.
(2.57 times for smaller batches and up to 4.63 times for a larger batch.)

### Jaeger compare output

`wc -l` shows that the output files have the exact same number of lines as the
concatenated output files I generated earlier. The main difference is that the
original TSV files have more columns and are therefore bigger.
In short, same output: use concatenated input and save a lot of time!
