#! /usr/bin/env bash

mkdir -p data/ATB

# Download metadata files for AllTheBacteria
# (https://allthebacteria.readthedocs.io/en/latest/)

# 1. ENA metadata
wget -O data/ATB/ena_metadata.0.2.20240606.tsv.gz https://osf.io/download/j47ug/
wget -O data/ATB/ena_metadata.20240801.tsv.gz https://osf.io/download/ks7yt/

# 2. Sample lists
wget -O data/ATB/sample_list.txt.gz https://osf.io/download/ev3sj/

# 3. Sylph (species abundances per sample)
wget -O data/ATB/sylph.tsv.gz  https://osf.io/download/nu5a6/

# 4. Assembly statistics
wget -O data/ATB/assembly-stats.tsv.gz https://osf.io/download/nbyqv/

# 5. CheckM2 results
wget -O data/ATB/checkm2.tsv.gz https://osf.io/download/289f5/

# 6. Species calls
wget -O data/ATB/species_calls.tsv.gz https://osf.io/download/7t9qd/

# Get a list of all files:
wget -O data/ATB/file_list.all.20240805.tsv.gz https://osf.io/download/dw3h7
