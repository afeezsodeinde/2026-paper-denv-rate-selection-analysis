#!/bin/bash

# DENV-1
augur filter \
  --metadata DENV1_70plus_genomes.csv \
  --min-date 1944 \
  --group-by date major_lineage country \
  --sequences-per-group 3 \
  --output-metadata DENV1_sampled.tsv \
  --output-strains DENV1_sampled.txt

# DENV-2
augur filter \
  --metadata DENV2_70plus_genomes.csv \
  --min-date 1944 \
  --group-by date major_lineage country \
  --sequences-per-group 3 \
  --output-metadata DENV2_sampled.tsv \
  --output-strains DENV2_sampled.txt

# DENV-3
augur filter \
  --metadata DENV3_70plus_genomes.csv \
  --min-date 1953 \
  --group-by date major_lineage country \
  --sequences-per-group 3 \
  --output-metadata DENV3_sampled.tsv \
  --output-strains DENV3_sampled.txt

# DENV-4
augur filter \
  --metadata DENV4_70plus_genomes.csv \
  --min-date 1956 \
  --group-by date major_lineage country \
  --sequences-per-group 3 \
  --output-metadata DENV4_sampled.tsv \
  --output-strains DENV4_sampled.txt
