#!/usr/bin/env bash

# script_run_cctyper.sh
#
# Identify CRISPR–Cas loci in genome FASTA files using CCTyper.
#
# Input:
# - genome FASTA files (e.g. <ID>.fna or <ID>.fasta)
#
# Assumptions:
# - you are in the repository root
# - genomes/ (or fna/) contains genome FASTA files named <ID>.fna
# - list.txt contains one genome ID per line (matching FASTA basenames)
# - cctyper_output/ exists (or create it before running)
#
# Notes:
# - CCTyper is run in "meta" mode for Prodigal, as used in the original analysis.
# - Thread usage is controlled via the --threads argument passed to CCTyper.

set -euo pipefail

# Activate conda environment with CCTyper installed
conda activate cctyper

# Run CCTyper on each genome listed in list.txt
cd fna/

cat ../list.txt | parallel --jobs 15 \
  'cctyper "{}.fna" "../cctyper_output/{}" --threads 25 --prodigal meta'

# END
