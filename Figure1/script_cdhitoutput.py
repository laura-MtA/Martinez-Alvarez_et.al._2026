#!/usr/bin/env python3

"""
script_cdhitoutput.py

Parse a CD-HIT .clstr file and write a long-format table:

  cluster_id <tab> protein_id

This script is used for Fig 1a (Venn overlap between PADLOC and DefenseFinder).
Protein IDs are expected to retain their original FASTA headers (e.g. PD_* or DF_* prefixes).

Usage:
  python script_cdhitoutput.py <input.clstr> <output.tsv>

Example:
  python script_cdhitoutput.py bac120_pd_df_c1.faa.clstr bac120_transformed.txt
"""

import sys


def transform_cdhit_output(cdhit_output, output_file):
    current_cluster = None

    with open(cdhit_output, "r") as f, open(output_file, "w") as out_f:
        for line in f:
            line = line.strip()

            # Cluster header line: ">Cluster 0"
            if line.startswith(">Cluster"):
                # Convert to "Cluster_0" (tab-safe, easier downstream)
                current_cluster = "Cluster_" + line.split()[1]
                continue

            
            if current_cluster is not None and ">" in line:
                protein_id = line.split(">")[1].split("...")[0]
                out_f.write(f"{current_cluster}\t{protein_id}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_cdhitoutput.py <input.clstr> <output.tsv>")
        sys.exit(1)

    cdhit_output = sys.argv[1]
    output_file = sys.argv[2]
    transform_cdhit_output(cdhit_output, output_file)

