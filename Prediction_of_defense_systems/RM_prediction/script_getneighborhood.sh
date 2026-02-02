#!/usr/bin/env bash

# script_getneighborhood.sh
#
# Extract genomic neighbourhoods (±5 proteins) around RM candidate proteins.
#
# Input:
# - allRM_ar53_nodef.txt        (RM candidates; protein.id is column 4; no header)
# - ar53.prot2RM.nodef.txt      (protein-to-genome table with RM fields; header present)
#
# Output:
# - neighborhood.txt            (subset of ar53.prot2RM.nodef.txt containing RM candidates
#                               and up to 5 rows of context above and below each match)
#
# Notes:
# - This script uses grep context (-C 5) to capture ±5 rows around each RM candidate protein ID.
# - The output is used as input for PredictRM.v3.py.

set -euo pipefail

# Extract RM candidate protein IDs (column 4 of allRM_ar53_nodef.txt)
cut -f 4 allRM_ar53_nodef.txt > RM_candidates.list

# Extract neighbourhoods (±5 rows) around each candidate protein ID
grep -C 5 -f RM_candidates.list ar53.prot2RM.nodef.txt > neighborhood.txt

# END
