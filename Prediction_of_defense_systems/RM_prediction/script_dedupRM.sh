#!/usr/bin/env bash

# script_dedupRM.sh
#
# Combine REBASE/MMseqs2 hit tables across RM component databases and assign a single RM annotation per protein.
#
# Priority rule:
# - Type IIG (IIG) and Type IV (IV) annotations override all other RM-component annotations. If a protein is annotated as IIG or IV, that label is retained and any other RM annotation for the same protein is discarded.
#
# Input:
# - Per-database MMseqs2 hit tables (*_dedup.txt), each already reduced to one best hit per target within that database.
#
# Output:
# - allRM_ar53_dedup.txt
#
# Column order in input and output tables:
#  1  module        (M / R / S / IV / IIG / ...)
#  2  type          (IM / IR / IIR / IIG / IV / ...)
#  3  ref           (REBASE reference entry)
#  4  protein.id    (target protein identifier)
#  5  evalue
#  6  pident
#  7  alen
#  8  qstart
#  9  qend
# 10  sstart
# 11  send
# 12  bitscore


###############################################################################
# 1) Merge all RM-component hit tables EXCEPT IIG and IV
###############################################################################

cat \
  C5M-ar53_dedup.txt \
  H-ar53_dedup.txt \
  IIIM-ar53_dedup.txt \
  IIIR-ar53_dedup.txt \
  IIM-ar53_dedup.txt \
  IIR-ar53_dedup.txt \
  IM-ar53_dedup.txt \
  IR-ar53_dedup.txt \
  IS-ar53_dedup.txt \
  N-ar53_dedup.txt \
  > temp_nonIIGIV.txt

# Deduplicate by target protein (col 4), keeping the best hit by bitscore (col 12), and applying a minimum bitscore threshold (>54).
sort -k4,4 -k12,12nr temp_nonIIGIV.txt  | awk '!seen[$4]++'  | awk '$12>54'  > temp_nonIIGIV_dedup.txt


###############################################################################
# 2) Process IIG and IV hit tables separately (overriding annotations)
###############################################################################

cat  IIG-ar53_dedup2.txt IV-ar53_dedup2.txt > temp_IIGIV.txt

# Deduplicate by target protein (col 4), keeping the best hit by bitscore and applying the same bitscore threshold.
sort -k4,4 -k12,12nr temp_IIGIV.txt | awk '!seen[$4]++' | awk '$12>54' > temp_IIGIV_dedup.txt


###############################################################################
# 3) Remove IIG/IV targets from the non-IIG/IV set, then merge
###############################################################################

# Extract protein IDs annotated as IIG/IV...
cut -f 4 temp_IIGIV_dedup.txt > IIGIV_targets.list

# and remove these from the non-IIG/IV pool so IIG/IV annotations take priority.
awk 'BEGIN { while ((getline x < "IIGIV_targets.list") > 0) rm[x] = 1; close("IIGIV_targets.list"); }
     !($4 in rm)' \
  temp_nonIIGIV_dedup.txt \
  > temp_nonIIGIV_filtered.txt

# Final merged RM annotation table (one annotation per protein)
cat temp_nonIIGIV_filtered.txt temp_IIGIV_dedup.txt | sort -k4,4 > allRM_ar53_dedup.txt

# END
