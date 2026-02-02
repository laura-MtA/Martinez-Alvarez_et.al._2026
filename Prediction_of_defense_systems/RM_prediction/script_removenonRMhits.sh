#!/usr/bin/env bash

# script_remove_nonRM_defense_hits.sh
#
# Remove RM candidate proteins that overlap with non-RM defense/DNA-modification-related hits identified by DefenseFinder and PADLOC.
#
# This script produces:
# - indfpd.list            : protein IDs to exclude (DefenseFinder + PADLOC-nonRM overlap)
# - allRM_ar53_nodef.txt   : REBASE-based RM candidate set after excluding those proteins
#
# Required inputs (must be prepared beforehand):
# 1) allRM_ar53_dedup.txt
#    - output of script_dedupRM.sh (REBASE/MMseqs2 RM candidates)
#    - protein.id is column 4
#
# 2) nonRMdf_orfs.list
#    - list of protein IDs predicted by DefenseFinder (one protein ID per line), excluding RM entries
#    - how this list was generated is not recorded here; any list of protein IDs is acceptable
#
# 3) ar53_padloc.txt
#    - PADLOC output table for the same dataset (used to define a PADLOC "non-RM, non-DMS" set)
#
# 4) ar53_genomes.fna and ../../all.faa
#    - nucleotide genomes used to extract PADLOC non-RM ORFs (for translation)
#    - protein FASTA for the full dataset (used to extract sequences for BLAST query)
#
# Notes:
# - PADLOC "non-RM" here is defined as entries not labelled RM_ or DMS_.
# - This script documents the procedure used; it is not intended as a general pipeline.

set -euo pipefail


###############################################################################
# DefenseFinder ORFs to exclude (non-RM systems)
###############################################################################

# ar53_defensefinder_allDS.txt is the raw DefenseFinder output table.
# We generated df_orfs.txt by retaining non-RM systems only (i.e. excluding rows labelled as RM / RM_Type_*) and extracting:
#   - system name (column 4)
#   - comma-separated protein IDs (column 8)
#
# This produces a two-column file:
#   <system>\t<protein1,protein2,...>

awk -F'\t' '$4 != "RM" && $5 !~ /^RM_/ {print $4 "\t" $8}' ar53_defensefinder_allDS.txt | awk 'NF==2 && $2 != ""' > df_orfs.txt

# Convert df_orfs.txt to a one-protein-ID-per-line list (df_orfs.list)
cut -f2 df_orfs.txt | tr ',' '\n' | sed '/^$/d' | sort -u > df_orfs.list

###############################################################################
# A) Build PADLOC "non-RM, non-DMS" ORF protein database for BLAST (noRM_pd_orfs.faa2)
###############################################################################

# Extract PADLOC rows NOT corresponding to RM systems or DMS entries
grep -v "RM_" ar53_padloc.txt | grep -v "DMS_" | cut -f 2,12,13,14 | awk '!seen[$0]++' > temp

# Convert to BED format
# (second awk shifts start coordinate by -1 for bedtools)
awk '{print $1 "\t" $2 "\t" $3 "\t" "region_" NR "\t" "0" "\t" $4}' temp | awk -v OFS='\t' '{$2 = $2 - 1; print}' > noRM_pd_orfs.bed

# Extract nucleotide regions and translate to proteins
bedtools getfasta -fi ar53_genomes.fna -bed noRM_pd_orfs.bed -fo noRM_pd_orfs.fasta -s
sed -i 's/:/_/g' noRM_pd_orfs.fasta
transeq -sequence noRM_pd_orfs.fasta -outseq noRM_pd_orfs.faa
sed -i 's/*//g' noRM_pd_orfs.faa

# Remove duplicated FASTA headers (script_removeduplicates.py outputs noRM_pd_orfs.faa2)
python3 script_removeduplicates.py


###############################################################################
# B) Build BLAST query FASTA: RM candidate proteins (RMorfs.faa)
###############################################################################

# Extract unique RM candidate protein IDs from allRM_ar53_dedup.txt (protein.id is column 4)
cut -f 4 allRM_ar53_dedup.txt | sort -u > RMorfs.list

# Extract those protein sequences from the full FASTA
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' RMorfs.list ../../all.faa > RMorfs.faa


###############################################################################
# C) Identify overlaps with PADLOC non-RM proteins by BLAST
###############################################################################

makeblastdb -in noRM_pd_orfs.faa2 -dbtype prot -title noRM_pd_orfsDB -parse_seqids -out noRM_pd_orfsDB

# Identify RM candidates that match PADLOC non-RM proteins
blastp -db noRM_pd_orfsDB -query RMorfs.faa -out blast_RM2pd_evalE-10.txt -evalue 1E-10 -outfmt 6 -max_hsps 1 -num_threads 15

# Keep high-identity matches (>=98% identity); extract query protein IDs (col 1)
awk '$3>=98' blast_RM2pd_evalE-10.txt | cut -f 1 | sort -u > inpd.list


###############################################################################
# D) Combine DefenseFinder ORFs + PADLOC-nonRM overlaps into a single exclude list
###############################################################################

cat inpd.list nonRMdf_orfs.list | sort -u > indfpd.list

###############################################################################
# E) Remove excluded proteins from REBASE-based RM candidate table
###############################################################################

# Remove proteins in indfpd.list from allRM_ar53_dedup.txt to create allRM_ar53_nodef.txt
awk 'BEGIN{while((getline x<"indfpd.list")>0)rm[x]=1} !($4 in rm)' allRM_ar53_dedup.txt > allRM_ar53_nodef.txt

# END
