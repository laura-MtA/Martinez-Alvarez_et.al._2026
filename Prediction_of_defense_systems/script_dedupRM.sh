# Log: post-processing REBASE/MMseqs2 hits to generate a unified RM annotation table
# Inputs: per-(sub)database deduplicated MMseqs2 hits (e.g. IM-ar53_dedup.m8)
# Output: allRM_ar53_dedup.txt and downstream mapping allRM2prot.txt

###############################################################################
# 1) Add two annotation columns to each deduplicated hit file
#    - component: M/R/S/etc.
#    - type-comp: IM/IR/IIR/IIG/IV/etc.
#
# Example: add "M" and "IM" as the first two columns
###############################################################################

awk '{print "M\tIM\t"$0}' filename > filenameED
# Repeat for each file, changing M/R/S and IM/IR/IIM/IIR/... as appropriate


###############################################################################
# 2) Combine Type I/II/III components and select the best hit per target protein
#
# Column order after adding the two cols:
# 1 component
# 2 type-comp
# 3 query
# 4 target
# 5 evalue
# 6 pident
# 7 alnlen
# 8 qstart
# 9 qend
# 10 tstart
# 11 tend
# 12 bits
###############################################################################

cat \
  TypeIIIM-ar53_dedupED.m8 \
  TypeIIIR-ar53_dedupED.m8 \
  TypeIIM-ar53_dedupED.m8 \
  TypeIIR-ar53_dedupED.m8 \
  TypeIM-ar53_dedupED.m8 \
  TypeIR-ar53_dedupED.m8 \
  TypeIS-ar53_dedupED.m8 \
  > temp

# Sort by target (col 4), then by bits (col 12) descending, then by type-comp (col 2)
# Deduplicate by target (col 4), then filter bits > 54
cat temp \
  | sort -k4,4 -k12,12nr -k2,2r \
  | awk '!seen[$4]++' \
  | awk '$12>54' \
  > temp2


###############################################################################
# 3) Process Type IV and Type IIG separately
#    (these were handled separately to avoid competing annotations)
###############################################################################

cat IV-ar53_dedupED.m8 TypeIIG-ar53_dedupED.m8 \
  | sort -k4,4 -k12,12nr \
  | awk '!seen[$4]++' \
  | awk '$12>54' \
  > temp3

# Remove any targets present in temp3 from temp2
cut -f4 temp3 > 2remove.list

awk 'BEGIN {
  while ((getline x < "2remove.list") > 0) rm[x]=1;
  close("2remove.list");
}
!($4 in rm)' temp2 > temp4

# temp4: RM genes without Type IIG/IV
# temp3: Type IIG/IV genes

cat temp4 temp3 | sort -k4,4 > allRM_ar53_dedup.txt


###############################################################################
# 4) Downstream: map RM-annotated proteins back to genome order using prot2gen
#
# allRM_ar53_nodef.txt: version after removing proteins labelled as part of
# other DNA-modification systems by PADLOC/DefenseFinder (see notRM.list).
#
# The merge is done in R to preserve the original protein order in prot2gen:
# Create allRM2prot.txt by merging allRM_ar53_nodef.txt with prot2gen.txt in R:
# Rscript script_Rmerge.R

###############################################################################
