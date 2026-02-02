# Fig 1a — overlap between PADLOC and DefenseFinder predictions (Venn diagram)
#
#   Proteins predicted by either tool were combined, clustered at 100% identity using CD-HIT and clusters were classified as:
#     - DefenseFinder-only
#     - PADLOC-only
#     - overlap (clusters containing proteins from both tools)
#
# Key inputs:
#   1) Combined FASTA of predicted proteins, with headers prefixed by source:
#        >PD_<protein_id> ...  (PADLOC)
#        >DF_<protein_id> ...  (DefenseFinder)
#
#      Example (PADLOC):
#        >PD_ABFD02000031.1_3057-4194(+)_1
#      Example (DefenseFinder):
#        >DF_CM000441.2_1514
#
#      Files used here:
#        - bac120_pd_df.faa(.gz)     (Bacteria)
#        - ar53_df_pd_all.faa        (Archaea)
#
#      Note: exact commands used to build these combined FASTA files are not retained. Conceptually, they were created by concatenating the predicted protein sets from each tool and prefixing each header with PD_ or DF_ to track provenance.
#
#   2) PADLOC system table filtered for comparison across tools:
#        - *_padloc-noDSM-VSPR-PDC_ED.txt
#      These tables exclude PADLOC entries labelled as:
#        - DMS (DNA modification systems; not complete defence systems)
#        - VSPR (not defence systems)
#        - PDC (not present in DefenseFinder, removed for tool-to-tool comparison)
#      The exact filtering command producing this file is not shown; the intent was to remove categories that are not comparable between PADLOC and DefenseFinder.
#
# Key outputs:
#   - CD-HIT cluster file: *.clstr
#   - transformed_output.txt (cluster_id + protein_id; derived from .clstr)
#   - Summary counts used for Venn (DF-only, PD-only, overlap)
#
# ------------------------------------------------------------------------------

###############################################################################
# 1) Cluster combined protein FASTA at 100% identity (CD-HIT)
###############################################################################

# Bacteria (bac120)
cd-hit -i bac120_pd_df.faa -o bac120_pd_df_c1.faa -c 1 -M 930

# Archaea (ar53)
cd-hit -i ar53_df_pd_all.faa -o ar53_df_pd_all_c1.faa -c 1


###############################################################################
# 2) Parse CD-HIT .clstr output to a long table (cluster_id <tab> protein_id)
###############################################################################

# script_cdhitoutput.py reads bac120_pd_df_c1.faa.clstr (or ar53_df_pd_all_c1.faa.clstr)
# and writes a two-column table:
#   >Cluster 0  <protein_id>
#   >Cluster 0  <protein_id>
#   >Cluster 1  <protein_id>
#   ...
#
# The exact version used in the analysis produced "transformed_output.txt" in this form.

python script_cdhitoutput.py

###############################################################################
# 3) Classify clusters by tool membership (DF-only, PD-only, both)
###############################################################################

# Clusters that contain at least one DefenseFinder protein (based on DF_ prefix)
grep "DF_" transformed_output.txt | cut -f 1 | sort | uniq > clustersWDF.temp

# Clusters that contain at least one PADLOC protein (based on PD_ prefix)
grep "PD_" transformed_output.txt | cut -f 1 | sort | uniq > clustersWPD.temp

# Clusters containing both DF and PD proteins
comm -12 <(sort clustersWDF.temp) <(sort clustersWPD.temp) > clustersWboth.temp

# DF-only clusters
comm -23 <(sort clustersWDF.temp) <(sort clustersWboth.temp) > clustersDF_only.temp

# PD-only clusters
comm -23 <(sort clustersWPD.temp) <(sort clustersWboth.temp) > clustersPD_only.temp

# Report cluster counts used in Fig 1a Venn
wc -l clustersWDF.temp
wc -l clustersWPD.temp
wc -l clustersWboth.temp
wc -l clustersDF_only.temp
wc -l clustersPD_only.temp


#END
