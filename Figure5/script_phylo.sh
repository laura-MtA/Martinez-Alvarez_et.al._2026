# ------------------------------------------------------------
# Indicative commands to build a phylogenetic tree for ONE system
# (replace filenames and chosen parameters as needed)
# ------------------------------------------------------------
# Steps 1-3 exemplify how the final set of sequences was selected.
# The final set of sequences for each tree is deposited in Supplementary Table 10 and can be used for steps 5-6 
# 1) Reduce redundancy (65% identity)
cd-hit -i system_raw.faa -o system_c0.65.faa -c 0.65

# 2) Preliminary tree
mafft --auto system_c0.65.faa > system_c0.65.aln
trimal -in system_c0.65.aln -out system_c0.65_gappyout.aln -gappyout
FastTree -lg -boot 100 system_c0.65_gappyout.aln > system_prelim.nwk

# 3) Prune preliminary tree to retain diversity and facilitate visualization
singularity exec ~/software/Treemmer/treemmer_sb/ python3 /Treemmer_v0.3.py system_prelim.nwk -X 250

# 4) Retrieve sequences retained by Treemmer (list → FASTA)
# (Treemmer produces a list file; replace the filename below accordingly)
# Example list name: system_prelim_X250.list
grep -v '^$' system_prelim_X250.list | sed 's/^>//g' > system_kept.ids
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' system_kept.ids system_c0.65.faa > system_pruned.faa

# These final trees can be reproduced by using the amino acid sequences deposited in Supplementary Table 10
# 5) Final alignment + trimming (high accuracy)
mafft --maxiterate 1000 --localpair --thread 10 system_pruned.faa > system_final.aln
trimal -in system_final.aln -gt 0.3 -out system_trim.aln

# 6) Final tree (IQ-TREE2)
iqtree2 -s system_trim.aln -m MFP -bb 1000 -alrt 1000 -bnni -T AUTO
