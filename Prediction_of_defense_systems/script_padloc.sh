# Uses PADLOC to identify defense systems in genomes

# Activate environment with PADLOC installed
conda activate padloc

# Assumes:
# - you are in the repo root
# - list.txt contains the names of the genome FASTA files (.fna), one per line
# - padloc_output/ exists (or is created before running)

parallel -j 15 padloc --fna {} --outdir ../padloc_output ::: $(cat list.txt)
