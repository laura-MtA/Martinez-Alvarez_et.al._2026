# Uses PADLOC to identify defense systems in genomes

# Activate environment with PADLOC installed
conda activate padloc

# Assumes:
# - you are in the repo root
# - list.txt contains the names of the genome FASTA files (.fna), one per line
# - padloc_out/ exists (or is created before running)

while read -r genome; do padloc --fna "fna/${genome}" --outdir padloc_out --cpu 10; done < fna_genomes.list
