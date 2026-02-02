# Uses CCTyper to identify CRISPR-Cas loci in genome FASTA files

# Activate environment with cctyper installed
conda activate cctyper

# Assumes:
# - you are in the repo root
# - genomes/ contains files named <ID>.fna (or .fasta)
# - list.txt contains one <ID> per line (matching the genome basenames)
# - cctyper_output/ exists (or create it before running)

cd fna/
cat ../list.txt | parallel --jobs 15 'cctyper "{}.fna" "../cctyper_output/{}.out"'
