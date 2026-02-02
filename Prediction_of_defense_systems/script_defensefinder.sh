#Uses DefenseFinder to identify defense systems in proteomes

#Activate environment with DefenseFinder installed
conda activate defensefinder

#Assumes:
# - you are in the repo root
# - faa/ contains files named <ID>.faa
# - infaa.list contains one <ID> per line (matching the .faa basenames)
# - df/ exists (or create it before running)

cd faa
cat ../infaa.list | parallel --jobs 40 'defense-finder run -o "../df/{}.out" "{}.faa"'
