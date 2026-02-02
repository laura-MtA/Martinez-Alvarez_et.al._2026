# Exemplary MMseqs2 workflow for identifying RM components using REBASE-derived protein sets
#
# Note:
# - REBASE protein entries were split into multiple sub-databases (by function/type)
#   and the steps below were repeated for each sub-database.
# - The list of input FASTA files from REBASE and the resulting MMseqs2 DBs used are noted below.
# - Paths and filenames must be adapted to your local setup.

# Databases created (REBASE-derived; names follow REBASE conventions)
# - V-DB        (Nickases)
# - C5M-DB      (cytosine-5 methyltransferases)
# - TypeIIM-DB  (Type II methyltransferases)
# - aMa-DB      (amino-MTases, subtype alpha)
# - aMb-DB      (amino-MTases, subtype beta)
# - aMg-DB      (amino-MTases, subtype gamma)
# - TypeIIR-DB  (Type II restriction enzymes)
# - TypeIV-DB   (Type IV)
# - TypeIM-DB   (Type I M subunit)
# - TypeIR-DB   (Type I R subunit)
# - TypeIS-DB   (Type I S subunit)
# - N-DB        (N genes)
# - TypeIIG-DB  (Type IIG enzymes)
# - TypeIIIM-DB (Type III M subunit)
# - TypeIIIR-DB (Type III R subunit)

# Activate environment with mmseqs2 installed
# conda activate mmseqs2

################################################################################
# A) Example: create ONE REBASE-derived MMseqs2 database (repeat per DB)
#
# Replace:
# - INPUT_FASTA with the corresponding REBASE protein FASTA/text file
# - DB_NAME with the desired MMseqs2 database name (e.g. aMa-DB, TypeIIR-DB)
################################################################################

INPUT_FASTA="All_Type_III_M_subunit_genes_Protein.txt"
DB_NAME="IIIM-ar53DB"
TMP_DIR="${DB_NAME}-tmp"
CLU_DB="${DB_NAME}-clu"

mkdir -p "${DB_NAME}" "${TMP_DIR}"

mmseqs createdb "${INPUT_FASTA}" "${DB_NAME}/${DB_NAME}"

mmseqs cluster \
  "${DB_NAME}/${DB_NAME}" \
  "${DB_NAME}/${CLU_DB}" \
  "${TMP_DIR}" \
  --min-seq-id 0.65 --cov-mode 0 -c 0.8

mmseqs createtsv \
  "${DB_NAME}/${DB_NAME}" \
  "${DB_NAME}/${CLU_DB}" \
  "${DB_NAME}/${CLU_DB}.tsv"


################################################################################
# B) Example: search archaeal proteins against ONE DB (repeat per DB)
#
# Replace:
# - REBASE_DB with the path to the chosen REBASE-derived MMseqs2 database
# - QUERY_DB with your archaeal protein DB (created separately with mmseqs createdb)
# - OUT_PREFIX with a meaningful output label per DB
################################################################################

REBASE_DB="IIIM-ar53DB/IIIM-ar53DB"
QUERY_DB="/path/to/ar53queryDB"          # e.g. mmseqs createdb archaeal_proteins.faa ar53queryDB
OUT_PREFIX="IIIM-ar53DB"
TMP_SEARCH="tmp"

mmseqs search "${REBASE_DB}" "${QUERY_DB}" "${OUT_PREFIX}DB" "${TMP_SEARCH}"

mmseqs convertalis \
  "${REBASE_DB}" \
  "${QUERY_DB}" \
  "${OUT_PREFIX}DB" \
  "${OUT_PREFIX}_aln.m8" \
  --format-output query,target,evalue,pident,alnlen,qstart,qend,tstart,tend,bits

# Best hit per target by highest bitscore
 sort -k2,2 -k10,10nr "${OUT_PREFIX}_aln.m8" > "${OUT_PREFIX}_sorted.m8"

# Remove duplicates by target, keeping the first (best-ranked) hit
awk '!seen[$2]++' "${OUT_PREFIX}_sorted.m8" > "${OUT_PREFIX}_dedup.m8"
