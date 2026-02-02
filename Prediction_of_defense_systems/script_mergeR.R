# Merge REBASE-based RM annotations with protein-to-genome mapping
# Input:
#   - allRM_ar53_nodef.txt : RM annotations after filtering PADLOC/DefenseFinder hits
#   - all.proteins.genomes.tab : protein-to-genome mapping table
# Output:
#   - ar53.prot2RM.nodef.txt

library(readr)
library(dplyr)

###############################################################################
# 1) Read RM annotations (no header in file)
###############################################################################

allRM <- read_tsv(
  "allRM_ar53_nodef.txt",
  col_names = FALSE,
  col_types = NULL,
  guess_max = 100
)

names(allRM) <- c(
  "module",      # M / R / S
  "type",        # IM / IR / IIR / IIG / IV / ...
  "ref",         # reference/query sequence
  "protein.id",  # target protein ID
  "eval",
  "pident",
  "alen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "bitscore"
)

###############################################################################
# 2) Read protein-to-genome mapping
###############################################################################

prot2gen <- read_tsv(
  "../../all.proteins.genomes.tab",
  col_names = TRUE,
  col_types = NULL,
  guess_max = 100
)

###############################################################################
# 3) Merge, keeping all proteins (RM annotations added where available)
#
# all.x = TRUE ensures that proteins without RM annotations are retained
###############################################################################

output <- merge(
  prot2gen,
  allRM,
  by = "protein.id",
  all.x = TRUE,
  sort = FALSE
)

###############################################################################
# 4) Order output and add a row index
###############################################################################

output <- output %>%
  arrange(genome.id, protein.id)

output$row.num <- seq_len(nrow(output))

###############################################################################
# 5) Reorder columns (row number first, then original fields)
###############################################################################

output <- output[, c("row.num", setdiff(names(output), "row.num"))]

###############################################################################
# 6) Write output
###############################################################################

write.table(
  output,
  "ar53.prot2RM.nodef.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
