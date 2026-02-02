#!/usr/bin/env Rscript

# script_Rmerge.R
#
# Merge REBASE-based RM candidate annotations (after removing overlaps with other systems) onto the full protein-to-genome mapping table.
#
# Input:
# - allRM_ar53_nodef.txt            (no header; RM candidates; protein.id in column 4)
# - ../../all.proteins.genomes.tab  (protein-to-genome mapping; must include protein.id and genome.id)
#
# Output:
# - ar53.prot2RM.nodef.txt          (all proteins, with RM fields filled or NA; includes row.num)

library(readr)
library(dplyr)

# Read RM candidates (no header in file)
allRM <- read_tsv("allRM_ar53_nodef.txt", col_names = FALSE, col_types = NULL, guess_max = 100)

names(allRM) <- c(
  "module",     # M / R / S / IV / IIG / ...
  "type",       # IM / IR / IIR / IIG / IV / ...
  "ref",        # REBASE reference entry
  "protein.id", # target protein identifier
  "eval",
  "pident",
  "alen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "bitscore"
)

# Read full protein-to-genome mapping
prot2gen <- read_tsv("../../all.proteins.genomes.tab", col_names = TRUE, col_types = NULL, guess_max = 100)

# Merge, retaining all proteins (RM annotations are added where available)
output <- merge(prot2gen, allRM, by = "protein.id", all.x = TRUE, sort = FALSE)

# Arrange and add a row index
output <- output %>% arrange(genome.id, protein.id)

output$row.num <- seq_len(nrow(output))

# Move row.num to the first column
output <- output[, c("row.num", setdiff(names(output), "row.num"))]

# Write output
write.table(output,"ar53.prot2RM.nodef.txt",  quote = FALSE,  row.names = FALSE,  sep = "\t")
