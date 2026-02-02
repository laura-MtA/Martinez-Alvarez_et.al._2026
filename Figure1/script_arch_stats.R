# Figure 1E–F — Archaea only
#
# Panel 1E:
#   Average number of defence systems per genome across archaeal phyla (mean ± SD).
#
# Panel 1F:
#   Average number of defence system families per genome across archaeal phyla (mean ± SD).
#
# Input:
#   padloc_ar53_noVSPRDMS
#     PADLOC output table for Archaea after filtering out categories not belonging to defense systems (e.g. VSPR and DMS; see repository Methods/notes).
#
# Assumptions about padloc_ar53_noVSPRDMS:
#   - Column "genome.acc" identifies genomes
#   - Column "phylum" provides archaeal phylum assignment per genome
#   - Column "system" is used here as a defence system family label (unique system types per genome)

library(readr)
library(dplyr)
library(ggplot2)

###############################################################################
# Panel 1E — Average number of defense systems per genome per phylum
###############################################################################

# Count number of PADLOC system calls per genome
arcount <- table(padloc_ar53_noVSPRDMS$genome.acc)
arcount <- as.data.frame(arcount)
names(arcount) <- c("genome.acc", "syst.freq")  # syst.freq = systems per genome

# Extract genome → phylum mapping (one row per genome)
tax.temp <- unique(padloc_ar53_noVSPRDMS[, c("genome.acc", "phylum")])

# Add phylum column to per-genome counts
arcount <- merge(arcount, tax.temp, by = "genome.acc", all.x = TRUE, sort = FALSE)

# Phylum abundance (number of genomes per phylum), used to order bars
ar53abun <- table(arcount$phylum)
ar53abun <- data.frame(phylum = names(ar53abun), count = as.numeric(ar53abun))
ar53abun <- arrange(ar53abun, desc(count))

# Mean and SD of systems per genome per phylum
ar53res_mean <- aggregate(syst.freq ~ phylum, data = arcount, FUN = mean)
ar53res_sd   <- aggregate(syst.freq ~ phylum, data = arcount, FUN = sd)

aggA_res2 <- merge(ar53res_mean, ar53res_sd, by = "phylum", all.x = TRUE, sort = FALSE)
colnames(aggA_res2) <- c("phylum", "mean", "sd")

# Add "All Archaea" summary row (overall mean/sd across genomes)
# (Row index is not fixed; using rbind keeps this robust if phyla change.)
aggA_res2 <- rbind(aggA_res2,data.frame(phylum = "All Archaea",mean = mean(arcount$syst.freq),sd   = sd(arcount$syst.freq)))

# Order phyla by abundance (and append "All Archaea" at the end)
aggA_res2$phylum <- as.character(aggA_res2$phylum)
aggA_res2$phylum <- factor(aggA_res2$phylum, levels = c(ar53abun$phylum, "All Archaea"))

# Plot: mean systems/genome ± SD per phylum
AavgDSnum <- ggplot(aggA_res2, aes(x = phylum, y = mean)) +
  geom_bar(stat = "identity", fill = "lightgray", color = "black") +
  geom_point(color = "red", size = 1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "black") +
  labs(title = "Average number of defence systems per genome per phylum",
       x = "Phylum",
       y = "Average number of defence systems per genome",
       caption = "Error bars show standard deviation within each phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Overall archaeal mean (dashed line)
  geom_hline(yintercept = mean(arcount$syst.freq), linetype = "dashed", color = "blue") +
  scale_y_continuous(breaks = seq(0, 80, by = 5))

AavgDSnum
ggsave("ar53_avgDSnum.pdf", plot = AavgDSnum, width = 15, height = 5)

###############################################################################
# Panel 1F — Average number of defense system families per genome per phylum
###############################################################################

# Define "family diversity" as the number of unique PADLOC "system" labels per genome
arcdiv <- distinct(padloc_ar53_noVSPRDMS, genome.acc, system)

# Count unique system families per genome
argendiv <- table(arcdiv$genome.acc)
argendiv <- as.data.frame(argendiv)
names(argendiv) <- c("genome.acc", "DSfams.freq")  # DSfams.freq = unique families per genome

# Add phylum to per-genome family counts
argen5 <- merge(argendiv, tax.temp, by = "genome.acc", all.x = TRUE, sort = FALSE)

# Mean and SD of family diversity per genome per phylum
ar53res4_mean <- aggregate(DSfams.freq ~ phylum, data = argen5, FUN = mean)
ar53res4_sd   <- aggregate(DSfams.freq ~ phylum, data = argen5, FUN = sd)

ar53res4 <- merge(ar53res4_mean, ar53res4_sd, by = "phylum", all.x = TRUE, sort = FALSE)
names(ar53res4) <- c("phylum", "DSfams.freq", "sd")

# Add "All Archaea" summary row
ar53res4 <- rbind(ar53res4, data.frame(phylum = "All Archaea", DSfams.freq = mean(argendiv$DSfams.freq),sd = sd(argendiv$DSfams.freq)))

# Order phyla by abundance (and append "All Archaea" at the end)
ar53res4$phylum <- as.character(ar53res4$phylum)
ar53res4$phylum <- factor(ar53res4$phylum, levels = c(ar53abun$phylum, "All Archaea"))

# Plot: mean families/genome ± SD per phylum
AdivDS <- ggplot(ar53res4, aes(x = phylum, y = DSfams.freq)) +
  geom_bar(stat = "identity", fill = "lightgray", color = "black") +
  geom_point(color = "red", size = 1) +
  geom_errorbar(aes(ymin = DSfams.freq - sd, ymax = DSfams.freq + sd), width = 0.2, color = "black") +
  labs(title = "Average number of defence system families per genome per phylum",
       x = "Phylum",
       y = "Average number of defence system families per genome",
       caption = "Error bars show standard deviation within each phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Overall archaeal mean family diversity (dashed line)
  geom_hline(yintercept = mean(argendiv$DSfams.freq), linetype = "dashed", color = "blue")

AdivDS
ggsave("ar53_divDSnum.pdf", plot = AdivDS, width = 15, height = 5)

# END
