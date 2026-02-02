# Script for Figure 3
# Taxonomic distribution of the archaeal core-defensome
#
# This script is shown as a working example for one system (RM).
# To reproduce Fig 3, repeat the "RM block" below for each core system,
# changing:
#   i)  which system column you select (df_SYS <- ...)
#   ii) the model formula (SYS ~ phylum)
#   iii) the overall_mean calculation
#   iv) output filenames and plot title

#Packages
library(brglm2)
library(emmeans)
library(dplyr)
library(ggplot2)
library(readr)

# Input table
# -----------------------------
# ar53_fun_v4 must already be in the environment (loaded earlier), or load it here (the table is deposited in the folder "Figure2" of this repository):
# ar53_fun_v4 <- read_tsv("ar53_fun_v4.txt", col_names = TRUE, col_types = NULL, guess_max = 100)


# Create presence/absence matrix from the counts table (counts to presence/absence)
ar53_fun_bin <- ar53_fun_v4
ar53_fun_bin[, 2:270] <- (ar53_fun_bin[, 2:270] > 0) * 1

# Define low-count phyla and phylum order
# Phyla with <10 genomes are excluded from the analysis (as per Methods).
phylum_counts <- table(ar53_fun_bin$phylum)
low_counts <- names(phylum_counts[phylum_counts < 10])

# phylum order for plotting (most abundant to least abundant)
levels.conf <- names(sort(phylum_counts[!(names(phylum_counts) %in% low_counts)], decreasing = TRUE))

# Example: RM system (working block)
# Here you select only genome id, phylum and the system column.

df_RM <- ar53_fun_bin[, c(1, 279, 270)] #Col 1 is the genome, Col 2 is the phylum and Col 3 is the system
df_RM <- as.data.frame(df_RM)

# Keep only phyla with >=10 genomes
df_RM <- df_RM[!(df_RM$phylum %in% low_counts), ]
df_RM$phylum <- factor(df_RM$phylum, levels = rev(levels.conf))

# Model
model <- glm(`RM` ~ phylum, data = df_RM, family = "binomial", method = "brglmFit")
# Estimated marginal means per phylum
em_resp <- emmeans(model, specs = "phylum", type = "response")
em_df <- as.data.frame(em_resp)
# Overall prevalence across all archaeal genomes
overall_mean <- mean(df_RM$`RM`, na.rm = TRUE)
# Significance flag
em_df <- em_df %>%
  mutate(diff_from_overall = prob - overall_mean)

em_df$significant <- em_df$asymp.LCL > overall_mean | em_df$asymp.UCL < overall_mean

# Save table used for plotting
write.table(em_df, "RM_vsAll_CI.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot (one system)
CIplotRM <- ggplot(em_df, aes(x = prob, y = phylum)) +
  geom_vline(xintercept = overall_mean, linetype = "dashed", color = "gray40") +
  geom_col(aes(fill = significant), width = 0.6, color = NA) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL), height = 0.2) +
  geom_point(aes(color = significant), size = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = c("TRUE" = "steelblue4", "FALSE" = "lightgray")) +
  scale_color_manual(values = c("TRUE" = "chocolate2", "FALSE" = "black")) +
  labs(
    x = "Estimated prevalence (probability)",
    y = NULL,
    title = "Prevalence of RM system across phyla",
    subtitle = paste0("Dashed line = overall prevalence: ", round(overall_mean, 2)),
    fill = "Significant",
    color = "Significant"
  ) +
  theme_minimal()

CIplotRM

# Combine plots (optional)
# The combined plot section assumes you created CIplotCas, CIplotSoFic, etc.

# library(patchwork)
# CIplotsCoreDef <- (
#   (CIplotRM | CIplotCas | CIplotSoFic) /
#   (CIplotAbiE | CIplotVip | CIplotPT) /
#   (CIplotCBASS | CIplotAgo)
# )
# print(CIplotsCoreDef)
# ggsave("CIplotsCoreDef2.pdf", plot = CIplotsCoreDef, width = 28, height = 15)
