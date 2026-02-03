# Script for Figure 4
# CRISPR-Cas prevalence across archaeal taxa using PADLOC output
#
#
# Required objects 
#   - ar53_fun_v4: see folder Figure2 of this repository
#   - padloc_ar53_noVSPRDMS: see folder Figure1 of this repository
#
# Packages
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Cas prevalence by CLASS
# ------------------------------------------------------------

# Count genomes per class
ar53abuncl <- table(ar53_fun_v4$class)
ar53abuncl <- as.data.frame(ar53abuncl)
names(ar53abuncl) <- c("class", "p.count")   # p.count = genomes per class

# Numerator: number of genomes with cas, aggregated at class + phylum
casclass <- aggregate(`cas` ~ class + phylum, data = ar53_fun_bin, sum)

# Merge denominators
casclass <- merge(casclass, ar53abuncl, by = "class", all.x = TRUE, sort = FALSE)

# Add "All Archaea" row 
all_row_class <- data.frame(
  class = "All Archaea",
  phylum = "All Archaea",
  cas = sum(casclass$cas, na.rm = TRUE),
  p.count = 7747   # total archaeal genomes in dataset
)
casclass <- rbind(casclass, all_row_class)

# Convert to percentage of genomes
casclass$rel.abundance <- (casclass$cas / casclass$p.count) * 100
casclass$full <- 100

# Ordering 
bar_order <- casclass$class
casclass$class <- factor(casclass$class, levels = rev(bar_order))

# Dashed reference line = "All Archaea" prevalence
line <- casclass[casclass$class == "All Archaea", "rel.abundance"]

# Plot
P30 <- ggplot(data = casclass[casclass$p.count > 9, ], aes(x = class)) +
  geom_bar(aes(y = full, fill = "Total Genomes"),
           stat = "identity", color = "black", width = 0.5, alpha = 1) +
  geom_bar(aes(y = rel.abundance, fill = "cas"),
           stat = "identity", color = "black", width = 0.5) +
  coord_flip() +
  geom_hline(yintercept = line, linetype = "dashed", color = "red") +
  scale_fill_manual(name = "Legend",
                    values = c("Total Genomes" = "lightgray", "cas" = "steelblue")) +
  labs(title = "cas",
       x = "Class",
       y = "Percentage of Genomes") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

ggsave("casclass.pdf", plot = P30, width = 9, height = 10)
P30


# ------------------------------------------------------------
# Cas prevalence by PHYLUM 
# ------------------------------------------------------------

# aggregate cas by phylum
Casphyl <- aggregate(`cas` ~ phylum, data = ar53_fun_bin, sum)

#phylum taxonomy
ar53abun <- table(ar53_fun_v4$phylum)
ar53abun <- as.data.frame(ar53abun)
names(ar53abuncl) <- c("phylum", "p.count")   # p.count = genomes per phylum

# Merge denominators
Casphyl <- merge(Casphyl, ar53abun, by = "phylum", all.x = TRUE, sort = FALSE)

# Add "All Archaea"
all_row_phyl <- data.frame(
  phylum = "All Archaea",
  cas = sum(Casphyl$cas, na.rm = TRUE),
  count = 7747
)
Casphyl <- rbind(Casphyl, all_row_phyl)

Casphyl$rel.abundance <- (Casphyl$cas / Casphyl$count) * 100
Casphyl$full <- 100

# Reference line
line <- Casphyl[Casphyl$phylum == "All Archaea", "rel.abundance"]

P3 <- ggplot(data = Casphyl[Casphyl$count > 9, ], aes(x = phylum)) +
  geom_bar(aes(y = full, fill = "Total Genomes"),
           stat = "identity", color = "black", width = 0.5, alpha = 1) +
  geom_bar(aes(y = rel.abundance, fill = "Cas"),
           stat = "identity", color = "black", width = 0.5) +
  coord_flip() +
  geom_hline(yintercept = line, linetype = "dashed", color = "red") +
  scale_fill_manual(name = "Legend",
                    values = c("Total Genomes" = "lightgray", "Cas" = "steelblue")) +
  labs(title = "Cas",
       x = "Phylum",
       y = "Percentage of Genomes") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

ggsave("casphylum.pdf", plot = P3, width = 9, height = 5)
P3


# ------------------------------------------------------------
# Cas type composition per PHYLUM
# ------------------------------------------------------------

# Subset Padloc table to cas systems
x <- filter(padloc_ar53_noVSPRDMS, system == "cas")

# Drop ambiguous subtypes and adaptation cassettes
x <- x[!x$system.subtype %in% c("cas_type_other", "cas_adaptation"), ]

# Extract cas "type" from PADLOC subtype strings
# Expected patterns like cas_type_I-A, cas_type_III-B, etc
x$cas.type <- sub("^cas_type_([A-Z]+)(-.*)?$", "type_\\1", x$system.subtype)

# NOTE: these column indices depend on your exact padloc table structure.
# If it breaks later, switch to selecting by name instead of index.
x <- x[, c(1, 2, 5, 8, 9, 14)]

# Denominators per phylum
x.totgen <- table(x$phylum)
x.totgen <- as.data.frame(x.totgen)
names(x.totgen) <- c("phylum", "total_genomes_in_phylum_cas_only")  # clarity label
x.totgen <- arrange(x.totgen, phylum)

# Count rows by phylum + cas.type
castype <- x %>%
  group_by(phylum, cas.type) %>%
  summarise(genome_count = n(), .groups = "drop")

# Merge denominators and compute proportions
castype <- merge(castype, x.totgen, by = "phylum", all.x = TRUE, sort = FALSE)

castype$percentage <- castype$genome_count / castype$total_genomes_in_phylum_cas_only
castype$percentage <- castype$percentage * 100

# Add "All Archaea" rows by cas.type
all_types <- c("type_I", "type_II", "type_III", "type_IV", "type_V")

all_archaea_rows <- lapply(all_types, function(tt) {
  data.frame(
    phylum = "All Archaea",
    cas.type = tt,
    genome_count = sum(castype$genome_count[castype$cas.type == tt], na.rm = TRUE),
    total_genomes_in_phylum_cas_only = sum(castype$genome_count, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
all_archaea_rows <- bind_rows(all_archaea_rows)
castype <- bind_rows(castype, all_archaea_rows)
castype$percentage <- castype$genome_count / castype$total_genomes_in_phylum_cas_only * 100

# Factor order for legend
castype$cas.type <- factor(castype$cas.type, levels = c("type_I", "type_II", "type_III", "type_IV", "type_V"))

P31 <- ggplot(castype, aes(x = phylum, y = percentage, fill = cas.type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.5) +
  labs(title = "Percentage of cas types per phylum",
       x = "Phylum",
       y = "Percentage",
       fill = "Cas Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

ggsave("castype.pdf", plot = P31, width = 9, height = 5)
P31


############################################################
#95% CIs per phylum
############################################################
library(dplyr)
library(binom)

# Keep only what we need
df_cas <- ar53_fun_bin %>%
  select(genome, phylum, cas)

# Remove phyla with too few genomes
phylum_counts <- df_cas %>%
  count(phylum, name = "n_genomes")

valid_phyla <- phylum_counts %>%
  filter(n_genomes >= 10) %>%
  pull(phylum)

df_cas <- df_cas %>%
  filter(phylum %in% valid_phyla)

cas_phylum <- df_cas %>%
  group_by(phylum) %>%
  summarise(
    n = n(),                   # number of genomes
    k = sum(cas, na.rm = TRUE) # genomes with cas
  ) %>%
  mutate(
    prevalence = k / n
  )

# Wilson confidence intervals
ci <- binom.confint(cas_phylum$k, cas_phylum$n, method = "wilson")

cas_phylum <- cas_phylum %>%
  mutate(
    lower = ci$lower,
    upper = ci$upper
  )

all_archaea <- df_cas %>%
  summarise(
    phylum = "All Archaea",
    n = n(),
    k = sum(cas, na.rm = TRUE),
    prevalence = k / n
  )

ci_all <- binom.confint(all_archaea$k, all_archaea$n, method = "wilson")

all_archaea <- all_archaea %>%
  mutate(
    lower = ci_all$lower,
    upper = ci_all$upper
  )

# Combine
cas_phylum_ci <- bind_rows(cas_phylum, all_archaea)

archaea_mean <- all_archaea$prevalence

cas_phylum_ci <- cas_phylum_ci %>%
  mutate(
    signif = case_when(
      lower > archaea_mean ~ "overrepresented",
      upper < archaea_mean ~ "underrepresented",
      TRUE                 ~ "not_significant"
    )
  )

cas_phylum_ci <- cas_phylum_ci %>%
  mutate(
    prevalence_pct = prevalence * 100,
    lower_pct = lower * 100,
    upper_pct = upper * 100
  )

library(ggplot2)
library(dplyr)

# Here: sort by prevalence and keep All Archaea at the top.
cas_phylum_ci_plot <- cas_phylum_ci %>%
  mutate(phylum = as.character(phylum)) %>%
  arrange(prevalence_pct) %>%
  mutate(phylum = factor(phylum, levels = c("All Archaea", setdiff(phylum, "All Archaea"))))

# Archaeal mean line (in percent)
archaea_mean_pct <- (cas_phylum_ci %>%
                       filter(phylum == "All Archaea") %>%
                       pull(prevalence_pct))[1]

# Plot
p_cas_phylum <- ggplot(cas_phylum_ci_plot, aes(x = prevalence_pct, y = phylum)) +
  geom_vline(xintercept = archaea_mean_pct, linetype = "dashed", color = "red") +
  
  geom_col(aes(fill = signif), width = 0.7, color = "black") +
  geom_errorbarh(aes(xmin = lower_pct, xmax = upper_pct), height = 0.2) +
  
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "overrepresented" = "steelblue4",
      "underrepresented" = "lightblue3",
      "not_significant" = "lightgray"
    ),
    name = NULL
  ) +
  labs(
    title = "cas prevalence across archaeal phyla",
    x = "Prevalence (% genomes with cas)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_cas_phylum

# Save
ggsave("Fig3_cas_phylum_CI.pdf", plot = p_cas_phylum, width = 8.5, height = 6)


############################################################
#95% CIs per class
############################################################
library(dplyr)
library(binom)

# Keep only what we need
df_cas_cl <- ar53_fun_bin %>%
  select(genome, class, cas)

# Remove phyla with too few genomes
class_counts <- df_cas_cl %>%
  count(class, name = "n_genomes")

valid_class <- class_counts %>%
  filter(n_genomes >= 10) %>%
  pull(class)

df_cas_cl <- df_cas_cl %>%
  filter(class %in% valid_class)

cas_class <- df_cas_cl %>%
  group_by(class) %>%
  summarise(
    n = n(),                   # number of genomes
    k = sum(cas, na.rm = TRUE) # genomes with cas
  ) %>%
  mutate(
    prevalence = k / n
  )

# Wilson confidence intervals
ci <- binom.confint(cas_class$k, cas_class$n, method = "wilson")

cas_class <- cas_class %>%
  mutate(
    lower = ci$lower,
    upper = ci$upper
  )

all_archaea <- df_cas_cl %>%
  summarise(
    class = "All Archaea",
    n = n(),
    k = sum(cas, na.rm = TRUE),
    prevalence = k / n
  )

ci_all <- binom.confint(all_archaea$k, all_archaea$n, method = "wilson")

all_archaea <- all_archaea %>%
  mutate(
    lower = ci_all$lower,
    upper = ci_all$upper
  )

# Combine
cas_class_ci <- bind_rows(cas_class, all_archaea)

archaea_mean <- all_archaea$prevalence

cas_class_ci <- cas_class_ci %>%
  mutate(
    signif = case_when(
      lower > archaea_mean ~ "overrepresented",
      upper < archaea_mean ~ "underrepresented",
      TRUE                 ~ "not_significant"
    )
  )

cas_class_ci <- cas_class_ci %>%
  mutate(
    prevalence_pct = prevalence * 100,
    lower_pct = lower * 100,
    upper_pct = upper * 100
  )

library(ggplot2)
library(dplyr)

# sort by prevalence
cas_class_ci_plot <- cas_class_ci %>%
  mutate(class = as.character(class)) %>%
  arrange(prevalence_pct) %>%
  mutate(class = factor(class, levels = c("All Archaea", setdiff(class, "All Archaea"))))

# Archaeal mean line (in percent)
archaea_mean_pct <- (cas_class_ci %>%
                       filter(class == "All Archaea") %>%
                       pull(prevalence_pct))[1]

# Plot
p_cas_class <- ggplot(cas_class_ci_plot, aes(x = prevalence_pct, y = class)) +
  geom_vline(xintercept = archaea_mean_pct, linetype = "dashed", color = "red") +
  
  geom_col(aes(fill = signif), width = 0.7, color = "black") +
  geom_errorbarh(aes(xmin = lower_pct, xmax = upper_pct), height = 0.2) +
  
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "overrepresented" = "steelblue4",
      "underrepresented" = "lightblue3",
      "not_significant" = "lightgray"
    ),
    name = NULL
  ) +
  labs(
    title = "cas prevalence across archaeal classes",
    x = "Prevalence (% genomes with cas)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_cas_class

# Save
ggsave("Fig3_cas_class_CI.pdf", plot = p_cas_class, width = 8.5, height = 6)
