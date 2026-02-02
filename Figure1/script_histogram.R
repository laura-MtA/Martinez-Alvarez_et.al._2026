# Fig 1d — Histogram of number of genomes vs number of defence systems per genome
#
# Input tables (provided with the repository):
#   - padloc_bac_noVSPRDMS.txt
#   - padloc_ar53_noVSPRDMS.txt
#
# These tables contain PADLOC calls after removing categories not corresponding to defense systems (e.g. DMS and VSPR). Each row is treated as one defence system call for a genome.
#
# Output:
#   - fig1d.svg  (overlay histogram: Bacteria vs Archaea)

library(readr)


# Read input tables

padloc_bac_noVSPRDMS <- read_tsv("padloc_bac_noVSPRDMS.txt", col_types = cols())
padloc_ar53_noVSPRDMS <- read_tsv("padloc_ar53_noVSPRDMS.txt", col_types = cols())

# Count defence systems per genome (frequency of genome.acc)

baccount <- as.data.frame(table(padloc_bac_noVSPRDMS$genome.acc))
names(baccount) <- c("genome.acc", "syst.freq")

arcount <- as.data.frame(table(padloc_ar53_noVSPRDMS$genome.acc))
names(arcount) <- c("genome.acc", "syst.freq")

#Adjust brakes
max_count <- max(c(baccount$syst.freq, arcount$syst.freq))
breaks_shared <- seq(-0.5, max_count + 0.5, by = 1)  # bins centred on integers

# Compute histograms without plotting 
hist_bac <- hist(baccount$syst.freq, breaks = breaks_shared, plot = FALSE)
hist_ar <- hist(arcount$syst.freq, breaks = breaks_shared, plot = FALSE)

# 4) Plot (overlay)

svg("fig1d.svg", width = 7, height = 6)

plot(hist_bac,
     main = "Distribution of defence systems across genomes",
     xlab = "Number of defence systems per genome",
     ylab = "Number of genomes",
     col = "skyblue",
     border = "black",
     xlim = c(0, max_count))

plot(hist_ar,
     col = "pink",
     border = "black",
     add = TRUE)

axis(side = 1, at = seq(0, max_count, by = 10))

legend("topright",
       legend = c("Archaea", "Bacteria"),
       fill = c("pink", "skyblue"),
       border = "black")

dev.off()
