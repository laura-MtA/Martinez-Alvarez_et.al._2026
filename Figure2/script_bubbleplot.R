# Figure 2 — Bubble plot of defence system prevalence (PADLOC)

# Inputs provided with repository:
#   - ar53_padloc-sumtax.txt
#   - b120_padloc-sumtax.txt
#   - ar53_fun_v4.txt
#   - bac_fun_v3.txt

###############################################################
# Generation of count tables ar53_fun_v4.txt and bac_fun_v3.txt
# IMPORTANT! Note that this section of the script is indicative to show how the tables were created. 
# These tables are the result of merging several metadata to the output from padloc. 
# These metadata are not provided in this repository.
# Please import the tables to continue with the analysis from the next section
###############################################################

#Get archaeal data
padloc_ar53<-read_tsv("ar53_padloc-sumtax.txt", col_names=FALSE, col_types=NULL,guess_max=100)
names(padloc_ar53)<-c("genome.acc", "system", "system.number","seqid","system.subtype","genome","domain","phylum","class","order","family","genus","species")
padloc_ar53_noVSPR<-subset(padloc_ar53, padloc_ar53[,2] != "VSPR")
padloc_ar53_noVSPRDMS<-subset(padloc_ar53_noVSPR, padloc_ar53_noVSPR[,2] != "DMS_other")
padloc_ar53_noVSPRDMS<-subset(padloc_ar53_noVSPRDMS, padloc_ar53_noVSPRDMS[,2] != "DMS")
padloc_ar53_noVSPRDMS$genome<-substr(padloc_ar53_noVSPRDMS$genome.acc,1,13)
write.table(padloc_ar53_noVSPRDMS,"padloc_ar53_noVSPRDMS.txt", quote=FALSE, sep="\t", row.names=FALSE)
nrow(padloc_ar53_noVSPRDMS) 

#Make a count table (systems present in either domain; used as columns in the counts matrices)

arbac_systems<-unique(c(unique(padloc_ar53_noVSPRDMS$system), unique(padloc_bac_noVSPRDMS$system)))
length(arbac_systems) #~269 systems
arbac_systems

ar53_genomes<-read_tsv("ar53_genomes.list", col_names=FALSE, col_types=NULL,guess_max=100)
names(ar53_genomes)<-"genome_id"
nrow(ar53_genomes) #7747
head(ar53_genomes)
ar53_genomes$genome<-substr(ar53_genomes$genome_id,1,13)

#Make counts table for archaea
convert_to_counts_matrix <- function(input_table) {
  # Create an empty counts matrix
  counts_matrix <- matrix(0, nrow = length(unique(ar53_genomes$genome)), ncol = length(arbac_systems))
  rownames(counts_matrix) <- unique(ar53_genomes$genome)
  colnames(counts_matrix) <- arbac_systems

  for (i in 1:nrow(input_table)) {
    genome <- input_table$genome[i]
    system <- input_table$system[i]
    counts_matrix[genome, system] <- counts_matrix[genome, system] + 1
  }
  return(counts_matrix)
}
ar53_fun_v4<-convert_to_counts_matrix(padloc_ar53_noVSPRDMS)

#Convert rownames to column
ar53_fun_v4 <- as.data.frame(ar53_fun_v4)
ar53_fun_v4 <- ar53_fun_v4 %>% tibble::rownames_to_column(var = "genome")
View(ar53_fun_v4)

#Add rebase data to the archaeal dataset (can only be added to ar53_fun type of table)
#This substitutes Padloc's RM prediction (see Methods in the manuscript for more information)
RM_rebase<-read_tsv("RMcountPerGenome.txt",col_names=FALSE, col_types=NULL,guess_max=100)
names(RM_rebase)<-c("RM_rebase","genome.acc")
RM_rebase$genome<-substr(RM_rebase$genome.acc,1,13)
ar53_fun_v4<-merge(ar53_fun_v4,RM_rebase[,c(3,1)],by="genome", all.x=TRUE,sort=F)
ar53_fun_v4$RM_rebase[is.na(ar53_fun_v4$RM_rebase)] <- 0
ar53_fun_v4<-ar53_fun_v4[,-4] #Remove padloc RM data

#Add genome size data
ar_sizes<-read_tsv("ar_genome_sizes.txt", col_names=TRUE, col_types=NULL,guess_max=100)
ar53_fun_v4<-merge(ar53_fun_v4,ar_sizes,by="genome", all.x=TRUE, sort=F)

#Add completeness
ar53_completeness<-read_tsv("ar53_completeness.txt", col_names=TRUE, col_types=NULL, guess_max=100)
ar53_fun_v4<-merge(ar53_fun_v4,ar53_completeness[,c(1,3)], by="genome", all.x=TRUE, sort=FALSE)

#Add temperature
ogt<-read_tsv("ogt_out.txt", col_names=FALSE, col_types=NULL, guess_max=100)
ar53_fun_v4<-merge(ar53_fun_v4, ogt[,c(1,3)], by="genome", all.x=TRUE, sort=F)

#Calculate abundance per genome and system diversity per genome
#NOTE: columns 2:270 correspond to defence system columns in this table
ar53_fun_v4$all<-rowSums(ar53_fun_v4[,2:270])
ar53_fun_v4$syst.div<-rowSums(ar53_fun_v4[,2:270]>0, na.rm=TRUE)
ar53_fun_v4$all_noRMCas<-rowSums(ar53_fun_v4[,c(2:6,8:269)])

#Add taxonomy (phylum and class)
artaxonomy<-read_tsv("ar53_taxonomy_0523ED.tsv",col_names=TRUE,col_types=NULL,guess_max=100)
ar53_fun_v4<-merge(ar53_fun_v4, artaxonomy[,c(1:4)], by="genome", all.x=TRUE, sort=F)

write.table(ar53_fun_v4,"ar53_fun_v4.txt", quote=FALSE, sep="\t",row.names = FALSE)

#Bacterial data  
padloc_bac<-read_tsv("b120_padloc-sumtax.txt", col_names=FALSE, col_types=NULL,guess_max=100)
names(padloc_bac)<-c("genome.acc", "system", "system.number","seqid","system.subtype","genome","domain","phylum","class","order","family","genus","species")
padloc_bac_noVSPR<-subset(padloc_bac, padloc_bac[,2] != "VSPR")
padloc_bac_noVSPRDMS<-subset(padloc_bac_noVSPR, padloc_bac_noVSPR[,2] != "DMS_other")
padloc_bac_noVSPRDMS<-subset(padloc_bac_noVSPRDMS, padloc_bac_noVSPRDMS[,2] != "DMS")
padloc_bac_noVSPRDMS$genome<-substr(padloc_bac_noVSPRDMS$genome.acc,1,13)
write.table(padloc_bac_noVSPRDMS,"padloc_bac_noVSPRDMS.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Make counts table for bacteria
bactaxonomy<-read_tsv("padlocGenomeList.txt",col_names=FALSE,col_types=NULL,guess_max=100) 
bactaxonomy<-as.data.frame(bactaxonomy)
bactaxonomy<-arrange(bactaxonomy,genome.acc)

convert_to_counts_matrix <- function(input_table) {
  counts_matrix <- matrix(0, nrow = length(bactaxonomy$genome.acc), ncol = length(arbac_systems))
  rownames(counts_matrix) <- substr(bactaxonomy$genome.acc,1,13)
  colnames(counts_matrix) <- arbac_systems

  for (i in 1:nrow(input_table)) {
    genome <- input_table$genome[i]
    system <- input_table$system[i]
    counts_matrix[genome, system] <- counts_matrix[genome, system] + 1
  }
  return(counts_matrix)
}

bac_fun<-convert_to_counts_matrix(padloc_bac_noVSPRDMS)

#Convert rownames to column
bac_fun <- as.data.frame(bac_fun)
bac_fun <- bac_fun %>% tibble::rownames_to_column(var = "genome")

#Add genome size
bac_sizes<-read_tsv("bac_sizes.txt", col_names=TRUE, col_types=NULL,guess_max=100)
bac_fun<-merge(bac_fun, bac_sizes,by="genome",all.x=TRUE,sort=F)

#Add taxonomy
taxonomy<-read_tsv("bac120_tax.txt", col_names=TRUE, col_types=NULL,guess_max=100)
bac_fun<-merge(bac_fun, taxonomy[,c(9,2,3,4)], by="genome",all.x=TRUE,sort=F)

#Add completeness
bac_fun<-merge(bac_fun,b120_completeness, by="genome", all.x=TRUE, sort=F)

#Calculate abundance per genome and system diversity per genome
#NOTE: columns 2:270 correspond to defence system columns in this table
bac_fun$all<-rowSums(bac_fun[,2:270])
bac_fun$syst.div<-rowSums(bac_fun[,2:270]>0, na.rm=TRUE)
bac_fun$all_noRMCas<-rowSums(bac_fun[,c(2:3,5:7,9:270)])

write.table(bac_fun,"bac_fun_v3.txt", quote=FALSE, sep="\t",row.names = FALSE)

###############################################################################
# Bubble plot: relative abundance table
# Continue from this section after importing ar53_fun_v4.txt and bac_fun_v3.txt
###############################################################################
ar53_fun_v4<-read_tsv("ar53_fun_v4.txt", col_names=FALSE, col_types=NULL,guess_max=100)
bac_fun_v3<-read_tsv("bac_fun_v3.txt", col_names=FALSE, col_types=NULL,guess_max=100)

Arelab <- ar53_fun_v4 %>%
  mutate(across(2:270, ~ as.integer(.x > 0))) %>%
  pivot_longer(cols = 2:270, names_to = "system", values_to = "present") %>%
  group_by(domain, system) %>%
  summarise(tot.genomes = n(),
            no.genomes = sum(present, na.rm = TRUE),
            rel.abundance = no.genomes / tot.genomes,
            .groups = "drop")

Brelab <- bac_fun %>%
  mutate(across(2:270, ~ as.integer(.x > 0))) %>%
  pivot_longer(cols = 2:270, names_to = "system", values_to = "present") %>%
  group_by(domain, system) %>%
  summarise(tot.genomes = n(),
            no.genomes = sum(present, na.rm = TRUE),
            rel.abundance = no.genomes / tot.genomes,
            .groups = "drop")

ABrelab<-rbind(Arelab, Brelab)

#Obtain top20 systems per domain
topA<-arrange(filter(ABrelab, domain=="d__Archaea"),desc(rel.abundance))
topA<-topA$system[1:20]
topB<-arrange(filter(ABrelab, domain=="d__Bacteria"),desc(rel.abundance))
topB<-topB$system[1:20]
topAB <- unique(c(topA, topB))

ABrelab<-ABrelab[ABrelab$system %in% topAB, ]
ABrelab$rel.abundance<-ABrelab$rel.abundance*100

#Manually add "No systems" (raw)
ABrelab$system<-as.character(ABrelab$system)
ABrelab <- rbind(
  ABrelab,
  data.frame(domain="d__Archaea", system="No systems", tot.genomes=7747, no.genomes=990, rel.abundance=(990/7747)*100),
  data.frame(domain="d__Bacteria", system="No systems", tot.genomes=40000, no.genomes=650, rel.abundance=(650/40000)*100)
)

ABrelab<-arrange(ABrelab,domain,system)

bubble_order<-topAB
bubble_order[29]<-"No systems"
bubble_palette <- colorRampPalette(c("#8c510a", "#d8b365", "#f6e8c3"))(10)

ABrelab$system <- factor(ABrelab$system,  levels = rev(bubble_order))

Bubbleplot<-ggplot(ABrelab, aes(x = domain, y = system, size = rel.abundance, color = rel.abundance)) +
  geom_point() +
  scale_size_continuous(range = c(0.5, 12), guide = "legend", name = "Value") +
  scale_color_gradientn(colors = bubble_palette, guide = "legend", name = "Value") + #FIX bubble_colors->bubble_palette
  labs(title = "Most abundant defense systems per domain",
       x = "Domain",
       y = "Defense systems") +
  theme_minimal() +
  geom_text(aes(label = round(rel.abundance, 2)),
            position = position_nudge(x = 0.5, y = 0.5),
            size = 3,
            color = "black") +
  scale_y_discrete(limits = levels(ABrelab$system))
Bubbleplot

#HQ dataset
ArelabHQ <- filter(ar53_fun_v4,checkm_completeness >=90) %>%
  mutate(across(2:270, ~ as.integer(.x > 0))) %>%
  pivot_longer(cols = 2:270, names_to = "system", values_to = "present") %>%
  group_by(domain, system) %>%
  summarise(tot.genomes = n(),
            no.genomes = sum(present, na.rm = TRUE),
            rel.abundance = no.genomes / tot.genomes,
            .groups = "drop")

BrelabHQ <- filter(bac_fun, checkm_completeness>=90) %>%
  mutate(across(2:270, ~ as.integer(.x > 0))) %>%
  pivot_longer(cols = 2:270, names_to = "system", values_to = "present") %>%
  group_by(domain, system) %>%
  summarise(tot.genomes = n(),
            no.genomes = sum(present, na.rm = TRUE),
            rel.abundance = no.genomes / tot.genomes,
            .groups = "drop")

ABrelabHQ<-rbind(ArelabHQ, BrelabHQ)
ABrelabHQ<-ABrelabHQ[ABrelabHQ$system %in% topAB, ]
ABrelabHQ$rel.abundance<-ABrelabHQ$rel.abundance*100

#Change domain name to avoid overlap with the raw dataset
ABrelabHQ$domain[ABrelabHQ$domain == "d__Archaea"]  <- "HQ-Archaea"
ABrelabHQ$domain[ABrelabHQ$domain == "d__Bacteria"]  <- "HQ-Bacteria"

#Manually add "No systems" (HQ)
ABrelabHQ <- rbind(
  ABrelabHQ,
  data.frame(domain="HQ-Archaea", system="No systems", tot.genomes=3460, no.genomes=191, rel.abundance=(191/3460)*100),
  data.frame(domain="HQ-Bacteria", system="No systems", tot.genomes=34273, no.genomes=280, rel.abundance=(280/34273)*100)
)

ABrelabHQ<-arrange(ABrelabHQ,domain,system)
ABrelabHQ$system <- factor(ABrelabHQ$system,  levels = rev(bubble_order))

HQbubbleplot<-ggplot(ABrelabHQ, aes(x = domain, y = system, size = rel.abundance, color = rel.abundance)) +
  geom_point() +
  scale_size_continuous(range = c(0.5, 12), guide = "legend", name = "Value") +
  scale_color_gradientn(colors = bubble_palette, guide = "legend", name = "Value") +
  labs(title = "HQ dataset",
       x = "Domain",
       y = "Defense systems") +
  theme_minimal() +
  geom_text(aes(label = round(rel.abundance, 2)),
            position = position_nudge(x = 0.5, y = 0.5),
            size = 3,
            color = "black") +
  scale_y_discrete(limits = levels(ABrelabHQ$system))
HQbubbleplot

RawHQ_bubble<- (Bubbleplot | HQbubbleplot)
RawHQ_bubble
ggsave("RawHQ_BubblePlot.pdf", plot= RawHQ_bubble, dpi=300)
