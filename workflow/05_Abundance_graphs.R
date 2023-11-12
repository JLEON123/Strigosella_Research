# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Looking at abundance of the bacterial communities
# Author: Josh Leon
# Software versions:  R v 4.1.1
#                     tidyverse v 1.3.2
#                     phyloseq v 1.36.0
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")

source("./scripts/theme2.R")

# Load cleaned phyloseq object ####
ps <- readRDS("./Final_phyloseq_object.RDS")

# Getting to know your phyloseq data ####

# number of taxa
ntaxa(ps)

# number of samples
nsamples(ps)

# sample names
sample_names(ps)
rank_names(ps)

# Find abundance of the community ###
any(taxa_sums(ps) == 0) # Checks if there are any unobserved OTUs
sum(taxa_sums(ps) == 0) # 0
any(sample_sums(ps) == 0) # No empty samples

# Look at reads by OTUs and Samples
readsumsdf = data.frame(nreads = sort(taxa_sums(ps), TRUE), sorted = 1:ntaxa(ps), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps), 
                                                        TRUE), sorted = 1:nsamples(ps), type = "Samples"))

title0 = "Total number of reads"
numreads <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") + 
  ggtitle(title0) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

ggsave(numreads, filename = "./graphs/number_of_reads.png", dpi = 300, width = 12, unit = "in")


# Phyla and genus abundances ###
set.seed(123)
# Shows structure on the y-axis and percentage of sequences on the x-axis
ps_RM = merge_samples(ps, "Structure")
ps_RM@sam_data$Structure <- row.names(ps_RM@sam_data)
ps_RM = transform_sample_counts(ps_RM, function(x) 100 * x/sum(x))

tax_composition <- 
  plot_bar(ps_RM, "Structure", fill = "Phylum") + 
  coord_flip() + 
  geom_bar(aes(color = Phylum, fill= Phylum), stat = "identity", position = "stack") +
  scale_color_manual(values = pal.abundance, aesthetics = c("color", "fill")) +
  scale_x_discrete(expand = c(0,0)) +
  ylim(0,100) +
  labs(x = "",
       y = "Relative Abundance",
       fill = "Phylum",
       color = "Phylum")  +
  guides(color = guide_legend(byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = .1, hjust = .5),
        legend.spacing.y = unit(.18, "cm"),
        legend.title = element_blank())

ggsave(tax_composition, filename = "./graphs/taxcomp.png", dpi= 300, width = 12, units = "in")

#actinobacteria genus
ps_actino <- subset_taxa(ps, Phylum == "Actinobacteria")
ps_actino_genera <- tax_glom(ps_actino, taxrank="Genus")  # 7 taxonomic ranks

# Get the top 10 genera
top <- names(sort(taxa_sums(ps_actino_genera), decreasing = TRUE))[1:10]

# Calculate relative abundance
ps_actino_genera.prop <- transform_sample_counts(ps_actino_genera, function(x) x / sum(x) )

ps_actino_genera.prop.top <- prune_taxa(top, ps_actino_genera.prop)

otu_prop.top <- data.frame(otu_table(ps_actino_genera.prop.top))
otu_prop.top[is.na(otu_prop.top)] <- 0
ps_actino_genera.prop.top <- 
  phyloseq(tax_table(ps_actino_genera.prop.top), 
           sample_data(ps_actino_genera.prop.top),
           otu_table(otu_prop.top, taxa_are_rows = FALSE))

ps_ARM = merge_samples(ps_actino_genera.prop.top, "Structure")
ps_ARM@sam_data$Structure <- row.names(ps_ARM@sam_data)
ps_ARM = transform_sample_counts(ps_ARM, function(x) 100 * x/sum(x))

tax_composition_genus <- 
  plot_bar(ps_ARM, "Structure", fill = "Genus") + 
  coord_flip() + 
  geom_bar(aes(color = Genus, fill= Genus), stat = "identity", position = "stack") +
  scale_color_manual(values = pal, aesthetics = c("color", "fill")) +
  scale_x_discrete(expand = c(0,0)) +
  ylim(0,100) +
  labs(x = "",
       y = "Relative Abundance",
       fill = "Genus",
       color = "Genus") + 
  theme(plot.title = element_text(face = "bold", color = "#000000", size = 35), 
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(face = 'bold', size = 20, color = "#000000"),
        axis.title = element_text(face = "bold", size = 25), 
        axis.text.x = element_text(angle = .5, hjust = .5),
        axis.text.y = element_text(angle = 90, hjust = .5, vjust = 0),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "#000000"),
        legend.text = element_text(face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 25),
        legend.title.align = 0.5)


ggsave(tax_composition_genus, filename = "./graphs/taxcompGenus.png", dpi= 300, width = 12, units = "in")

# Get relative abundance of actinobacteria for each sample type
tax_table(ps)[,2] %>% table()
2628 / 7058 # 37% of all samples were in the phylum actinobacteria
1910 / 7078 # 27% of all samples were in the phylum proteo bacteria
(1910 + 2628) / 7058 # 64% of all samples were in either acino bacteria or proteobacteri

ps_root <- subset_samples(ps, Structure=="Root")
ps_root <-  prune_taxa(taxa_sums(ps_root) > 0, ps_root)
tax_table(ps_root)[,2] %>% table()
2150 / 5654 # 38% root sequences are in the phylum Actinobacteria

ps_shoot <- subset_samples(ps, Structure=="Shoot")
ps_shoot <-  prune_taxa(taxa_sums(ps_shoot) > 0, ps_shoot)
tax_table(ps_shoot)[,2] %>% table()
555 / 1415 # 39% shoot sequences are in the phylum Actinobacteria

ps_soil <- subset_samples(ps, Structure=="Soil")
ps_soil <-  prune_taxa(taxa_sums(ps_soil) > 0, ps_soil)
tax_table(ps_soil)[,2] %>% table()
628 / 1569 # 40% Soil sequences are in the phylum Actinobacteria

tax_table(ps)[,2] %>% table()
(1910 + 2628) / 7058 # 64% of all sequences are in Actinobacteria or Proteobacteria
2628 / 7058 # 37% in actinobacteria
1910 / 7058 # 27% in proteobacteria


ps_location <- 
  merge_samples(ps, "Location", fun = sum)
sam <- ps_location@sam_data
# repair metadata using sample_names
x <- sample_names(ps_location)

ps_location@sam_data$Location <- 
  x %>% 
  map_chr(1)

ps_location@sam_data$Location <- 
  ps_location@sam_data$Location %>% factor(levels=c("Capitol Reef","Manning Canyon"))

ps_location <- prune_taxa(taxa_sums(ps_location) > 0, ps_location)
ps_location <- prune_samples()

ps_location %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_heatmap(taxa.label = "Phylum",sample.label = "Location") + 
  facet_wrap(~Location,scales = "free_x")
