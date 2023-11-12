# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Looking at abundance of the bacterial communities
# Author: Josh Leon
# Software versions:  R v 4.1.1
#                     tidyverse v 1.3.2
#                     phyloseq v 1.36.0
# -----------------------------------------------------------------------------#

library(phyloseq)
library(tidyverse)

theme_set(theme_bw())

ps <- readRDS("./Final_phyloseq_object.RDS")
ps_ra <- transform_sample_counts(ps, fun = function(x){x/sum(x)})
# Heatmaps ###
# Family ###
ps_family <- tax_glom(ps, "Family")
(p <- plot_heatmap(ps_family, "NMDS", 
                   "bray", 
                   "Structure", 
                   "Family",
                   low="#FFFFCC", 
                   high="#6A0213",
                   na.value="white",
                   sample.order = "Structure"))
p
# Phylum ###
ps_phylum <- tax_glom(ps_ra, "Phylum")
(p2 <- plot_heatmap(ps_phylum, "NMDS", "bray", "Structure", "Phylum",
                   low="#FFFFCC", 
                   high="#6A0213",
                   na.value="white",
                   sample.order = "Structure"))
p2

# Actinobacteria ###
ps_actino <- subset_taxa(ps, Phylum == "Actinobacteria")

top_10_act_genera <- 
  data.frame(ps_actino@tax_table) %>% 
  group_by(Genus) %>% 
  summarize(N = n(),
            Relative_Abundance = (N / 2659) * 100) %>% 
  filter(N >= 51)

top10overall <- top_10_act_genera$Genus

all_other_act_genera <- 
  tax <- data.frame(ps_actino@tax_table) %>% 
  group_by(Genus) %>% 
  summarize(N = n()) %>% 
  filter(N < 51)

ps_root <- subset_samples(ps_actino, Structure=="Root")
ps_root <-  prune_taxa(taxa_sums(ps_root) > 0, ps_root)

ps_shoot <- subset_samples(ps_actino, Structure=="Shoot")
ps_shoot <-  prune_taxa(taxa_sums(ps_shoot) > 0, ps_shoot)

ps_soil <- subset_samples(ps_actino, Structure=="Soil")
ps_soil <-  prune_taxa(taxa_sums(ps_soil) > 0, ps_soil)

top_10R_act_genera <- 
  data.frame(ps_root@tax_table) %>% 
  group_by(Genus) %>% 
  summarize(N = n(),
            Relative_Abundance = (N / 2150) * 100) %>% 
  filter(N >= 46) %>% 
  filter(Genus %in% top10overall)

top_10SH_act_genera <- 
  data.frame(ps_shoot@tax_table) %>% 
  group_by(Genus) %>% 
  summarize(N = n(),
            Relative_Abundance = (N / 555) * 100) %>% 
  filter(N >= 14) %>% 
  filter(Genus %in% top10overall)

top_10SO_act_genera <- 
  data.frame(ps_soil@tax_table) %>% 
  group_by(Genus) %>% 
  summarize(N = n(),
            Relative_Abundance = (N / 628) * 100) %>% 
  filter(N >= 14) %>% 
  filter(Genus %in% top10overall)

Root_abundances <- c(0.00,top_10R_act_genera$Relative_Abundance)
Shoot_abundances <- top_10SH_act_genera$Relative_Abundance
Shoot_abundances <- append(Shoot_abundances, 0.00, after = 4)
Shoot_abundances <- append(Shoot_abundances, 0.00, after = 9)
Soil_abundances <- top_10SO_act_genera$Relative_Abundance
Soil_abundances <- append(Soil_abundances, 0.00, after = 6)

top_10_genera <- top_10_act_genera %>% 
  mutate(Root = Root_abundances,
         Shoot = Shoot_abundances,
         Soil = Soil_abundances)

top_10_genera$N <- NULL
top_10_genera$Relative_Abundance <- NULL
Abundance_table <- top_10_genera %>% 
  pivot_longer(cols = c("Root","Shoot", "Soil"),names_to = "Sample_Type",values_to = "Percentage") %>% 
  filter(!is.na(Genus))

matrix_abundance_table <- top_10_genera %>% filter(!is.na(Genus))
names <- matrix_abundance_table$Genus
matrix_abundance_table$Genus <- NULL
row.names(matrix_abundance_table) <- names

heatmap(as.matrix(matrix_abundance_table), scale = "column", Colv = NA, Rowv = NA)

Abundance_table[Abundance_table == 0] <- NA

Abundance_table <- Abundance_table %>% 
  mutate_if(is.numeric, round, digits = 2) %>% 
  mutate(Relative_Abundance = Percentage / 100)

HeatmapGenus <- 
  Abundance_table %>% 
  ggplot(aes(x = Sample_Type, y = Genus, fill = Relative_Abundance)) +
  geom_tile() +
  scale_fill_gradient(low="#FFFFCC", high="#6A0213", na.value = "#FFFFFF") +
  labs(fill = "Relative\nAbundance") +
  theme(panel.border = element_rect(color = "#000000"),
        axis.text = element_text(face = 'bold', size = 20, color = "#000000"),
        axis.title = element_text(face = "bold", size = 25), 
        axis.text.x = element_text(angle = .5, hjust = .5),
        axis.text.y = element_text(angle = .1, hjust = .5, vjust = 0),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 25),
        legend.title.align = 0.5,
        legend.key.width=unit(2,"cm"),
        legend.key.height = unit(2, "cm"),
        legend.key = element_rect(color = "#000000"))

ggsave(HeatmapGenus, filename = "./graphs/HeatmapGenus.png", dpi= 300, width = 12, units = "in")

# 10 Most abundant Species ###
# Never mind 98% of all ASVs could not be identified to the species level

# Network
set.seed(123)

plot_net(ps, maxdist = .9, point_label = "Structure")
plot_net(ps, maxdist = .9, color = "Structure", shape="Location")
plot_net(ps, maxdist = .9, color = "Structure")

pal.network = c("#6A0213", "#359B73", "#FF7F52")

ig <- make_network(ps, max.dist=.95, dist.fun="bray")
net1 <- 
  plot_network(ig, ps, 
             color = "Structure",
             shape = "Location", 
             label = NULL, 
             line_weight = 0.4) +
  geom_point(aes(color = Structure)) +
  scale_color_manual(values = pal.network, aesthetics = "color")

ggsave(net1, filename ="./graphs/SampleLocation_network.png", dpi = 300)

net2 <- 
  plot_network(ig, ps, 
               color = "Structure",
               label = NULL, 
               line_weight = 0.4) +
  geom_point(aes(color = Structure)) +
  scale_color_manual(values = pal.network, aesthetics = "color")

ggsave(net2, filename ="./graphs/Sample_network.png", dpi = 300)


# Tree   
head(phy_tree(ps)$node.label, 10)
ntaxa(ps)

fultree <- 
  plot_tree(ps,color="Phylum", base.spacing=0.03, ladderize = TRUE) +
  scale_color_manual(values = pal2, aesthetics = "color")

top50otus = names(sort(taxa_sums(ps), TRUE)[1:50])

ps50 <- prune_taxa(top50otus, ps)
tree50 <- 
  plot_tree(ps50,color="Structure",shape = "Phylum", base.spacing=0.03, ladderize = TRUE, label.tips = "Genus") +
  scale_color_manual(values = pal.network, aesthetics = "color")

ggsave(fultree, filename = "./graphs/fulltree.png", dpi = 300)
ggsave(tree50, filename = "./graphs/top50ASVs_tree.png", dpi = 300)

