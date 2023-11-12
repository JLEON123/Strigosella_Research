#-----------------------------------------------#
### Strigosella africana Minion 16S           ###
### Exploring the Cleaned Data                ###
### Date: 02/10/2023                          ###
### Author: Josh Leon                         ###
### Software Versions: R v 4.1.1              ###
###                    tidyverse v 1.3.2      ###
#-----------------------------------------------#

# IMPORTS ####
library(tidyverse)
library(vegan)
library(reshape2)
library(betapart)

source("./scripts/theme.R")
source("./scripts/random_discretecolors.R")


# SETUP ####
df <- read_rds("./data/cleaned_strigosella")
df$Structure <- as.factor(df$Structure)

# phylum stuff
phylum_df <- 
  df %>% 
  group_by(Structure, phylum) %>% 
  summarize(N = n())

phylum_df2 <- 
  phylum_df %>% 
  group_by(Structure) %>% 
  summarize(Totals = sum(N))

phylum_full <- 
  merge(phylum_df, phylum_df2)

phylum_full <- 
  phylum_full %>% 
  mutate(Relative_Abundance = N / Totals)

# Genus stuff
genus_df <- 
  df %>% 
  group_by(Structure, genus) %>% 
  summarize(N = n())

genus_df2 <- 
  genus_df %>% 
  group_by(Structure) %>% 
  summarize(Totals = sum(N))

genus_full <- 
  merge(genus_df, genus_df2)

genus_full <- 
  genus_full %>% 
  mutate(Relative_Abundance = N / Totals)

genus_full_root <- 
  genus_full %>% 
  filter(Structure == "Root")

genus_full_shoot <- 
  genus_full %>% 
  filter(Structure == "Shoot")

genus_full_soil <- 
  genus_full %>% 
  filter(Structure == "Soil")

shannonstuff <- dcast(df,barcode~genus)
shannonstuff <- 
  shannonstuff %>% select(!barcode)
  
H <- diversity(shannonstuff)
# Observed Richness
richness <- specnumber(shannonstuff)  

# Pielou's Evenness
evenness <- H/log(shannonstuff)

alpha <- cbind(shannon = H, richness = richness, evenness = evenness)
head(alpha)

alpha$Structure <- c("Shoot","Shoot", "Shoot","Shoot", "Shoot","Root","Root","Root","Root","Root","Soil","NC")

alpha <- 
  alpha %>% 
  filter(Structure != "NC")

plot.shan <- ggplot(alpha, aes(x = Structure, y = shannon, colour = Structure)) +
  geom_point(size = 5, show.legend = FALSE) +
  scale_colour_viridis_d(option = "E", begin = 0.2, end = 0.8) +
  ylab("Shannon Diversity") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.text.y = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.title.y = element_text(hjust = .5, vjust = 0.4, face = 'bold', size = 15))

plot.rich <-ggplot(alpha, aes(x = Structure, y = richness, colour = Structure)) +
  geom_point(size = 5, show.legend = FALSE) +
  scale_colour_viridis_d(option = "E", begin = 0.2, end = 0.8) +
  ylab("Genera Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.text.y = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.title.y = element_text(hjust = .5, vjust = 0.4, face = 'bold', size = 15))

ggsave("./figures/shannondiv.png", plot = plot.shan, dpi = 300)
ggsave("./figures/richness.png", plot = plot.rich, dpi = 300)

# Species stuff
species_df <- 
  df %>% 
  group_by(Structure, species) %>% 
  summarize(N = n())

species_df2 <- 
  species_df %>% 
  group_by(Structure) %>% 
  summarize(Totals = sum(N))

species_full <- 
  merge(species_df, species_df2)

species_full <- 
  species_full %>% 
  mutate(Relative_Abundance = N / Totals)

shannonstuff2 <- dcast(df,barcode~species)
shannonstuff2 <- 
  shannonstuff2 %>% select(!barcode)

H2 <- diversity(shannonstuff2)
# Observed Richness
richness2 <- specnumber(shannonstuff2)  

# Pielou's Evenness
evenness2 <- H2/log(shannonstuff2)

alpha2 <- cbind(shannon = H2, richness = richness2, evenness = evenness2)
head(alpha2)

alpha2$Structure <- c("Shoot","Shoot", "Shoot","Shoot", "Shoot","Root","Root","Root","Root","Root","Soil","NC")

alpha2 <- 
  alpha2 %>% 
  filter(Structure != "NC")

plot.shan2 <- ggplot(alpha2, aes(x = Structure, y = shannon, colour = Structure)) +
  geom_point(size = 5, show.legend = FALSE) +
  scale_colour_viridis_d(option = "E", begin = 0.2, end = 0.8) +
  ylab("Shannon Diversity") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.text.y = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.title.y = element_text(hjust = .5, vjust = 0.4, face = 'bold', size = 15))

plot.rich2 <-ggplot(alpha2, aes(x = Structure, y = richness, colour = Structure)) +
  geom_point(size = 5, show.legend = FALSE) +
  scale_colour_viridis_d(option = "E", begin = 0.2, end = 0.8) +
  ylab("Species Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.text.y = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13),
        axis.title.y = element_text(hjust = .5, vjust = 0.4, face = 'bold', size = 15))

ggsave("./figures/shannondiv2.png", plot = plot.shan2, dpi = 300)
ggsave("./figures/richness2.png", plot = plot.rich2, dpi = 300)
# Graphs ####
phylum_full$Structure <- 
  phylum_full$Structure %>% factor(levels = c("Soil", "Root", "Shoot", "Negative Control"))

Phyla_abundance <- 
  phylum_full %>% 
  filter(Structure != "Negative Control",
         !is.na(phylum)) %>% 
  ggplot(aes(x = Structure, y = Relative_Abundance, fill = phylum)) +
  geom_col() +
  scale_fill_manual(values = col_vector[56:74]) +
  labs(x = "",
       y = "",
       fill = "Phyla") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13,
                                   color = "black"),
        axis.text.y = element_text(hjust = 1, vjust = 0.4, face = 'bold', size = 13, color = "black"),
        panel.grid = element_blank())

ggsave(filename = "./figures/phyla_abundance.png",plot = Phyla_abundance, dpi = 300)




# Beta Diversity ####
beta <- 
  alpha2 %>% 
  select(!contains("Structure"))

bray <- bray.part(x = beta)

bd<-betadisper(bray[[3]], group)

plot(bd)

dist<-beta.pair(presabs, index.family="jaccard")
NMDS <- 
  alpha2 %>% transform
  
  
  NMDS <- 
  ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS", distance = "bray")




